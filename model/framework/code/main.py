# imports
import os
import sys
import csv
import json
import yaml
import pickle
import pandas as pd

# current file directory (needed before patching)
root = os.path.dirname(os.path.abspath(__file__))
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))

# Patch chebifier to load chebi_graph.pkl from local checkpoints instead of HuggingFace.
# Must happen before cli_adapted is imported (which triggers base_ensemble import).
import chebifier.ensemble.base_ensemble as _base_ensemble
_orig_load_chebi_graph = _base_ensemble.load_chebi_graph
def _local_load_chebi_graph(filename=None):
    local = os.path.join(checkpoints_dir, "chebi_graph.pkl")
    if filename is None and os.path.isfile(local):
        print("Loading ChEBI graph from local checkpoints...")
        return pickle.load(open(local, "rb"))
    return _orig_load_chebi_graph(filename)
_base_ensemble.load_chebi_graph = _local_load_chebi_graph

_orig_get_disjoint_files = _base_ensemble.get_disjoint_files
def _local_get_disjoint_files():
    local = [
        os.path.join(checkpoints_dir, "disjoint_chebi.csv"),
        os.path.join(checkpoints_dir, "disjoint_additional.csv"),
    ]
    if all(os.path.isfile(f) for f in local):
        return local
    return _orig_get_disjoint_files()
_base_ensemble.get_disjoint_files = _local_get_disjoint_files

# Patch ChEBILookupPredictor.get_smiles_lookup to use an absolute path in checkpoints
# instead of a path relative to CWD (which breaks when CWD is /tmp).
from chebifier.prediction_models.chebi_lookup import ChEBILookupPredictor as _ChEBILookupPredictor

def _local_get_smiles_lookup(self_inner):
    local = os.path.join(checkpoints_dir, "smiles_lookup.json")
    if os.path.isfile(local):
        print("Loading SMILES lookup from local checkpoints...")
        with open(local, "r", encoding="utf-8") as f:
            return json.load(f)
    print("Building SMILES lookup (first run, may take a few minutes)...")
    smiles_lookup = self_inner.build_smiles_lookup()
    with open(local, "w", encoding="utf-8") as f:
        json.dump(smiles_lookup, f, indent=4)
    return smiles_lookup

_ChEBILookupPredictor.get_smiles_lookup = _local_get_smiles_lookup

from criteria import *
from cli_adapted import predict


mapping_features_to_functions = {
    "aminocoumarins": aminocoumarins,
    "aminoglycosides_aminocyclitols": aminoglycosides_aminocyclitols,
    "ansamycins_rifamycins_macrolides": ansamycins_rifamycins_macrolides,
    "b_lactams_carbapenems": b_lactams_carbapenems,
    "b_lactams_cephalosporins_cephems": b_lactams_cephalosporins_cephems,
    "b_lactams_penicillins": b_lactams_penicillins,
    "b_lactams_monobactams": b_lactams_monobactams,
    "b_lactamase_inhibitors": b_lactamase_inhibitors,
    "b_lactams_all": b_lactams_all,
    "aminopyrimidines_trimethoprim_der": aminopyrimidines_trimethoprim_der,
    "glycopeptides": glycopeptides,
    "lipopeptides": lipopeptides,
    "nitrofurans": nitrofurans,
    "nucleosides": nucleosides,
    "oxazolidinones": oxazolidinones,
    "phenols": phenols,
    "polyketides_type_I": polyketides_type_I,
    "polyketides_type_II": polyketides_type_II,
    "polypeptides": polypeptides,
    "fluoroquinolones": fluoroquinolones,
    "quinolines": quinolines,
    "quat_ammonium_cpds": quat_ammonium_cpds,
    "sulfonamides": sulfonamides,
    "arsenic_cpds": arsenic_cpds,
    "amino_acid_der": amino_acid_der,
    "phosphonic_acid_der": phosphonic_acid_der,
    "hydrazides": hydrazides,
    "mc_fatty_acids": mc_fatty_acids,
    "lc_fatty_acids": lc_fatty_acids,
    "pyrazines": pyrazines,
    "biguanides": biguanides,
    "depsipeptides": depsipeptides,
    "peroxides": peroxides,
    "pyridinium": pyridinium,
    "antifungal": antifungal,
    "quinolones": quinolones,
    "penams": penams,
    "anthracyclines": anthracyclines,
}

# parse arguments
input_file = os.path.abspath(sys.argv[1])
output_file = os.path.abspath(sys.argv[2])
tmp_file = os.path.abspath(output_file.replace(".csv", '_tmp.csv'))
output_file_chebified = os.path.abspath(output_file.replace(".csv", "_chebified.json"))

# read smiles and create tmp file
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

with open(tmp_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["smiles"])
    for s in smiles_list:
        writer.writerow([s])

# generate ensemble config pointing to local checkpoint files (avoids HuggingFace downloads)
local_ensemble_config = {
    "electra": {
        "type": "electra",
        "ckpt_path": os.path.join(checkpoints_dir, "14ko0zcf_epoch=193.ckpt"),
        "target_labels_path": os.path.join(checkpoints_dir, "electra_classes.txt"),
        "classwise_weights_path": os.path.join(checkpoints_dir, "metrics_electra_14ko0zcf_80-10-10_short.json"),
    },
    "resgated": {
        "type": "resgated",
        "ckpt_path": os.path.join(checkpoints_dir, "0ps1g189_epoch=122.ckpt"),
        "target_labels_path": os.path.join(checkpoints_dir, "resgated_classes.txt"),
        "classwise_weights_path": os.path.join(checkpoints_dir, "metrics_0ps1g189_80-10-10_short.json"),
        "molecular_properties": [
            "chebai_graph.preprocessing.properties.AtomType",
            "chebai_graph.preprocessing.properties.NumAtomBonds",
            "chebai_graph.preprocessing.properties.AtomCharge",
            "chebai_graph.preprocessing.properties.AtomAromaticity",
            "chebai_graph.preprocessing.properties.AtomHybridization",
            "chebai_graph.preprocessing.properties.AtomNumHs",
            "chebai_graph.preprocessing.properties.BondType",
            "chebai_graph.preprocessing.properties.BondInRing",
            "chebai_graph.preprocessing.properties.BondAromaticity",
            "chebai_graph.preprocessing.properties.RDKit2DNormalized",
        ],
    },
    "chemlog_peptides": {"type": "chemlog_peptides"},
    "chemlog_element": {"type": "chemlog_element"},
    "chemlog_organox": {"type": "chemlog_organox"},
    "c3p": {
        "type": "c3p",
        "classwise_weights_path": os.path.join(checkpoints_dir, "c3p_trust.json"),
    },
    "chebi_lookup": {"type": "chebi_lookup"},
}
local_ensemble_config_path = "/tmp/ensemble_config_local.yml"
with open(local_ensemble_config_path, "w") as f:
    yaml.dump(local_ensemble_config, f)

# change working directory to a writable location before running the model
# (chemlog_extra writes relative path data/ from cwd at init time,
#  which fails if cwd is inside a read-only container image)
os.chdir("/tmp")

# run the model
predict(
    ensemble_config=local_ensemble_config_path,
    smiles=(),  # none inline
    smiles_file=tmp_file,
    output=output_file_chebified,
    ensemble_type="wmv-f1",
    use_confidence=True,
    resolve_inconsistencies=True
)


# read columns file
columns = pd.read_csv(os.path.join(root, "..", "fit", 'data', "Final_column_criteria.csv"))
features_to_classes = {i: [j, k] for i,j,k in zip(columns['Column name'], columns["Class"], columns['CHEBI ids'])}
features = sorted(features_to_classes)

# read chebifier output and map to antibiotic classes
with open(output_file_chebified) as f:
    chebifier_output = json.load(f)

output_content = []
for smi in smiles_list:
    chebi_predictions = chebifier_output.get(smi)
    if chebi_predictions is None:
        fp = [None] * len(features)
    else:
        chebi_ids = set("CHEBI:" + o for o in chebi_predictions)
        fp = [mapping_features_to_functions[i](chebi_ids) for i in features]
    output_content.append(fp)

output_content = pd.DataFrame(output_content, columns=features)
output_content.to_csv(output_file, sep=',', index=False)

# remove tmp files
os.remove(output_file_chebified)
os.remove(tmp_file)

