# imports
import os
import sys
import csv
import json
import pandas as pd
from criteria import *


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
    "benzenesulfonyls": benzenesulfonyls,
    "peroxides": peroxides,
    "pyridinium": pyridinium,
    "antifungal": antifungal,
    # "heterocyclic_antibiotic": heterocyclic_antibiotic,
    "quinolones": quinolones,
    "penams": penams,
    "anthracyclines": anthracyclines,
}

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# import chebifier
sys.path.insert(0, os.path.join(root, "python-chebifier"))
from chebifier.cli_adapted import predict

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
tmp_file = output_file.replace(".csv", '_tmp.csv')
output_file_chebified = output_file.replace(".csv", "_chebified.json")

# read smiles and create tmp file
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

with open(tmp_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["input"])
    for s in smiles_list:
        writer.writerow([s])

# change working directory before running the model
os.chdir(os.path.join(root, "python-chebifier"))

# run the model
predict(
    ensemble_config=os.path.join(root, "..", "..", "checkpoints", "ensemble_config.yml"),
    smiles=(),  # none inline
    smiles_file=os.path.join(root, "..", tmp_file),
    output=os.path.join(root, "..", output_file_chebified),
    ensemble_type="wmv-f1",
    chebi_version=241,
    use_confidence=True,
    resolve_inconsistencies=True
)


# read columns file
columns = pd.read_csv(os.path.join(root, "..", "fit", 'data', "Final_column_criteria.csv"))
features_to_classes = {i: [j, k] for i,j,k in zip(columns['Column name'], columns["Class"], columns['CHEBI ids'])}
features = sorted(features_to_classes)

# for feature in features:
#     print(",".join([feature, "integer", "high", "Presence of ChEBI predicted parents associated with " + 
#                     str(features_to_classes[feature][0])]))  # + " - '" + str(features_to_classes[feature][1]) + "'"]))

# read json output
output = json.load(open(os.path.join(root, "..", output_file_chebified)))
output_content = []
for smi in smiles_list:

    # chebifier
    r = ["CHEBI:" + o for o in sorted(output[smi])]
    r = set(r)

    # eos19mt
    eos19mt_output = []
    fp = [mapping_features_to_functions[i](r) for i in features]
    output_content.append(fp)

output_content = pd.DataFrame(output_content, columns=features)
output_content.to_csv(os.path.join(root, "..", output_file), sep=',', index=False)

# remove json output
os.remove(os.path.join(root, "..", output_file_chebified))
os.remove(os.path.join(root, "..", tmp_file))

