# Antibiotic classes prediction

An automated classification of chemicals in the ChEBI ontology based on a neuro-symbolic AI technique that harnesses the ontology itself to create the learning system and enables the classification of coumpounds into GARDP antibiotic classes.

This model was incorporated on 2025-08-26.Last packaged on 2025-08-29.

## Information
### Identifiers
- **Ersilia Identifier:** `eos19mt`
- **Slug:** `chebifier-antibiotic`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Activity prediction`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Antimicrobial activity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1`
- **Output Consistency:** `Fixed`
- **Interpretation:** Presence (1) or absence (0) of ChEBI predicted parents associated with pre-defined GARDP-inspired antibiotic classes

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| amino_acid_der | integer | high | Presence of ChEBI predicted parents associated with Amino acid derivatives |
| aminocoumarins | integer | high | Presence of ChEBI predicted parents associated with Aminocoumarins |
| aminoglycosides_aminocyclitols | integer | high | Presence of ChEBI predicted parents associated with Aminoglycosides (incl. aminocyclitols) |
| aminopyrimidines_trimethoprim_der | integer | high | Presence of ChEBI predicted parents associated with Aminopyrimidines / trimethoprim derivatives |
| ansamycins_rifamycins_macrolides | integer | high | Presence of ChEBI predicted parents associated with Ansamycins (incl. Rifamycins) (macrolides) |
| anthracyclines | integer | high | Presence of ChEBI predicted parents associated with Anthracyclines |
| antifungal | integer | high | Presence of ChEBI predicted parents associated with Antifungals |
| arsenic_cpds | integer | high | Presence of ChEBI predicted parents associated with Arsenic compounds |
| b_lactamase_inhibitors | integer | high | Presence of ChEBI predicted parents associated with β-lactamase inhibitors |
| b_lactams_all | integer | high | Presence of ChEBI predicted parents associated with β-lactams (all) |

_10 of 40 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `Internal`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos19mt](https://hub.docker.com/r/ersiliaos/eos19mt)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos19mt.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos19mt.zip)

### Resource Consumption
- **Model Size (Mb):** `1735`
- **Environment Size (Mb):** `7657`
- **Image Size (Mb):** `8558.11`

**Computational Performance (seconds):**
- 10 inputs: `61.34`
- 100 inputs: `101.93`
- 10000 inputs: `-1`

### References
- **Source Code**: [https://github.com/ChEB-AI/python-chebifier](https://github.com/ChEB-AI/python-chebifier)
- **Publication**: [https://pubs.rsc.org/en/content/articlehtml/2024/dd/d3dd00238a](https://pubs.rsc.org/en/content/articlehtml/2024/dd/d3dd00238a)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2025`
- **Ersilia Contributor:** [arnaucoma24](https://github.com/arnaucoma24)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [GPL-3.0-only](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos19mt
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos19mt
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
