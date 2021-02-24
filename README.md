# isomiR_biogenesis
## Environment
To run the code the following packages should be installed on your machine: NumPy, Pandas, skikit-learn, SciPy.
## How to prepare miRBase and TCGA data
Run the following commands:
    cd scripts/prepare
    python3 parse_miRBase.py
    python3 extract_canonical_isomiRNA.py
    python3 agregate_fraction_table_by_TCGA.py
