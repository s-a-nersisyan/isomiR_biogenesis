import pandas as pd
import numpy as np
import sys
import os

PATH_TO_MAIN_ISOMIRS = "../../processed_data/TCGA/main_isomiRs.tsv"
OUTPUT_DIR = "../../processed_data/TCGA"

# fraction of TCGA projects that should contain a particular
# main isoform to make it main throughout all projects
project_fraction = 0.5
isoform_treshold = 0.95

frequency_tb = pd.read_csv(PATH_TO_MAIN_ISOMIRS, sep="\t")
miRNAs = set(frequency_tb["miRNA"])

# TCGA project merged data
merged_tb = pd.DataFrame(columns=[
    "miRNA",
    "shift_5",
    "heterogenous_5",
    "shift_3",
    "heterogenous_3"
])

for miRNA in miRNAs:
    # extract all TCGA projects containing isoforms of a particular miRNA
    miRNA_tb = frequency_tb[frequency_tb["miRNA"] == miRNA]

    # extrac TCGA maximum spreaded shifts of the miRNA - the main consensus isomiR
    shifts_5, counts_5 = np.unique(list(miRNA_tb["shift_5"]),
                                   return_counts=True)
    shifts_3, counts_3 = np.unique(list(miRNA_tb["shift_3"]),
                                   return_counts=True)
    max_shift_5 = shifts_5[np.argmax(counts_5)]
    max_count_5 = np.max(counts_5)
    max_shift_3 = shifts_3[np.argmax(counts_3)]
    max_count_3 = np.max(counts_3)


    # count of TCGA projects having the main consensus isomiR as the main
    isof_count_5 = len(miRNA_tb[
        (miRNA_tb["shift_5"] == max_shift_5) &
        (miRNA_tb["fraction_5"] <= isoform_treshold)
    ])

    # count of TCGA projects having the main consensus isomiR as the main
    isof_count_3 = len(miRNA_tb[
        (miRNA_tb["shift_3"] == max_shift_3) &
        (miRNA_tb["fraction_3"] <= isoform_treshold)
    ])

    if (isof_count_5 / max_count_5 < project_fraction):
        heter_5 = 0
    else:
        heter_5 = 1

    if (isof_count_3 / max_count_3 < project_fraction):
        heter_3 = 0
    else:
        heter_3 = 1

    merged_tb.loc[len(merged_tb)] = [
        miRNA,
        max_shift_5,
        heter_5,
        max_shift_3,
        heter_3
    ]

merged_tb.to_csv(os.path.join(
    OUTPUT_DIR, "main_consensus_isomiRs.tsv"), sep="\t", index=None)
