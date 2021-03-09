import pandas as pd
import numpy as np
import sys
import os

PATH_TO_MAIN_ISOMIRS = "../../processed_data/TCGA/main_isomiRs.tsv"
OUTPUT_DIR = "../../processed_data/TCGA"

# fraction of TCGA projects that should contain a particular
# main isoform to make it main throughout all projects
projects_fraction = 0.5
# miRNA has heterogeneous cleavage pattern, if
# main isoform fraction is less then this threshold
hetero_cleavage_threshold = 0.95

df_all = pd.read_csv(PATH_TO_MAIN_ISOMIRS, sep="\t")
miRNAs = set(df_all["miRNA"])

# TCGA project merged data
df_merged = pd.DataFrame(columns=[
    "miRNA",
    "shift_5",
    "heterogeneous_5",
    "shift_3",
    "heterogeneous_3"
])

for miRNA in miRNAs:
    # extract all TCGA projects containing isoforms of a particular miRNA
    miRNA_tb = df_all[df_all["miRNA"] == miRNA]

    # Extract most frequent 5'- and 3'-shifts of miRNA across TCGA projects
    shifts_5, counts_5 = np.unique(miRNA_tb["shift_5"], return_counts=True)
    shifts_3, counts_3 = np.unique(miRNA_tb["shift_3"], return_counts=True)

    max_shift_5 = shifts_5[np.argmax(counts_5)]
    max_count_5 = np.max(counts_5)
    
    max_shift_3 = shifts_3[np.argmax(counts_3)]
    max_count_3 = np.max(counts_3)

    # Number of TCGA projects with the most frequent 5'-shift of
    # miRNA and heterogeneous cleavage pattern
    isof_count_5 = len(miRNA_tb[
        (miRNA_tb["shift_5"] == max_shift_5) &
        (miRNA_tb["fraction_5"] <= hetero_cleavage_threshold)
    ])

    # Number of TCGA projects with the most frequent 3'-shift of
    # miRNA and heterogeneous cleavage pattern
    isof_count_3 = len(miRNA_tb[
        (miRNA_tb["shift_3"] == max_shift_3) &
        (miRNA_tb["fraction_3"] <= hetero_cleavage_threshold)
    ])

    if (isof_count_5 / max_count_5 < projects_fraction):
        hetero_5 = 0
    else:
        hetero_5 = 1

    if (isof_count_3 / max_count_3 < projects_fraction):
        hetero_3 = 0
    else:
        hetero_3 = 1

    df_merged.loc[len(df_merged)] = [
        miRNA,
        max_shift_5,
        hetero_5,
        max_shift_3,
        hetero_3
    ]

df_merged = df_merged.sort_values("miRNA")
df_merged.to_csv(
    os.path.join(OUTPUT_DIR, "main_consensus_isomiRs.tsv"), 
    sep="\t", index=None
)
