import numpy as np
import pandas as pd
import operator
import sys
import os

# fraction of TCGA projects that should contain a particular
# canonical isoform to make it canon. through all projects
FRACTION = 0.5
ISOFORM_TRESHOLD = 0.95

INPUT_DIR = "../../data/"
PATH_TO_FREQUENCY_TABLE = os.path.join(INPUT_DIR, "canonical_isomiR.tsv")

frequency_table = pd.read_csv(PATH_TO_FREQUENCY_TABLE, sep="\t", index_col=0)
miRNAs = set(frequency_table["miRNA"])

# TCGA project agregated data
agreg_data = {
    "miRNA" : [],
    "5shift" : [],
    "5isof" : [],
    "3shift" : [],
    "3isof" : []
}

isof_count5 = 0
isof_count3 = 0
nonisof_count5 = 0
nonisof_count3 = 0

for miRNA in miRNAs:
    # extract all TCGA projects containing isoforms of a particular miRNA
    miRNA_table = frequency_table[frequency_table["miRNA"] == miRNA]

    # get counting dicts (5prime and 3prime) of canonical isoforms
    shift5_freq_dict = dict()
    shift3_freq_dict = dict()
    for index, row in miRNA_table.iterrows():
        shift5, frac5 = row["5shift"], row["5frac"]
        shift3, frac3 = row["3shift"], row["3frac"]
        if (shift5 in shift5_freq_dict):
            shift5_freq_dict[shift5] += 1
        else:
            shift5_freq_dict[shift5] = 1

        if (shift3 in shift3_freq_dict):
            shift3_freq_dict[shift3] += 1
        else:
            shift3_freq_dict[shift3] = 1

    # get aggregated canonicial (ac) isoforms
    max_shift5 = max(shift5_freq_dict.items(), key=operator.itemgetter(1))[0]
    max_shift3 = max(shift3_freq_dict.items(), key=operator.itemgetter(1))[0]

    # count of ac isoform in TCGA projects
    count5 = miRNA_table[miRNA_table["5shift"] == max_shift5].shape[0]
    # count of ac isoform in TCGA projects,
    # s.t. its frequncy is greater than ISOFORM_TRESHOLD
    count5_freq = miRNA_table[(miRNA_table["5shift"] == max_shift5) &
                              (miRNA_table["5frac"] >= ISOFORM_TRESHOLD)].shape[0]

    # the same thing for 3prime canonical isoform
    count3 = miRNA_table[miRNA_table["3shift"] == max_shift3].shape[0]
    count3_freq = miRNA_table[(miRNA_table["3shift"] == max_shift3) &
                              (miRNA_table["3frac"] >= ISOFORM_TRESHOLD)].shape[0]

    if (count5_freq / count5 < FRACTION):
        isof_count5 += 1
        agreg_data["5isof"].append(1)
    else:
        nonisof_count5 += 1
        agreg_data["5isof"].append(0)
    if (count3_freq / count3 < FRACTION):
        isof_count3 += 1
        agreg_data["3isof"].append(1)
    else:
        nonisof_count3 += 1
        agreg_data["3isof"].append(0)

    agreg_data["miRNA"].append(miRNA)
    agreg_data["5shift"].append(max_shift5)
    agreg_data["3shift"].append(max_shift3)

table = pd.DataFrame(agreg_data, columns=["miRNA", "5shift", "5isof", "3shift", "3isof"])
table.to_csv(os.path.join(INPUT_DIR, "isoform_presence.tsv"), sep="\t")
print("5 isof: {}, nonisof: {}".format(isof_count5, nonisof_count5))
print("3 isof: {}, nonisof: {}".format(isof_count3, nonisof_count3))
