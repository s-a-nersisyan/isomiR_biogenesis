import pandas as pd
import numpy as np
import tqdm
import sys
import os

TCGA_PATH = "/home/dude/huge/bulk/TCGA/"
TOP = 0.2

OUT_PATH = "../../data/"

organs = [organ for organ in os.listdir(TCGA_PATH)
          if os.path.isdir(os.path.join(TCGA_PATH, organ))]
organs = [organ for organ in organs
          if os.path.exists(os.path.join(TCGA_PATH, organ, "miRNA", "isomiR_CPM.tsv"))]

result = {
    "TCGA":[],
    "miRNA":[],
    "5frac": [],
    "3frac": [],
    "5shift": [],
    "3shift": [],
}

for organ in tqdm.tqdm(organs):
    # tumor patients indexes
    tumor_idx = list(pd.read_csv(os.path.join(
        TCGA_PATH, organ, "miRNA", "tumor_samples.tsv"), sep="\t")["Tumor ID"])

    # extract isomiRNA expressions of tumor patients in TCGA-"organ" projects
    isomiR_CPM = pd.read_csv(os.path.join(
        TCGA_PATH, organ, "miRNA", "isomiR_CPM.tsv"), sep="\t", index_col=0)[tumor_idx]

    # normalize expressions
    isomiR_CPM = 2 ** isomiR_CPM - 1

    # get TOP expressioned isomiRNAs
    quantile = isomiR_CPM.median(axis=1).quantile(q=1 - TOP)
    miRNAs = list(isomiR_CPM[isomiR_CPM.median(axis=1) >= quantile].index)

    # get TOP expressioned miRNAs cutting isomiRNAs shifts
    miRNAs = set([index.split("|")[0] for index in miRNAs])

    # get isoforms of TOP expressioed miRNAs
    new_indexes = []
    for index in isomiR_CPM.index:
        if index.split("|")[0] in miRNAs:
            new_indexes.append(index)
    isomiR_CPM = isomiR_CPM.loc[new_indexes]

    # get canonicial isoform and its fraction for each TOP miRNA
    for miRNA in miRNAs:
        # isoforms of the particular miRNA
        isomiRs = [isomiR for isomiR in isomiR_CPM.index if isomiR.startswith(miRNA)]

        shifts_5 = set()
        shifts_3 = set()
        for isomiR in isomiRs:
            shift_5 = isomiR.split("|")[1]
            shift_3 = isomiR.split("|")[2]

            shifts_5.add(shift_5)
            shifts_3.add(shift_3)

        shifts_5 = sorted(list(shifts_5))
        shifts_3 = sorted(list(shifts_3))

        # compute medians (by tumor patients) of 5 shifted miRNA isoform, summing
        # this variant by all available 3 shifts
        medians_5 = []
        for shift in shifts_5:
            median_val = isomiR_CPM.loc[isomiR_CPM.index.str.startswith(
                miRNA + "|" + shift)].sum(axis=0).median()
            medians_5.append(median_val)

        # compute medians (by tumor patients) of 5 (5prime) shifted miRNA isoform, summing
        # this variant by all available 3 (3prime) shifts
        medians_3 = []
        for shift in shifts_3:
            median_val = isomiR_CPM.loc[(isomiR_CPM.index.str.startswith(
                miRNA) & isomiR_CPM.index.str.endswith(shift))].sum(axis=0).median()
            medians_3.append(median_val)

        # compute canonical 5prime and 3prime isoforms
        miRNA_max_shift_5 = shifts_5[np.argmax(medians_5)]
        miRNA_max_shift_3 = shifts_3[np.argmax(medians_3)]

        # compute fracton of 5prime and 3prime canonincal forms
        miRNA_exp = isomiR_CPM.loc[(isomiR_CPM.index.str.startswith(miRNA))].sum(axis=0)
        miRNA_exp5 = isomiR_CPM.loc[(isomiR_CPM.index.str.startswith(
            miRNA + "|" + miRNA_max_shift_5))].sum(axis=0)
        miRNA_exp3 = isomiR_CPM.loc[(isomiR_CPM.index.str.startswith(
            miRNA) & isomiR_CPM.index.str.endswith(miRNA_max_shift_3))].sum(axis=0)

        fraction_5 = (miRNA_exp5 / miRNA_exp).median()
        fraction_3 = (miRNA_exp3 / miRNA_exp).median()

        # store information
        result["TCGA"].append(organ)
        result["miRNA"].append(miRNA)

        result["5shift"].append(miRNA_max_shift_5)
        result["3shift"].append(miRNA_max_shift_3)

        result["5frac"].append(fraction_5)
        result["3frac"].append(fraction_3)

data = pd.DataFrame(result, columns=["TCGA", "miRNA", "5shift", "3shift", "5frac", "3frac"])
data = data.sort_values(by=["TCGA", "miRNA", "5shift", "3shift", "5frac", "3frac"])
data.to_csv(os.path.join(OUT_PATH, "canonical_isomiR.tsv"), sep="\t")
