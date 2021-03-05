import pandas as pd
import numpy as np
import tqdm
import sys
import os

TCGA_path = "../../raw_data/TCGA"
# This quantile will be used to determine
# whether isomiR is "expressed" or not
expression_quantile = 0.8

TCGA_projects = [filename[7:][:-4] for filename in os.listdir(TCGA_path)]

res = pd.DataFrame(columns=[
    "Project", "miRNA",
    "shift_5", "median_expr_5", "fraction_5",
    "shift_3", "median_expr_3", "fraction_3"
])

dropped_miRNA = {
    "hsa-miR-1269a": "hsa-miR-1269a-3p",
    "hsa-miR-137": "hsa-miR-137-3p",
    "hsa-miR-217": "hsa-miR-217-5p",
    "hsa-miR-320a": "hsa-miR-320a-3p",
    "hsa-miR-375": "hsa-miR-375-3p",
    "hsa-miR-429": "hsa-miR-429-3p",
    "hsa-miR-451a": "hsa-miR-451a-5p"
}

for project in tqdm.tqdm(TCGA_projects):
    # Load isomiR expression table
    isomiR_CPM = pd.read_csv(
        os.path.join(TCGA_path, "isomiR_{}.tsv".format(project)),
        sep="\t", index_col=0
    )

    # Exponentiate (inverse log)
    isomiR_CPM = 2 ** isomiR_CPM - 1

    # Get highly expressed isomiRs
    quantile = isomiR_CPM.median(axis=1).quantile(expression_quantile)
    top_isomiRs = isomiR_CPM[isomiR_CPM.median(axis=1) >= quantile].index

    # Top isomiRs -> top miRNAs
    top_miRNAs = np.unique([isomiR.split("|")[0] for isomiR in top_isomiRs])

    # Get all isomiRs of top miRNAs
    top_isomiRs = [isomiR for isomiR in isomiR_CPM.index if isomiR.split("|")[0] in top_miRNAs]
    isomiR_CPM = isomiR_CPM.loc[top_isomiRs]

    # Identify the most expressed isoform and its fraction for each miRNA
    for miRNA in top_miRNAs:
        isomiRs = [isomiR for isomiR in isomiR_CPM.index if isomiR.startswith(miRNA)]
        # All possible shifts from 5'- and 3'-ends
        shifts_5 = np.unique([isomiR.split("|")[1] for isomiR in isomiRs])
        shifts_3 = np.unique([isomiR.split("|")[2] for isomiR in isomiRs])

        # For each 5'-shift, compute median expression over patients
        # by summing up over all possible 3'-shifts
        medians_5 = []
        for shift in shifts_5:
            median_expr = isomiR_CPM.loc[
                isomiR_CPM.index.str.startswith(miRNA + "|" + shift)
            ].sum(axis=0).median()
            medians_5.append(median_expr)

        # For each 3'-shift, compute median expression over patients
        # by summing up over all possible 5'-shifts
        medians_3 = []
        for shift in shifts_3:
            median_expr = isomiR_CPM.loc[
                (isomiR_CPM.index.str.startswith(miRNA)) & \
                (isomiR_CPM.index.str.endswith(shift))
            ].sum(axis=0).median()
            medians_3.append(median_expr)

        # Identify offsets of the most expressed 5'- and 3'-isomiRs
        max_shift_5 = shifts_5[np.argmax(medians_5)]
        max_shift_3 = shifts_3[np.argmax(medians_3)]

        # Calculate fractons of the most highly expressed isomiRs
        total_expr = isomiR_CPM.loc[(isomiR_CPM.index.str.startswith(miRNA))].sum(axis=0)
        expr_5 = isomiR_CPM.loc[
            isomiR_CPM.index.str.startswith(miRNA + "|" + max_shift_5)
        ].sum(axis=0)
        expr_3 = isomiR_CPM.loc[
            (isomiR_CPM.index.str.startswith(miRNA)) & \
            (isomiR_CPM.index.str.endswith(max_shift_3))
        ].sum(axis=0)

        fraction_5 = (expr_5 / total_expr).median()
        fraction_3 = (expr_3 / total_expr).median()

        if (miRNA in dropped_miRNA):
            miRNA = dropped_miRNA[miRNA]

        res.loc[len(res)] = [
            project, miRNA,
            max_shift_5, np.log2(np.max(medians_5) + 1), fraction_5,
            max_shift_3, np.log2(np.max(medians_3) + 1), fraction_3
        ]

res = res.sort_values(["Project", "miRNA"])
res.to_csv("../../processed_data/TCGA/main_isomiRs.tsv", sep="\t", index=None)
