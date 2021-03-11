import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

OUTPUT_DIR = "../../../plots/main_isomiRs/shifts/"

term = int(sys.argv[1])
side = int(sys.argv[2])
tb = pd.read_csv("../../../processed_data/TCGA/main_isomiRs.tsv", sep="\t")
tb = tb[tb["miRNA"].str.endswith("{}p".format(side))]
# tb["miRNA"] = tb["miRNA"].str.lstrip("hsa-miR-")

tb = tb.pivot_table(
        index = "shift_{}".format(term),
        columns="miRNA", values="Project",
        aggfunc=len,
        fill_value=0
)

shifts = list(range(-2, 3))
# adding missed shifts
for shift in shifts:
    if shift not in list(tb.index):
        tb.loc[shift] = 0

for miRNA in tb.columns:
    tb.loc[:, miRNA] = tb.loc[:, miRNA] / tb.loc[:, miRNA].sum()

# drop too large shifts
tb = tb.drop(list(set(tb.index) - set(shifts)))

# make lexicographical order
tb = tb.sort_index()
tb = tb.reindex(sorted(tb.columns), axis=1)


cmap = sns.color_palette("Blues", as_cmap=True)
#p = sns.clustermap(
#    tb.T, dendrogram_ratio=(.05, .2), figsize=(7.5, 8),
#    cmap=cmap, linecolor="lightgrey", linewidth=0.05,
#    cbar_kws={"label": "fraction of TCGA projects"},
#    col_cluster=False
#)

#plt.setp(p.ax_heatmap.get_xticklabels(), rotation=0)
#plt.setp(p.ax_heatmap.get_yticklabels(), rotation=0)
#p.ax_row_dendrogram.set_visible(False)
#p.ax_heatmap.set_xlabel("")
#p.ax_heatmap.set_ylabel("")

#plt.tight_layout()
#plt.savefig(
#    OUTPUT_DIR + "{}term_{}p.png".format(term, side)
#)

fig, ax = plt.subplots(figsize=(8,8))
ax = sns.heatmap(
    tb.T,
    cmap=cmap, linecolor="lightgrey", linewidth=0.05,
    cbar_kws={"label": "Fraction of TCGA projects"},
    ax=ax
)

ax.set_xlabel("")
ax.set_ylabel("")
plt.tight_layout()
plt.savefig(
    OUTPUT_DIR + "{}term_{}p.png".format(term, side)
)




