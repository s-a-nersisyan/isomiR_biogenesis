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
tb["miRNA"] = tb["miRNA"].str.lstrip("hsa-miR-")

tb = tb.pivot_table(
        index = "shift_{}".format(term),
        columns="miRNA", values="Project",
        aggfunc=len,
        fill_value=0
)

for miRNA in tb.columns:
    tb.loc[:, miRNA] = tb.loc[:, miRNA] / tb.loc[:, miRNA].sum()

cmap = sns.color_palette("Blues")
p = sns.clustermap(
    tb.T, dendrogram_ratio=(.05, .2),
    cmap=cmap, linecolor="lightgrey", linewidth=0.1,
    cbar_kws={"label": "fraction of TCGA projects"},
    col_cluster=False
)

plt.setp(p.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
plt.setp(p.ax_heatmap.get_yticklabels(), rotation=0)
p.ax_heatmap.set_xlabel("")
p.ax_heatmap.set_ylabel("")
# p.ax_heatmap.set_yticklabels("")
# plt.title("{}term_{}p.png".format(term, side), loc="left", fontdict={"fontsize": "xx-large", "fontweight": "bold"})

plt.tight_layout()
plt.savefig(
    OUTPUT_DIR + "{}term_{}p.png".format(term, side)
)
