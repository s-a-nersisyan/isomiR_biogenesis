import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

OUTPUT_DIR = "../../../plots/main_isomiRs/fractions/"

term = int(sys.argv[1])
side = int(sys.argv[2])
tb = pd.read_csv("../../../processed_data/TCGA/main_isomiRs.tsv", sep="\t")
tb = tb[tb["miRNA"].str.endswith("{}p".format(side))]
# tb["miRNA"] = tb["miRNA"].str.lstrip("hsa-miR-")

tb = tb.pivot_table(
    index = "Project", columns="miRNA",
    values = "fraction_{}".format(term),
    fill_value=0
)

# make lexicographical order
tb = tb.sort_index()
tb = tb.reindex(sorted(tb.columns), axis=1)

cmap = sns.color_palette("Blues", as_cmap=True)
#p = sns.clustermap(
#    tb.T, dendrogram_ratio=(.05, .2),
#    cmap=cmap, linecolor="lightgrey", linewidth=0.1,
#    cbar_kws={"label": "fraction of the main isomiR"},
#)

#plt.setp(p.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
#plt.setp(p.ax_heatmap.get_yticklabels(), rotation=0)
#p.ax_heatmap.set_xlabel("")
#p.ax_heatmap.set_ylabel("")

fig, ax = plt.subplots(figsize=(8,8))
ax = sns.heatmap(
    tb.T,
    cmap=cmap, linecolor="lightgrey", linewidth=0.05,
    cbar_kws={"label": "Fraction of the main isomiR"},
    ax=ax
)

plt.setp(ax.get_xticklabels(), rotation=90, ha="center")
plt.setp(ax.get_yticklabels(), rotation=0)
ax.set_xlabel("")
ax.set_ylabel("")
plt.tight_layout()

plt.savefig(
    OUTPUT_DIR + "{}term_{}p.png".format(term, side)
)
