import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

fraction_quantile = 0.99
fill_value = -10

termini = int(sys.argv[1])
prime = int(sys.argv[2])
tb = pd.read_csv("../../processed_data/TCGA/main_isomiRs.tsv", sep="\t")
tb = tb[tb["miRNA"].str.endswith("{}p".format(termini))]

shift_tb =  tb.pivot_table(
    index="Project",
    columns="miRNA",
    values="shift_{}".format(prime),
    fill_value=fill_value
)
shift_tb = shift_tb.astype(int)

for miRNA in shift_tb.columns:
    shifts, counts = np.unique(shift_tb[miRNA], return_counts=True)
    fill_value_ind = 0
    while (fill_value_ind < len(shifts) and
           shifts[fill_value_ind] != fill_value):
        fill_value_ind += 1
    if (fill_value_ind < len(shifts)):
        shifts = np.delete(shifts, fill_value_ind)
        counts = np.delete(counts, fill_value_ind)
    # print(shifts, counts, miRNA)

    max_shift_fraction = max(counts) / sum(counts)
    if (max_shift_fraction < fraction_quantile):
        projects = np.array(shift_tb.index)[(shift_tb[miRNA] != fill_value)]
        print(miRNA, len(projects), shifts, counts, max_shift_fraction)


