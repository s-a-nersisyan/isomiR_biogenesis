import numpy as np
import pandas as pd
import tqdm
import sys
import os

import scipy.stats as stats

FEATURE_TABLE_DIR = "../../../processed_data/TCGA/local_cleavage_nucleotides/"
OUTPUT_DIR = "../../../data/subseq_fisher/"
SUBSEQ_MAX_LENGTH = 4
RADIUS = 4
PVALUE_TRESHOLD = 0.05

# geneartes all strictly monotone sequences
# of size "length" with values in [0, "vmax"]
def sequence_iter(vmax, length):
    sequence = [-1]
    while (len(sequence) > 0):
        sequence[-1] += 1
        if (len(sequence) >= length):
            yield sequence

        if (sequence[-1] >= vmax):
            sequence.pop()
        elif (len(sequence) < length):
            sequence.append(sequence[-1])

termini = sys.argv[1]
side = sys.argv[2]
gap_flag = int(sys.argv[3])
if (gap_flag):
    gap_flag = "gap"
else:
    gap_flag = "ungap"

# table of local cleavage nucleotides
tb_path = os.path.join(FEATURE_TABLE_DIR, "{}term_{}p_{}_{}.tsv".format(
    termini[0], side[0], gap_flag, RADIUS
))
tb = pd.read_csv(tb_path, sep="\t")

# merge source and opposite local cleavage nucleotides
tb["sequence"] = (tb["source"] + tb["opposite"]).str.lower()
sequences = np.array([list(seq) for seq in tb["sequence"]])
for length in range(1, SUBSEQ_MAX_LENGTH + 1):
    out_tb = pd.DataFrame(columns=[
        "indexes",
        "nucleotides",
        "heter_with",
        "heter_without",
        "nonheter_with",
        "nonheter_without",
        "pvalue"
    ])

    for indexes in tqdm.tqdm(sequence_iter(4 * RADIUS - 1, length)):
        # get all nucleotide subsequences with nucleotides on indexes position
        tb["subsequence"] = ["".join(subseq) for subseq in sequences[:, indexes]]
        for subseq in set(tb["subsequence"]):
            a = len(tb[(tb["heterogenous"] == 1) & (tb["subsequence"] == subseq)])
            b = len(tb[tb["heterogenous"] == 1]) - a
            c = len(tb[(tb["heterogenous"] == 0) & (tb["subsequence"] == subseq)])
            d = len(tb[tb["heterogenous"] == 0]) - c

            _, pv = stats.fisher_exact([[a, b], [c, d]])
            out_tb.loc[len(out_tb)] = [
                ",".join(str(i) for i in indexes),
                subseq,a, b, c, d, pv
            ]

    out_tb.to_csv(os.path.join(OUTPUT_DIR, "{}term_{}p_{}_{}.tsv".format(
        termini[0], side[0], gap_flag, length
    )), sep="\t", index=None)
