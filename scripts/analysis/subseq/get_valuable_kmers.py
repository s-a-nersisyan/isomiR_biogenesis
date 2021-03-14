import numpy as np
import pandas as pd
import tqdm
import sys
import os

import scipy.stats as stats

FEATURE_TABLE_DIR = "../../../processed_data/TCGA/local_cleavage_nucleotides/"
OUTPUT_DIR = "../../../data/kmers_fisher/"
KMERS_MAX_LENGTH = 4
RADIUS = 4
PVALUE_TRESHOLD = 0.05

# extract kmers from all strings in "sequences"
def extract_kmers(sequences, k):
    kmers = set()
    for seq in sequences:
        for start in range(len(seq) - k):
            kmers.add(seq[start:start + k])
    return sorted(list(kmers))

# generate patterns by kmer
def generate_patterns(kmer):
    for i in range(1, len(kmer)):
        pattern = kmer[:i] + "*" + kmer[i:]
        yield pattern
    pattern = "*".join(list(kmer))
    yield pattern

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

for k in range(1, KMERS_MAX_LENGTH + 1):
    kmers = extract_kmers(tb["sequence"], k)
    out_tb = pd.DataFrame(columns=[
        "pattern",
        "heter_with",
        "heter_without",
        "nonheter_with",
        "nonheter_without",
        "pvalue"
    ])

    for kmer in tqdm.tqdm(kmers):
        for pattern in generate_patterns(kmer):
            a = len(tb[
                (tb["heterogenous"] == 1) &
                (tb["sequence"].str.contains(pattern, na=False, regex=True))
            ])
            b = len(tb[tb["heterogenous"] == 1]) - a
            c = len(tb[
                (tb["heterogenous"] == 0) &
                (tb["sequence"].str.contains(pattern, na=False, regex=True))
            ])
            d = len(tb[tb["heterogenous"] == 0]) - c

            _, pv = stats.fisher_exact([[a, b], [c, d]])
            out_tb.loc[len(out_tb)] = [
                pattern, a, b, c, d, pv
            ]

    out_tb.to_csv(os.path.join(OUTPUT_DIR, "{}term_{}p_{}_{}.tsv".format(
        termini[0], side[0], gap_flag, k
    )), sep="\t", index=None)
