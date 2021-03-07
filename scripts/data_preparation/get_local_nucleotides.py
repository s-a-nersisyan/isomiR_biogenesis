import pandas as pd
import numpy as np
import pickle
import sys
import os

OUTPUT_DIR = "../../processed_data/TCGA/isomiR_local_cleavage_nucleotides/"
MIRBASE_PATH = "../../processed_data/miRBase/miRBase.pkl"
with open(MIRBASE_PATH, 'rb') as handle:
    miRBase = pickle.load(handle)

# for pri_miRNA in miRBase:
#     print(pri_miRNA)

MAIN_CONSENSUS_MIRNA_PATH = "../../processed_data/TCGA/main_consensus_isomiRs.tsv"
consensus_tb = pd.read_csv(MAIN_CONSENSUS_MIRNA_PATH, sep="\t")

def parse_description(description):
    miRNA = description.split(" ")[0]
    side_count = len(description.split("[")) - 1

    sides = []
    for ind in range(side_count):
        side_desc = description.split("[")[1 + ind].split("]")[0]
        side = side_desc.split(":")[0]
        indexes = [int(el) - 1 for el in side_desc.split(":")[1].split("-")]
        sides.append((side, indexes[0], indexes[1]))
    return sides

def merge_cycle_stem(cycle, stem):
    result = ""
    for ind in range(len(stem)):
        if not (stem[ind] == " "):
            result += stem[ind]
        else:
            result += cycle[ind]
    return result

def get_local_nucleotides(pri_miRNA, side, termini, shift, radius, gap_flag=False):
    # rotate second structure so that side termini position
    # is on the upper left part of miRNA structure
    second_structure = miRBase[pri_miRNA]
    if (side == "3p"):
        second_structure[0],second_structure[-1] = second_structure[-1], second_structure[0]
        second_structure[1], second_structure[-2] = second_structure[-2], second_structure[1]
        second_structure = [line[::-1] for line in second_structure]
    if (termini == "3prime"):
        second_structure = [line[::-1] for line in second_structure]

    upp_side = merge_cycle_stem(second_structure[0], second_structure[1])
    dwn_side = merge_cycle_stem(second_structure[-1], second_structure[-2])

    # find index of the canonicial cleavage position
    canon_cleavage_pos = -1
    for ind in range(len(upp_side)):
        if upp_side[ind].isupper():
            canon_cleavage_pos = ind
            break
    if (canon_cleavage_pos == -1):
        # return -1
        return "cleavage position is not found"

    # find index of the main (shifted) cleavage position
    main_cleavage_pos = canon_cleavage_pos
    if (abs(shift) > 0):
        direction = shift // abs(shift)
    nucleotides_count = 0
    while (nucleotides_count < abs(shift)):
        if (upp_side[main_cleavage_pos] != "-"):
            nucleotides_count += 1
        main_cleavage_pos += direction

    # extract cleavage local nucleotides with gaps
    upp_gap_wrd = upp_side[max(0, main_cleavage_pos - radius):
                           min(len(upp_side), main_cleavage_pos + radius)]
    dwn_gap_wrd = dwn_side[max(0, main_cleavage_pos - radius):
                           min(len(upp_side), main_cleavage_pos + radius)][::-1]

    # extract cleavage local nucleotides without gaps for the upp side
    start_pos = main_cleavage_pos - 1
    nucleotides_count = 0
    while (start_pos >= 0 and nucleotides_count < radius):
        if (upp_side[start_pos] != "-"):
            nucleotides_count += 1
        start_pos -= 1
    start_pos += 1
    if (start_pos < 0):
        return "upp side start position is not founded"

    end_pos = main_cleavage_pos
    nucleotides_count = 0
    while (end_pos < len(upp_side) and nucleotides_count < radius):
        if (upp_side[end_pos] != "-"):
            nucleotides_count += 1
        end_pos += 1
    if (end_pos > len(dwn_side)):
        return "upp side end position is not founded"
    upp_ungap_wrd = upp_side[start_pos:end_pos]
    upp_ungap_wrd = "".join([nucl for nucl in upp_ungap_wrd if (nucl != "-")])

    # extract cleavage local nucleotides without gaps for the dwn side
    start_pos = main_cleavage_pos - 1
    nucleotides_count = 0
    while (start_pos >= 0 and nucleotides_count < radius):
        if (dwn_side[start_pos] != "-"):
            nucleotides_count += 1
        start_pos -= 1
    start_pos += 1
    if (start_pos < 0):
        return "dwn side start position is not founded"

    end_pos = main_cleavage_pos
    nucleotides_count = 0
    while (end_pos < len(dwn_side) and nucleotides_count < radius):
        if (dwn_side[end_pos] != "-"):
            nucleotides_count += 1
        end_pos += 1
    if (end_pos > len(dwn_side)):
        return "dwn side end position is not founded"
    dwn_ungap_wrd = dwn_side[start_pos:end_pos][::-1]
    dwn_ungap_wrd = "".join([nucl for nucl in dwn_ungap_wrd if (nucl != "-")])

    return [upp_gap_wrd, dwn_gap_wrd,
            upp_ungap_wrd, dwn_ungap_wrd]

terminis = ["5prime", "3prime"]
radiuses = list(range(1, 5))

for radius in radiuses:
    for termini in terminis:
        table_upp_gap = pd.DataFrame(columns=[
            "pri_miRNA",
            "miRNA",
            "sequence",
            "heterogenous"
        ])

        table_dwn_gap = pd.DataFrame(columns=[
            "pri_miRNA",
            "miRNA",
            "sequence",
            "heterogenous"
        ])

        table_upp_ungap = pd.DataFrame(columns=[
            "pri_miRNA",
            "miRNA",
            "sequence",
            "heterogenous"
        ])

        table_dwn_ungap = pd.DataFrame(columns=[
            "pri_miRNA",
            "miRNA",
            "sequence",
            "heterogenous"
        ])

        for index, miRNA_row in consensus_tb.iterrows():
            # extracting pri_miRNA identifier
            miRNA = miRNA_row["miRNA"]
            shift = miRNA_row["shift_" + termini[0]]
            heterogenous = miRNA_row["heterogenous_" + termini[0]]
            if (miRNA.rfind("-3p")):
                side = "3p"
            elif (miRNA.rfind("-5p")):
                side = "5p"
            else:
                print("miRNA is not sided")
            pri_miRNA = miRNA[:miRNA.rfind("-")].lower()

            # if pri_miRNA in miRBase then extract its
            # cleavage local nucleotides
            if (pri_miRNA in miRBase):
                features = get_local_nucleotides(pri_miRNA, side, termini, shift, radius)
                table_upp_gap.loc[len(table_upp_gap)] = [
                    pri_miRNA,
                    miRNA,
                    features[0],
                    heterogenous
                ]

                table_dwn_gap.loc[len(table_dwn_gap)] = [
                    pri_miRNA,
                    miRNA,
                    features[1],
                    heterogenous
                ]

                table_upp_ungap.loc[len(table_upp_ungap)] = [
                    pri_miRNA,
                    miRNA,
                    features[2],
                    heterogenous
                ]

                table_dwn_ungap.loc[len(table_dwn_ungap)] = [
                    pri_miRNA,
                    miRNA,
                    features[3],
                    heterogenous
                ]

                continue

            # otherwise miRBase consider other stem pri_miRNA variants of miRNA
            # extrect them all
            for variant in miRBase:
                if variant.startswith(pri_miRNA + "-"):
                    features = get_local_nucleotides(variant, side, termini, shift, radius)
                    table_upp_gap.loc[len(table_upp_gap)] = [
                        pri_miRNA,
                        miRNA,
                        features[0],
                        heterogenous
                    ]

                    table_dwn_gap.loc[len(table_dwn_gap)] = [
                        pri_miRNA,
                        miRNA,
                        features[1],
                        heterogenous
                    ]

                    table_upp_ungap.loc[len(table_upp_ungap)] = [
                        pri_miRNA,
                        miRNA,
                        features[2],
                        heterogenous
                    ]

                    table_dwn_ungap.loc[len(table_dwn_ungap)] = [
                        pri_miRNA,
                        miRNA,
                        features[3],
                        heterogenous
                    ]

        table_upp_gap.to_csv(OUTPUT_DIR + "{}'_upp_gap_{}.tsv".format(
            termini[0], radius
        ), sep="\t", index=None)

        table_dwn_gap.to_csv(OUTPUT_DIR + "{}'_dwn_gap_{}.tsv".format(
            termini[0], radius
        ), sep="\t", index=None)

        table_upp_ungap.to_csv(OUTPUT_DIR + "{}'_upp_ungap_{}.tsv".format(
            termini[0], radius
        ), sep="\t", index=None)

        table_dwn_ungap.to_csv(OUTPUT_DIR + "{}'_dwn_ungap_{}.tsv".format(
            termini[0], radius
        ), sep="\t", index=None)
