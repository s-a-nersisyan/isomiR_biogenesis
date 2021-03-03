import pickle
import numpy as np


lines = open("../../raw_data/miRBase/miRNA.str").readlines()

pri_miRNA_info = dict()
i = 0
while i < len(lines):
    line = lines[i]
    # Find header of the next pri_miRNA
    while not line.startswith(">"):
        i += 1
        line = lines[i]

    # Remove trailing characters
    line = line[1:-1]

    pri_miRNA = line.split(" ")[0]
    i += 2

    # Extract secondary structure of pri_miRNA
    #line = lines[i]
    #secondary_structure = [line]
    secondary_structure = []
    line = lines[i]
    while line != "\n":
        if "|" not in line:
            secondary_structure.append(line[:-1])

        i += 1
        line = lines[i]

    pri_miRNA_info[pri_miRNA] = secondary_structure
    i += 1

with open("../../processed_data/miRBase/miRBase.pkl", "wb") as out_file:
    pickle.dump(pri_miRNA_info, out_file)
