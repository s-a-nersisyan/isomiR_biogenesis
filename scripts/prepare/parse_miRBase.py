import numpy as np
import pickle

MIRBASE_PATH = "../../miRBase/miRNA.str"

hierarchy = dict()

mir_file = open(MIRBASE_PATH, "r")
lines = mir_file.readlines()

index = 0
while (index < len(lines)):
    line = lines[index]
    # find header of the next miRNA
    while not (line.startswith(">")):
        index += 1
        line = lines[index]

    # extract the second structure of miRNA
    line = line[1:-1]

    miRNA = line.split(" ")[0]
    description = line
    index += 2


    line = lines[index]
    second_structure = [description]
    while not (line.startswith("\n")):
        if not ("|" in line):
            second_structure.append(line[:-1])
        index += 1
        line = lines[index]

    hierarchy[miRNA] = second_structure
    index += 1

with open("../../miRBase/parsed_miRBase.pkl", "wb") as handle:
    pickle.dump(hierarchy, handle, protocol=pickle.HIGHEST_PROTOCOL)
