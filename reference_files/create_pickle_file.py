import pickle

import argparse
import pickle

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', help='fasta input file')
args = parser.parse_args()

junction_id_to_seq = {}
with open(args.fasta, "r") as f:
    while True:
        line1 = f.readline()
        if not line1:
            break
        line2 = f.readline()

        junction_id_to_seq[line1.strip()] = line2.strip()

pickle.dump(junction_id_to_seq, open("known_fusions.pickle", "wb"))
