import argparse
import os
import random

parser = argparse.ArgumentParser()
parser.add_argument("-if", "--input_file", required=True, help="input_file")
parser.add_argument("-p", "--probability", required=True, help="probability that the read will survive")
parser.add_argument("-r", "--reference_name", required = "True", help = "reference name")
args = parser.parse_args()

p = float(args.probability)
counter = 0
MaxInDel = 5
f1 = open(args.input_file, mode = "rU")

for i in range(1, int(MaxInDel) + 1):
    f1.seek(0)

    fout = open(args.reference_name + "_reg_indels_" + str(i) + ".fa", mode="w")
    print "writing indels" + str(i) + ".fa"

    for line_raw in f1:

        counter += 1

        if line_raw[0] == ">":
            JunctionName = line_raw.strip()
            JunctionSeq = ""

        else:
            JunctionSeq += line_raw.strip()

        if len(JunctionSeq) > 290:
            LeftExon = JunctionSeq[0:150]
            RightExon = JunctionSeq[150:300]
            #        print JunctionName
            # print JunctionSeq

            if random.random()<p:
                fout.write(JunctionName + "|DEL" + str(i) + "\n")
                fout.write(LeftExon[0:-i] + RightExon[i:] + "\n")

            if random.random()<p:
                InsertN = "N" * (i * 2)
                fout.write(JunctionName + "|INS" + str(i) + "\n")
                fout.write(LeftExon + InsertN + RightExon + "\n")

        if counter == 5000:
            fout.flush()
    fout.close()

f1.close()