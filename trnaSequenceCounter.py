#!/usr/bin/env python3
# Created by Dennis Mulligan

import sqlite3
import numpy as np
from subprocess import call

def getTrnaSeqs():

    conn = sqlite3.connect("trna.db")
    c = conn.cursor()

    # TABLE trna
    # [0] "(accession text, genus text, species text, chromosome text, "
    # [4] "number int, anticodon text, aminoacid text, "
    # [7] "begin int, end int, strand text, intron text, "
    # [11] "preflank text, sequence text, postflank text, fullsequence text, "
    # [15] "notes text, id int)")

    accessionlist = "'JX270811.1', 'MH559320.1', 'HG975452.1', 'AM087200.3'"

    outfilename = "temp.fasta"
    with open(outfilename,"w") as outfile:

        for row in c.execute("SELECT genus,species,accession,chromosome,aminoacid,postflank "
                             "FROM trna "
                             "WHERE "
                             "anticodon='GCA' "
                             "AND "
                             "(chromosome IN ('chloroplast','chromosome1') )" ):
                             #" OR (accession IN (" + accessionlist + ")))" ):
                             #"species IN ('lycopersicum', 'pennellii', 'annuum', 'baccatum')"):


            outfile.write(">" + ".".join(row[:-1]) + "\n")
            outfile.write(row[-1] + "\n")

    call(["clustalo", "-i", "temp.fasta",
          "-o", "temp-clustalo.fasta",
          "--guidetree-out=temp-tree.txt","--force"])


def consensusCounter():
    with open("temp-clustalo.fasta") as file:
        sequence = {}
        for line in file:
            if line[0] == ">":
                sequencename = line[1:]
                sequence[sequencename] = []
            else:
                sequence[sequencename].append(line.strip("\n"))

    first = True
    numSeqs = 0
    for key,value in sequence.items():

        numSeqs += 1
        seq = "".join(value)
#        print(seq)
        if first:
            first = False
            bases = ["A","C","G","T","-"]
            baseIndex =  {"A":0,"C":1,"G":2,"T":3,"-":4}
            countarray = np.zeros((len(seq),5),dtype=int)

        for i,b in enumerate(seq):
            countarray[i,baseIndex[b]] += 1

    consensus = "".join([bases[m] for m in np.argmax(countarray,axis=1)])
    print(consensus)
    np.savetxt("countarray.txt",countarray,delimiter="\t")







def main():
    getTrnaSeqs()
    consensusCounter()



if __name__ == "__main__":
    main()
