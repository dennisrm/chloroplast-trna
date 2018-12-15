#!/usr/bin/env python3
# Created by Dennis Mulligan

import sqlite3
from subprocess import call

import argparse

def main():
    anticodons = set("TCT GCA GGT GCC GGA GAA CAT GAC ACG TAG GTT CAA CAT GAC ACG GTT "
                     "CAA CAT TGG CCA TGT CAT TGA TTC GTA GTC GCT TTG GTG".split())


    conn = sqlite3.connect("trna.db")
    c = conn.cursor()

    accessionlist = "JX270811.1, MH559320.1, HG975452.1, AM087200.3"


    for xxx in anticodons:
        outfilename = "temp.fasta"
        with open(outfilename,"w") as outfile:
            speciesset = set()
            for row in c.execute("SELECT id,genus,species,accession,chromosome,aminoacid,fullsequence "
                                 "FROM trna "
                                 "WHERE anticodon='"+xxx+"' "
                                 "AND (chromosome IN ('chloroplast','bacteria') "
                                 "OR (accession IN ('JX270811.1', 'MH559320.1', 'HG975452.1', 'AM087200.3')))" ):
                                 #"CASE chromosome "
                                 #"WHEN 'chloroplast' THEN 0 "
                                 #"WHEN 'bacteria' THEN 1 "
                                 #"WHEN NOT ('bacteria' OR 'chloroplast') THEN 2 "
                                 #"END"):

                if row[2] not in speciesset:
                    outfile.write(">" + ".".join([str(x) for x in [row[1],row[2]]]) + "\n")
                    outfile.write(row[-1] + "\n")
                    speciesset.add(row[2])
        call(["clustalo", "-i", "temp.fasta",
              "-o", "codons/clustalo-"+xxx+".fasta",
              "--guidetree-out=codons/tree-"+xxx+".txt","--force"])


    conn.close()



if __name__ == "__main__":
    main()