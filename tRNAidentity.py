#!/usr/bin/env python3
# Created by Dennis Mulligan

import sqlite3
import numpy as np
from subprocess import call

def getTrnaSeqs(xxx,aaa):

    conn = sqlite3.connect("trna.db")
    c = conn.cursor()

    # TABLE trna
    # [0] "(accession text, genus text, species text, chromosome text, "
    # [4] "number int, anticodon text, aminoacid text, "
    # [7] "begin int, end int, strand text, intron text, "
    # [11] "preflank text, sequence text, postflank text, fullsequence text, "
    # [15] "notes text, id int)")

    accessionlist = "'JX270811.1', 'MH559320.1', 'HG975452.1', 'AM087200.3'"


    names,preflanks,sequences,postflanks = [],[],[],[]

    for row in c.execute("SELECT genus,species,accession,chromosome,aminoacid,intron,preflank,sequence,postflank "
                         "FROM trna "
                         "WHERE "
                         "anticodon='"+
                         xxx+"' "
                             "AND aminoacid='"+
                         aaa+"' "
                             "AND "
                             "(chromosome IN ('bacteria','chromosome1','chromosome2','chromosome3','chromosome4',"
                             "'chromosome5','chromosome6','chromosome7','chromosome8','chromosome9','chromosome10',"
                             "'chromosome11','chromosome12') "
                             " OR (accession IN (" + accessionlist + "))) "
                                                                     "ORDER BY chromosome='chloroplast' DESC"):
        #"species IN ('lycopersicum', 'pennellii', 'annuum', 'baccatum')"):

        if row[-4] == "0:0":
            names.append(".".join(row[:-3]))
            preflanks.append(row[-3])
            sequences.append(row[-2])
            postflanks.append(row[-1])

    ##########################
    outfilename = "temp.fasta"
    with open(outfilename, "w") as outfile:
        for i in range(len(preflanks)):
            outfile.write(">" + names[i] + "\n")
            outfile.write(preflanks[i] + "\n")

    call(["clustalo", "-i", "temp.fasta",
          "-o", "temp-clustalo-pre.fasta",
          "--guidetree-out=temp-tree-pre.txt","--force"])

    ###########################
    outfilename = "temp.fasta"
    with open(outfilename, "w") as outfile:
        for i in range(len(sequences)):
            outfile.write(">" + names[i] + "\n")
            outfile.write(sequences[i] + "\n")

    call(["clustalo", "-i", "temp.fasta",
    "-o", "temp-clustalo-seq.fasta",
    "--guidetree-out=temp-tree-seq.txt", "--force"])

    ############################
    outfilename = "temp.fasta"
    with open(outfilename, "w") as outfile:
        for i in range(len(postflanks)):
            outfile.write(">" + names[i] + "\n")
            outfile.write(postflanks[i] + "\n")

    call(["clustalo", "-i", "temp.fasta",
    "-o", "temp-clustalo-post.fasta",
    "--guidetree-out=temp-tree-post.txt", "--force"])


def chloroplastConsensus(seqtype,aa,anticodon):
    with open("temp-clustalo-"+seqtype+".fasta") as file:
        sequence = {}
        i = 0
        for line in file:
            if line[0] == ">":
                i += 1
                sequencename = line[1:].strip("\n")+"."+aa+anticodon+str(i)
                sequence[sequencename] = []
            else:
                sequence[sequencename].append(line.strip("\n"))

    first = True
    numSeqs = 0
    for key,value in sequence.items():
        if "chloroplast" in key:
            numSeqs += 1
            seq = "".join(value)
            if first:
                first = False
                bases = ["A","C","G","T","-"]
                baseIndex =  {"A":0,"C":1,"G":2,"T":3,"-":4}
                countarray = np.zeros((len(seq),5),dtype=int)

            for i,b in enumerate(seq):
                countarray[i,baseIndex[b]] += 1

    consensus = "".join([bases[m] for m in np.argmax(countarray,axis=1)])
    return sequence,consensus

def identityCounter(sequence,consensus,outdic):
    """Writes the percent identity of sequences vs chloroplast consensus"""
    first = True
    id = 0
    for key, value in sequence.items():
        seq = "".join(value)
        if first:
            first = False
            bases = ["A", "C", "G", "T", "-"]
            baseIndex = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}
            countarray = np.zeros((len(seq), 5), dtype=int)

        identity = 0
        for i, b in enumerate(seq):
            if b == consensus[i]:
                identity += 1
        identity = identity/(i+1)
        if key in outdic:
            outdic[key].append(str(identity))
        else:
            outdic[key] = [str(identity)]



def main():
    anticodonlist = {('TCT','Arg'),
                     ('GCA','Cys'),
                     ('GGT','Thr'),
                     ('GCC','Gly'),
                     ('GGA','Ser'),
                     ('GAA','Phe'),
                     ('CAT','Met'),
                     ('GAC','Val'),
                     #('ACG','Arg'),
                     ('TAG','Leu'),
                     ('GTT','Asn'),
                     ('CAA','Leu'),
                     ('CAT','Ile2'),
                     ('GAC','Val'),
                     #('ACG','Arg'),
                     ('GTT','Asn'),
                     ('CAA','Leu'),
                     ('CAT','Ile2'),
                     ('TGG','Pro'),
                     #('CCA','Prp'),
                     ('TGT','Thr'),
                     ('CAT','fMet'),
                     ('TGA','Ser'),
                     ('TTC','Glu'),
                     ('GTA','Tyr'),
                     ('GTC','Asp'),
                     ('GCT','Ser'),
                     ('TTG','Gln') }
                     #'GTG','His')}



    outdic = {}

    for anticodon,aa in anticodonlist:
        print(anticodon,aa)
        getTrnaSeqs(anticodon,aa)

        preSeq,preConsensus = chloroplastConsensus("pre",aa,anticodon)
        seqSeq,seqConsensus = chloroplastConsensus("seq",aa,anticodon)
        postSeq,postConsensus = chloroplastConsensus("post",aa,anticodon)

        identityCounter(preSeq,preConsensus,outdic)
        identityCounter(seqSeq,seqConsensus,outdic)
        identityCounter(postSeq,postConsensus,outdic)

    outfilename = "temp-identity.txt"
    with open(outfilename,"w") as outfile:
        outfile.write("preid\tseqid\tpostid\tgenus\tspecies\taccesion\tdot\tchromosome\tamino\tintron\ti\n")
        for key,value in outdic.items():
            outfile.write("\t".join(value+key.split(sep="."))+"\n")




if __name__ == "__main__":
    main()
