#!/usr/bin/env python3
# Created by Dennis Mulligan

# Process tRNAscan-SE tables, construct tRNA databse,

import sqlite3


## Object for storing tRNAscan-SE row information representing a tRNA gene
class tRNA:
    def __init__(self,line,seqdic,chromosome):
        self.accession = line[0]
        self.number = line[1]
        self.begin = int(line[2])
        self.end = int(line[3])
        if self.begin < self.end:
            self.strand = "+"
            self.plus = True
        else:
            self.strand = "-"
            self.plus = False
        self.aminoacid = line[4]
        self.anticodon = line[5]
        self.intron = line[6]+":"+line[7]
        self.score = line[8]
        if len(line) >= 10:
            self.notes = line[9]
        else:
            self.notes = "0"

        self.chromosome = chromosome
        self.organism,accessionseq = seqdic[self.accession]

        self.genus,self.species = self.organism

        if self.plus and self.intron == "0:0":
            self.seq = accessionseq[self.begin-1:self.end]
            self.preflank = accessionseq[self.begin-101:self.begin]
            self.postflank = accessionseq[self.end:self.end+101]
            self.fullseq = self.preflank + self.seq + self.postflank
        elif self.intron == "0:0":
            self.seq = revcomp(accessionseq[self.end - 1:self.begin])
            self.preflank = revcomp(accessionseq[self.begin:self.begin+101])
            self.postflank = revcomp(accessionseq[self.end-101:self.end])
            self.fullseq = self.preflank + self.seq + self.postflank
        else:
            self.seq,self.preflank,self.postflank,self.fullseq = "intron","intron","intron","intron"

def revcomp(seq):
    compdic = {"A": "T", "T": "A", "G": "C", "C": "G","Y":"R","R":"Y","N":"N"}
    newseq = [compdic[x] for x in seq]
    newseq.reverse()
    return "".join(newseq)

def fastaReader(filename,seqdic):
    with open(filename) as f:
        first = True
        for line in f:
            if line[0] == ">":
                if first:
                    first = False
                else:
                    seqdic[accession] = [organism,"".join(seqlines)]
                accession = line.split()[0][1:]
                organism = line.split()[1:3]
                seqlines = []
            else:
                seqlines.append(line[:-1])
        seqdic[accession] = [organism,"".join(seqlines)]
        return seqdic


def readTable(filename,trnalist,seqdic,chromosome):
    with open(filename) as f:
        table = f.read().splitlines()
        for line in table[3:]:
            line = line.split()
            trnalist.append(tRNA(line,seqdic,chromosome))
    return trnalist



def getTrna():
    bacteriaList = ["Cyanobacterium_aponium","Arthrospira_platensis","Chloroflexus_aurantiacus","E_coli"]
    chloroplastList = "Vassobia Saracha Scopolia Przewalskia Physalis Lycium Datura Iochroma Eriolarynx Dunalia Atropa Acnistus Nicotiana Capsicum Solanum".split()
    chr8List=["Capsicum_chinense_chr8"]
    genomeList = ["Capsicum_annuum","Capsicum_baccatum","Solanum_lycopersicum","Solanum_pennellii"]

    # seqdic = {}
    trnalist = []

    print("bacteria...")
    for bacteria in bacteriaList:
        seqdic = fastaReader(bacteria+".fasta",{})
        trnalist = readTable("trnatable_"+bacteria+".txt",trnalist,seqdic, "bacteria")

    print("chloroplasts...")
    for genus in chloroplastList:
        seqdic = fastaReader(genus+"_cp_genbank.fasta",{})
        trnalist = readTable("trnatable_"+genus+".txt",trnalist,seqdic, "chloroplast")

    for species in genomeList:
        print(species+"...")
        for num in [str(x) for x in range(1,13)]:
            seqdic = fastaReader(species+"_chr"+num+".fasta", {})
            trnalist = readTable("trnatable_" + species+"_chr"+num+".txt", trnalist, seqdic, "chromosome"+num)
    return trnalist



def dbMaker(trnalist):

    conn = sqlite3.connect("trna.db")
    c = conn.cursor()
    c.execute("CREATE TABLE trna"
              "(accession text, genus text, species text, chromosome text, "
              "number int, anticodon text, aminoacid text, "
              "begin int, end int, strand text, intron text, "
              "preflank text, sequence text, postflank text, fullsequence text, "
              "notes text, id int)")

    rowlist = []
    i = 0
    for tr in trnalist:
        i += 1
        rowlist.append((tr.accession,tr.genus,tr.species,tr.chromosome,
                        tr.number,tr.anticodon,tr.aminoacid,
                        tr.begin,tr.end,tr.strand,tr.intron,
                        tr.preflank,tr.seq,tr.postflank,tr.fullseq,
                        tr.notes, i))

    c.executemany("INSERT INTO trna VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",rowlist)

    conn.commit()
    conn.close()



def main():
    trnalist = getTrna()
    dbMaker(trnalist)

if __name__ == "__main__":
    main()

