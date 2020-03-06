import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches


class Genome:

    def __init__(self, header, sequence, chromosome, getMiniMap=False):
        self.header = header
        hdsplit = header.split()
        self.accession = hdsplit[0][1:]
        self.genus = hdsplit[1]
        self.species = hdsplit[2]
        self.sequence = sequence
        self.chromosome = chromosome
        if getMiniMap:
            self.miniMap = MiniMap(sequence)
            print(len(self.miniMap))


class MiniMap(dict):

    def __init__(self, string, wkt=None):
        """Finds minimizerMap from string using wkt=(window, kmer length)"""
        # dict.__init__(self)

        w, k, t = (30, 13, 0)
        lastMinimerIndex = -1
        for i in range(len(string) - w + 1):
            minimer = "Z"
            for j in range(w - k + 1):
                kmer = string[i + j: i + j + k]
                if kmer < minimer:
                    minimer = kmer
                    currentMinimerIndex = i + j

            if lastMinimerIndex != currentMinimerIndex:
                lastMinimerIndex = currentMinimerIndex
                if minimer in self:
                    self[minimer].append(currentMinimerIndex)
                else:
                    self[minimer] = [currentMinimerIndex]

    def checkMap(self, other):
        for minimer, positionList in self.items():
            if minimer in other:
                for position in positionList:
                    yield (position, tuple(other[minimer]))


def faReader(filename, getMiniMap=False):
    """Save each sequence from fasta file as ChloroplastGenome objects"""
    with open(filename, "r") as fasta:
        sequence = []
        genlist = []
        chromosome = filename.split("_")[-1].split(".")[0]
        if chromosome == "genbank":
            chromosome = "cp"
        for line in fasta.readlines():
            if line[0] == ">":
                if sequence:
                    genlist.append(Genome(header, "".join(sequence), chromosome, getMiniMap))
                    sequence = []
                header = line
            else:
                sequence.append(line)
        genlist.append(Genome(header, "".join(sequence), chromosome, getMiniMap))
        return genlist


class tRNA:
    def __init__(self, row):
        attributes = "accession,genus,species,chromosome,number,anticodon,aminoacid,begin,end,strand,intron,preflank,sequence,postflank,fullsequence,notes,id"
        for i, attr in enumerate(attributes.split(",")):
            setattr(self, attr, row[i])


import sqlite3


def gettRNAdb():
    """  """
    # Connecting to sqlite databse
    conn = sqlite3.connect("trna.db")
    c = conn.cursor()

    # SELECT from this table:
    #          "CREATE TABLE trna"
    #          "(accession text, genus text, species text, chromosome text, "
    #          "number int, anticodon text, aminoacid text, "
    #          "begin int, end int, strand text, intron text, "
    #          "preflank text, sequence text, postflank text, fullsequence text, "
    #          "notes text, id int)"

    trnaList = []

    for row in c.execute("SELECT * FROM trna "):
        trnaList.append(tRNA(row))

    conn.close()

    return trnaList


def checkConnection(clust1, clust2, L=60):
    for x1, y1 in clust1:
        for x2, y2 in clust2[-1:]:
            if abs(x2 - x1) <= L and abs(y2 - y1) <= L:
                return True
    return False


def cluster(small, big):
    cludic = {}
    incount, outcount = 0, 0
    intotal, outtotal = 0, 0
    for position, entries in small.miniMap.checkMap(big.miniMap):
        if len(entries) < 40:
            incount += 1
            intotal += len(entries)
            x = position

            for y in entries:
                yslash = y // 1e7
                if yslash not in cludic:
                    cludic[yslash] = []
                cludic[yslash].append([(x, y)])
        else:
            outcount += 1
            outtotal += len(entries)

    print("in: ", incount, intotal)
    print("out: ", outcount, outtotal)
    print(cludic.keys())

    components = []
    for key in cludic.keys():
        clusters = cludic[key]
        while clusters:
            clust1 = clusters.pop()

            for clust2 in clusters:
                if checkConnection(clust1, clust2):
                    clust2.extend(clust1)
                    break
            else:
                if len(clust1) >= 2:
                    components.append(clust1)
                    print(len(clusters), " : ", len(components), " : ", len(clust1))

    xs = []
    ys = []
    for com in components:
        for x, y in com:
            xs.append(x)
            ys.append(y)
    return (xs, ys)

