#!/usr/bin/env python3
# Created by Dennis Mulligan

from subprocess import call

def calltRNAscanSE():
    #genusList = "Vassobia Saracha Scopolia Przewalskia Physalis Lycium Datura Iochroma Eriolarynx Dunalia Atropa Acnistus Nicotiana Capsicum Solanum".split()
    #genusList = ["Cyanobacterium_aponium","Arthrospira_platensis","Chloroflexus_aurantiacus","E_coli"]
    #genusList = ["Capsicum_chinense_chr8","Capsicum_baccatum_chr8","Solanum_pennellii_chr8"]
    #
    genusList = ["Capsicum_baccatum_chr"+str(i) for i in range(10,13)]

    #genusList = ["Capsicum_annum_chr" + str(i) for i in [2,5]]

    searchmode = "-B"

    unixcalls = [["tRNAscan-SE",x+".fasta","-o","trnatable_"+x+".txt", searchmode] for x in genusList]
    for unixcall in unixcalls:
        print(unixcall)
        call(unixcall)


def main():
    calltRNAscanSE()

if __name__ == "__main__":
    main()