#
"""
Get minimizer maps from genomes and save to sqlite database
Dennis Mulligan
"""
import time
start = time.time()

from minimapdb import *
from soltools import *


def speciesFaToMapDb(genus, species, create=True):
    tabname = genus[0]+species[:2]
    with ConnectDB("db6.db") as mapdb:

        print("time: ", time.time() - start)

        if create:
            genomes = faReader(f"./genomes/{genus}_cp_genbank.fasta", True)
            print("time: ", time.time() - start)
            mapdb.create(tabname,"minimer")
            for gen in genomes:
                if gen.species==species:
                    name = "".join([gen.genus[0], gen.species[:2], "_", gen.chromosome, "_", gen.accession.replace(".", "")])
                    print(name)
                    mapdb.addmap2(tabname, name, gen.miniMap)
                    mapdb.select("*",tabname)

            print("time: ", time.time() - start)

        for number in [str(i) for i in range(1,13)]:
            print (tabname,number)
            genomes = faReader(f"./genomes/{genus}_{species}_chr{number}.fasta",True)
            print("time: ", time.time() - start)
            for gen in genomes:
                name = "".join([gen.genus[0], gen.species[:2], "_", gen.chromosome, "_", gen.accession.replace(".", "")])
                print(name)
                mapdb.addmap2(tabname, name, gen.miniMap, False)
                mapdb.select("*",tabname)
                print("time: ", time.time() - start)

print("time: ", time.time() - start)
speciesFaToMapDb("Capsicum","baccatum", create=True)
print("time: ", time.time() - start)

