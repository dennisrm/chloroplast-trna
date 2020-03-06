#
"""
Get minimizer maps from genomes and save to sqlite database
Dennis Mulligan
"""
import time
start = time.time()

from minimapdb import *
from soltools import *

print("time: ", time.time() - start)

with ConnectDB("db6.db") as mapdb:

    #genomes = faReader("./genomes/Solanum_cp_genbank.fasta", True)

    print("time: ", time.time() - start)

    #mapdb.create("sly","minimer")
    #for gen in genomes:
    #    if gen.species=="lycopersicum":
    #        name = "".join([gen.genus[0], gen.species[:2], "_", gen.chromosome, "_", gen.accession.replace(".", "")])
    #        print(name)
    #        mapdb.addmap2("sly", name, gen.miniMap)
    #        mapdb.select("*","sly")

    print("time: ", time.time() - start)

    for number in [str(i) for i in range(2,13)]:
        genomes = faReader(f"./genomes/Solanum_lycopersicum_chr{number}.fasta",True)
        print("time: ", time.time() - start)
        for gen in genomes:
            name = "".join([gen.genus[0], gen.species[:2], "_", gen.chromosome, "_", gen.accession.replace(".", "")])
            print(name)
            mapdb.addmap2("sly", name, gen.miniMap, False)
            mapdb.select("*","sly")
        print("time: ", time.time() - start)

print("time: ", time.time() - start)

