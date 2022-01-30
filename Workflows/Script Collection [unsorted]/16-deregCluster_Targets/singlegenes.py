import os

smest2smesg = dict()
with open("/Work2/Planarian_genome_2021/Planarian_genome/SMEST_to_SMESG.tsv") as inpu:
    for line in inpu:
        smest2smesg[line.split()[0]] = line.split()[1].strip()
listoftranscripts = """SMEST013545001.1
SMEST022612001.1
SMEST022691001.1
SMEST032715001.1
SMEST080845001.1
SMEST048748001.1
SMEST059031001.1
SMEST072041001.1
SMEST077398001.1
SMEST018801001.1
SMEST040961001.1
SMEST014843005.1"""

listofgenes = set([smest2smesg[x] for x in listoftranscripts.splitlines()])

with open("smes_v2_hconf_FullAnnotation.gtf") as gtf, open("smes_v2_hconf_onlyderegcluster_smesgs.gtf", "w") as gtfout:
    for line in gtf:
        if line.startswith("##"):
            continue
        if line.split("\t")[-1].split()[1].strip('";') in listofgenes:
            gtfout.write(line)


