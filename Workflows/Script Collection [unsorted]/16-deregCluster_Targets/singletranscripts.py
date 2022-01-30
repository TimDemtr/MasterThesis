import os

header = """##gff-version 2
##source-version rtracklayer 1.54.0
##date 2022-01-26"""

with open("smes_v2_hconf_onlyderegcluster_smests.gtf", "r") as gtfall:
	for line in gtfall:
		if line.startswith("##"):
			continue
		with open("Transcripts_GTFs/"+line.split("\t")[-1].split()[1].strip('";')+".gtf", "w") as outp:
			outp.write(header)
			outp.write("\n")
			outp.write(line)

with open("smes_v2_hconf_onlyderegcluster_smesgs.gtf", "r") as gtfall:
	for line in gtfall:
		if line.startswith("##"):
			continue
		with open("Genes_GTFs/"+line.split("\t")[-1].split()[1].strip('";')+".gtf", "w") as outp:
			outp.write(header)
			outp.write("\n")
			outp.write(line)


