with open("smes_v2_hconf_genes_SMESG.bed", "r") as handle: 
	f = handle.readlines()

with open("smes_v2_hconf_genes_SMESG_geneonly.bed", "w") as handle: 
	for line in f: 
		if line.split("\t")[7]=="gene": 
			handle.write(f"{line}")
