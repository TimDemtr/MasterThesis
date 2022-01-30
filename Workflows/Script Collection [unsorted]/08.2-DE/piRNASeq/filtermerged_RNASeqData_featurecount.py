import sys

try:
	filtercount = int(sys.argv[1])
except:
	filtercount = 50 

with open("./RNA_SEQ_BEDs/merged_RNASeqData.bed", "r") as inputfile: 
	lines = inputfile.readlines()

with open(f"merged_RNASeqData_{filtercount}features.bed", "w") as outputfile:
	for line in lines: 
		if int(line.split("\t")[3]) >= filtercount:
			outputfile.write(line)


