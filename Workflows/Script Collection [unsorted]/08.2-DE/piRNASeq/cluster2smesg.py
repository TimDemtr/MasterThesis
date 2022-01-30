import os
import json

output_folder = os.getcwd()
try:
	os.mkdir("combined")
	os.mkdir("sense")
	os.mkdir("antisense")
except:
	pass
def cluster2smesg(file: str, folder: str, mapped: str="sense"):
	output_dir = dict()
	os.chdir(folder)
	with open(file, "r") as handle: 
		f = handle.readlines()
	
	for line in f: 
		attr = line.split("\t")[-1].split(";")
		cluster=attr[0].split()[-1].strip('"')
		smesgs=attr[-1].split()[-1].strip('"').split(",")
		
		output_dir[cluster]=smesgs

	json.dump(output_dir,open(output_folder+f"/{mapped}/"+file.replace("gtf","json"), "w"),indent=4)
	os.chdir(output_folder)


#Sense
for file in [x for x in os.listdir("../0-sense_mapped") if x.endswith(".gtf")]:
	cluster2smesg(file, "../0-sense_mapped", "sense")
#Antisense
for file in [x for x in os.listdir("../0-antisense_mapped") if x.endswith(".gtf")]:
	cluster2smesg(file, "../0-antisense_mapped", "antisense")

#combined
for file in [x for x in os.listdir("sense") if x.endswith(".json")]:
	sense = json.load(open(f"./sense/{file}", "r")) 
	antisense = json.load(open(f"./antisense/{file.replace('sense', 'antisense')}", "r"))
	combined = {"sense":sense, "antisense":antisense}
	json.dump(combined,open(output_folder+f"/combined/"+file.replace("sense","combined"), "w"),indent=4)
