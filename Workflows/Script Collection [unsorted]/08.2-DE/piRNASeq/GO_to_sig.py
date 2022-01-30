import os
import json

filepairs = dict()
os.chdir("GO_out_graphs/TSVs")

try:
    os.mkdir("GOtoGene_X_sigGenes")
except:
    pass

os.chdir("GOtoGene")
for txt in [x for x in os.listdir() if x.endswith(".txt")]: 
	with open(txt, "r") as inputtxt: 
		readlines = inputtxt.readlines()
		diggtri = dict()
	for line in readlines:
		if line.startswith("$"):
			current = line.replace("$", "").replace("`", "").strip()
			diggtri[current] = list()
			continue
		elif line=="\n":
			continue
		else:
			diggtri[current] += [x.replace('"', '') for x in line.split()[1:]]
		
		
	json.dump(diggtri, open(f"{txt.strip('.txt')}.json", "w"), indent=4)
	os.remove(txt)

os.chdir("..")

for file in [x for x in os.listdir("GOtoGene") if x.endswith(".json")]:
    print(file)
    with open(f"./sigGenes/{'.'.join(file.split('.')[:-2])}.sigGenes.txt") as siggenes_file:
        siggenes_list = [x.strip() for x in siggenes_file.readlines()]
    GOtoGene_json = json.load(open(f"./GOtoGene/{file}"))
    
    GOtoGene_x_sig = dict()
    
    for GO, genes in GOtoGene_json.items(): 
        GOtoGene_x_sig[GO] = list()
        for siggene in siggenes_list:
            if siggene in genes: 
                GOtoGene_x_sig[GO].append(siggene)
    GOtoGene_x_sig_filtered = dict()
    for GO, genes in GOtoGene_x_sig.items():
        if len(GOtoGene_x_sig.get(GO)): 
            GOtoGene_x_sig_filtered[GO] = genes
          
    json.dump(GOtoGene_x_sig_filtered, open(f"./GOtoGene_X_sigGenes/{'.'.join(file.split('.')[:-2])}.GOtoGene_X_sigGenes.json", "w"), indent=4)


