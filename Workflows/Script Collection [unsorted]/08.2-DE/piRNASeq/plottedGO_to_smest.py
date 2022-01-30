import os
import json

if __name__ == "__main__":
    path  = "./GO_out_graphs/TSVs/Plotted_Terms"
    for file in os.listdir(path):
        go2sig = "./GO_out_graphs/TSVs/GOtoGene_X_sigGenes/"+'.'.join(file.split('.')[:-3])+".GOtoGene_X_sigGenes.json"
        go2sig = json.load(open(go2sig))
        full_path = os.path.join(path,file)
        GO_terms = dict()
        with open(full_path) as GO_terms_plotf:
            GO_terms_plotted = [x.split()[0] for x in GO_terms_plotf.readlines()[1:]]
        for GO_term in GO_terms_plotted:
            GO_terms[GO_term] = go2sig.get(GO_term)

        json.dump(GO_terms, open(full_path.replace(".tsv", ".json"), "w"), indent=4)