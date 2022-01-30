import os
import json


if __name__ == "__main__":
    path  = "./GO_out_graphs/TSVs/Plotted_Terms"
    for file in os.listdir(path):
        full_path = os.path.join(path,file)
        GO_terms = list()
        with open(full_path) as GO_terms_plotf:
            GO_terms_plotted = [x.split()[1] for x in GO_terms_plotf.readlines()[1:]]
        for GO_term in GO_terms_plotted:
            GO_terms.append(GO_term)


        print(GO_terms)
        exit()