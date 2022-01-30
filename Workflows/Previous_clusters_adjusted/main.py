import pandas as pd
from intermine.webservice import Service
import json
import matplotlib.pyplot as plt

# def gene_output(search_query):
#     service = Service("https://planmine.mpibpc.mpg.de:443/planmine/service")
#     query = service.new_query("Gene")
#     query.add_view("sequence.length", "primaryIdentifier")
#
#     query.add_constraint("Gene", "LOOKUP", search_query, "", code = "A")
#
#     return query
#
#
# traf_genes = pd.read_csv(open("traf_by_seb"), sep=" ", names=["Gene", "Blast", "mean x1 weighted counts", "mean xins weighted counts", "mean epidermis weighted counts",
#                                                                "X1 [tpm]", "Xins [tpm]", "Epidermis [tpm]"])
#
# traf_genes_weighted_counts = traf_genes[["Gene", "mean epidermis weighted counts"]].to_dict(orient="records")
#
# for traf_gene in traf_genes_weighted_counts:
#     Gene = traf_gene.get("Gene")
#     output = gene_output(Gene)
#     weighted_counts = traf_gene.get("mean epidermis weighted counts")
#     for row in output.rows():
#         length = row["length"]
#     traf_gene["adjusted mean epidermis weighted counts"] = weighted_counts / length
#     print(f"{Gene}\t{length}\t{weighted_counts}\t{weighted_counts/length}")
#     traf_gene["length"] = length/1000
#
# json.dump(traf_genes_weighted_counts, open("traf_by_seb_adjusted.json", "w"), indent=4)
#

traf_genes_weighted_counts = json.load(open("traf_by_seb_adjusted.json", "r"))

fig, ax = plt.subplots()
plt.xticks(rotation=90)
values = list()
values_nonadjusted = list()

width = 0.3

for gene in traf_genes_weighted_counts:
    if gene.get("adjusted mean epidermis weighted counts") > 0.1:
        continue
    values.append(gene.get("adjusted mean epidermis weighted counts"))
    values_nonadjusted.append(gene.get("mean epidermis weighted counts"))

labels = [x for x in range(0,len(values))]
ax.bar(labels, values, color="#76BBDB", label='Genes')
plt.savefig(f"adjusted.svg", dpi=600, format="svg")

fig, ax = plt.subplots()

ax.bar(labels, values_nonadjusted, color="#76BBDB", label='Genes')
plt.savefig(f"non_adjusted.svg", dpi=600, format="svg")