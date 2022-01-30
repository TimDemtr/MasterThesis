import os
import json
import pprint
import matplotlib.pyplot as plt
import seaborn as sns

pp = pprint.PrettyPrinter()
stats = {}

for file in [x for x in os.listdir() if x.endswith("log")]:
    name = file.split("_")[3]
    stats[name] = {}

for file in [x for x in os.listdir() if x.endswith("log")]:
    name = file.split("_")[3]
    replicate = file.split("_")[-3]
    stats[name][replicate] = {"total":0, "tRNA":0, "passed":0}
    print(stats)
    with open(file, "r") as f:
        lines = f.readlines()
        stats[name][replicate]["total"] = int(lines[0].split()[-1])
        stats[name][replicate]["tRNA"] = int(lines[1].split()[-2])
        stats[name][replicate]["passed"] = int(lines[2].split()[-2])

json.dump(stats, open("stats.json", "w"), indent=4)

stats = json.load(open("stats.json", "r"))
pp.pprint(stats)
sns.set_palette("deep")
for sampletype in stats.keys():
    labels = ["r1", "r2", "r3"]
    barplot_tRNA = [stats[sampletype][f"r{key}"].get("tRNA") / stats[sampletype][f"r{key}"].get("total")*100 for key in range(1, 4)]
    barplot_passed = [stats[sampletype][f"r{key}"].get("passed") / stats[sampletype][f"r{key}"].get("total")*100 for key in range(1, 4)]
    plt.rcParams["figure.figsize"] = (6, 7)
    fig, ax = plt.subplots()
    # ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
    width = 0.5
    ax.bar(labels, barplot_passed, bottom=barplot_tRNA,
           width=width, label='Passed')
    ax.bar(labels, barplot_tRNA, width=width, label='tRNA')
    ax.set_ylabel('Percentage')
    ax.set_title(f'{sampletype} - UMI Filtering')
    plt.ylim((0,100))
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.04),ncol=5)
    plt.savefig(f"{sampletype}_stats.svg", dpi=600, format="svg")
