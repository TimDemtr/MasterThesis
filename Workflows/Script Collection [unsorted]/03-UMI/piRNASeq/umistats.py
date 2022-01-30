import os
import json
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
import pprint

pp = pprint.PrettyPrinter(indent=4)


def write_into_dict(itype, file):
    global stats
    print(f"Processing {itype} of {file}")
    rep = f"{file.split('_')[-2]}"
    with open(file) as f:
        dup = int(len(f.readlines())/4)
    if file.startswith("P"):
        stype="PolyIC"
    elif file.startswith("E"):
        stype="Ecoli"
    stats[stype][rep][itype] = dup

def calc_total(stats:dict):
    total = {}
    for key in stats.keys():
        total[key] = {"r1":0, "r2":0, "r3":0}
        print(stats[key])
        for n in range(1,4):
            total[key][f"r{n}"]=(sum(stats[key][f"r{n}"].values()))
    return total

if __name__ == "__main__":
    skipreadin=True
    if not skipreadin:
        duplicates = [x for x in os.listdir() if x.endswith("_dups.fq")]
        imp = [x for x in os.listdir() if x.endswith("_impropUMI.fq")]
        single = [x for x in os.listdir() if x.endswith("_single.fq")]
        stats = {"Ecoli": {}, "PolyIC": {}}
        for x in range(1, 4):
            stats["Ecoli"][f"r{x}"] = {"dup": 0, "single": 0, "imp": 0}
            stats["PolyIC"][f"r{x}"] = {"dup": 0, "single": 0, "imp": 0}
        print(stats)
        for file in duplicates:
            write_into_dict("dup", file)
        for file in imp:
            write_into_dict("imp", file)
        for file in single:
            write_into_dict("single", file)
        json.dump(stats, open("UMI_stats.json", "w"), indent=4)


    stats = json.load(open("UMI_stats.json", "r"))
    pp.pprint(stats)
    total = calc_total(stats)
    for sampletype in stats.keys():
        labels = ["r1", "r2", "r3"]
        barplot_dup = [stats[sampletype][f"r{key}"].get("dup")/total[sampletype].get(f"r{key}")*100 for key in range(1, 4)]
        barplot_single = [stats[sampletype][f"r{key}"].get("single")/total[sampletype].get(f"r{key}")*100 for key in range(1, 4)]
        barplot_imp = [stats[sampletype][f"r{key}"].get("imp")/total[sampletype].get(f"r{key}")*100 for key in range(1, 4)]
        print(barplot_dup)
        print(barplot_single)
        print(barplot_imp)
        plt.rcParams["figure.figsize"] = (6, 7)
        fig, ax = plt.subplots()
        # ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
        width = 0.5
        ax.bar(labels, barplot_single, color="#76BBDB", bottom=[i+j for i,j in zip(barplot_imp, barplot_dup)], width=width, label='Single')
        ax.bar(labels, barplot_imp, color="#DBAD76", bottom=barplot_dup, width=width, label='Improper')
        ax.bar(labels, barplot_dup, color="#8F6431", width=width, label='Duplicates')
        ax.set_ylabel('Percentage')
        ax.set_title(f'{sampletype} - UMI Filtering')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.04),ncol=5)
        plt.ylim((0, 100))
        plt.savefig(f"{sampletype}_stats.svg", dpi=600, format="svg")
