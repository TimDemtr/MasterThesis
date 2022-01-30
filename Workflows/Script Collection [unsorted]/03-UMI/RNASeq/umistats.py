import os
import json
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
import pprint

pp = pprint.PrettyPrinter(indent=4)

def write_into_dict(file):
    global stats
    print(f"Processing of {file}")
    rep = f"{file.split('_')[-1].strip('.log')}"
    with open(file) as f:
        lines = f.readlines()
        total_reads = lines[-5].split()[-1]
        read_matches = lines[-4].split()[-1]
        no_match = lines[-2].split()[-1]
    if file.startswith("P"):
        stype="PolyIC"
    elif file.startswith("E"):
        stype="Ecoli"
    stats[stype][rep]["total"] = int(total_reads)
    stats[stype][rep]["passed"] = int(read_matches)
    stats[stype][rep]["imp"] = int(no_match)



if __name__ == "__main__":
    skipreadin=False
    if not skipreadin:
        #duplicates = [x for x in os.listdir() if x.endswith("_dups.fq")]
        readin = [x for x in os.listdir() if x.endswith(".log")]
        stats = {"Ecoli": {}, "PolyIC": {}}
        for x in range(1, 4):
            stats["Ecoli"][f"r{x}"] = {"total": 0, "passed": 0, "imp": 0}
            stats["PolyIC"][f"r{x}"] = {"total": 0, "passed": 0, "imp": 0}
        print(stats)
        #for file in duplicates:
        #    write_into_dict("dup", file)
        for file in readin:
            write_into_dict(file)
        json.dump(stats, open("UMI_stats.json", "w"), indent=4)


    stats = json.load(open("UMI_stats.json", "r"))
    pp.pprint(stats)
    for sampletype in stats.keys():
        labels = ["r1", "r2", "r3"]
        barplot_passed = [stats[sampletype][f"r{key}"].get("passed")/stats[sampletype][f"r{key}"].get("total")*100 for key in range(1, 4)]
        barplot_imp = [stats[sampletype][f"r{key}"].get("imp")/stats[sampletype][f"r{key}"].get("total")*100 for key in range(1, 4)]
        print(barplot_passed)
        print(barplot_imp)
        plt.rcParams["figure.figsize"] = (6, 7)
        fig, ax = plt.subplots()
        # ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
        width = 0.5
        ax.bar(labels, barplot_passed, color="#76BBDB", bottom=barplot_imp, width=width, label='Passed')
        ax.bar(labels, barplot_imp, color="#DBAD76", width=width, label='Improper')
        ax.set_ylabel('Percentage')
        ax.set_title(f'{sampletype} - UMI Filtering')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.04),ncol=5)
        plt.ylim((0, 100))
        plt.savefig(f"{sampletype}_stats.svg", dpi=600, format="svg")
