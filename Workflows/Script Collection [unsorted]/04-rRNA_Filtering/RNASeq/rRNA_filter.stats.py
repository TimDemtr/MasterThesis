import os
import json
import pprint
import matplotlib.pyplot as plt
import seaborn as sns

pp = pprint.PrettyPrinter()
stats = {}

lookout_for = """    Total reads =
    Total reads passing E-value threshold =
    Total reads failing E-value threshold =""".splitlines()

for file in [x for x in os.listdir() if x.endswith("log")]:
    name = file.split("_")[3]
    stats[name] = {}

for file in [x for x in os.listdir() if x.endswith("log")]:
    name = file.split("_")[3]
    replicate = file.split("_")[-4]
    stats[name][replicate] = {"total":0, "rRNA":0, "passed":0}
    with open(file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(lookout_for[0]):
                stats[name][replicate]["total"] = int(line.split()[-1])
            if line.startswith(lookout_for[1]):
                stats[name][replicate]["rRNA"] = int(line.split()[-2])
            if line.startswith(lookout_for[2]):
                stats[name][replicate]["passed"] = int(line.split()[-2])

json.dump(stats, open("stats.json", "w"), indent=4)
pp.pprint(stats)
stats = json.load(open("stats.json", "r"))
sns.set_palette("deep")
for sampletype in stats.keys():
    labels = ["r1", "r2", "r3"]
    barplot_rRNA = [stats[sampletype][f"r{key}"].get("rRNA") / stats[sampletype][f"r{key}"].get("total") * 100 for key in range(1, 4)]
    barplot_passed = [stats[sampletype][f"r{key}"].get("passed") / stats[sampletype][f"r{key}"].get("total")*100 for key in range(1, 4)]
    plt.rcParams["figure.figsize"] = (6, 7)
    fig, ax = plt.subplots()
    # ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
    width = 0.5
    ax.bar(labels, barplot_passed,
           width=width, label='Passed')
    ax.bar(labels, barplot_rRNA,width=width, label='rRNA', bottom=barplot_passed)
    ax.set_ylabel('Percentage')
    plt.ylim((0,100))
    ax.set_title(f'{sampletype} - UMI Filtering')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.04), ncol=3)
    plt.savefig(f"{sampletype}_stats.svg", dpi=600, format="svg")
