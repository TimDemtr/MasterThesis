import matplotlib.pyplot as plt
import json
import pandas as pd
import seaborn as sns
import os
import sys

def sumzip(*items):
    return [sum(values) for values in zip(*items)]


if __name__ == "__main__":
    try:
        folder = sys.argv[1]
    except:
        folder = "./"

    os.chdir(folder)

    ecoli_dict = {"UP": {}, "DOWN": {}}
    polyIC_dict = {"UP": {}, "DOWN": {}}

    cond_dict = json.load(open("Analysis/de_unique_transcriptsonly.json", "r"))
    for entry, updown in cond_dict.items():
        for state, inbetween_info in updown.items():
            for inbetween, transcripts in inbetween_info.items():
                locals()[f"{entry}_dict"][state][inbetween] = len(transcripts)

    plt.rcParams["figure.figsize"] = (8, 10)
    plt.style.use("seaborn-deep")
    fig, ax = plt.subplots()

    up_ecoli = pd.DataFrame.from_dict(ecoli_dict.get("UP"), orient="index", columns=["count", ])
    down_ecoli = pd.DataFrame.from_dict(ecoli_dict.get("DOWN"), orient="index", columns=["count", ])
    up_polyIC = pd.DataFrame.from_dict(polyIC_dict.get("UP"), orient="index", columns=["count", ])
    down_polyIC = pd.DataFrame.from_dict(polyIC_dict.get("DOWN"), orient="index", columns=["count", ])

    up_ecoli_plt = up_ecoli.plot(kind="bar", legend=False, color=sns.color_palette("deep"))
    down_ecoli_plt = down_ecoli.plot(kind="bar", legend=False, color=sns.color_palette("deep"))
    up_polyIC_plt = up_polyIC.plot(kind="bar", legend=False, color=sns.color_palette("deep"))
    down_polyIC_plt = down_polyIC.plot(kind="bar", legend=False, color=sns.color_palette("deep"))

    up_ecoli_plt.figure.savefig(f"ecoli_up.svg", dpi=600, format="svg")
    down_ecoli_plt.figure.savefig(f"ecoli_down.svg", dpi=600, format="svg")
    up_polyIC_plt.figure.savefig(f"polyIC_up.svg", dpi=600, format="svg")
    down_polyIC_plt.figure.savefig(f"polyIC_down.svg", dpi=600, format="svg")