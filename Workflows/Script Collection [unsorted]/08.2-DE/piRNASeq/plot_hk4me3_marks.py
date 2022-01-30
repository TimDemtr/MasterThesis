import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(style="whitegrid", palette="deep")

RNASeq = pd.read_csv(open("K4me3_x_RNASeq_1000-13000-0.txt"), sep="\t", names=["position", "signal"])
RNASeq["Data Set"] = "RNASeq"
hconfRNASeq = pd.read_csv(open("K4me3_x_hconfRNASeq_1000-13000-0.txt"), sep="\t",names=["position", "signal"])
hconfRNASeq["Data Set"] = "RNASeq & Transcripts"
plt.rcParams['figure.figsize'] = 11, 7
mixed = pd.concat([RNASeq, hconfRNASeq]).reset_index(drop=True)
# Plot the responses for different events and regions
ab = sns.lineplot(x="position", y="signal", hue="Data Set",
             data=mixed,legend="brief")
plt.xlabel("Position relative to Cluster Starting Site")
plt.ylabel("HK4me3 ChiP signal")
plt.title("HK4me3 marks - Aggregated over all piRNA Clusters")
plt.xlim(left=-1000, right=13000)
plt.ylim(bottom=0)
plt.xticks([-1000, 0]+list(range(1000, 13001, 2000)))
plt.legend(title="Clusters intersected with...")
l1 = ab.lines[0]
l2 = ab.lines[1]
x1 = l1.get_xydata()[:,0]
y1 = l1.get_xydata()[:,1]
x2 = l2.get_xydata()[:,0]
y2 = l2.get_xydata()[:,1]
ab.fill_between(x1,y1, alpha=0.2)
ab.fill_between(x2,y2, alpha=0.2)

plt.savefig("HK4me3_marks_aggr_piRNA Clusters.svg", dpi=600, )