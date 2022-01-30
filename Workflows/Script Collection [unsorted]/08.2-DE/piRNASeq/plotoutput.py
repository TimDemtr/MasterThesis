import matplotlib.pyplot as plt
import sys

inputfile = sys.argv[1]

cutoff = 10000000

with open(inputfile) as inputf: 
	inputf = inputf.readlines()


positions = [int(x.split()[0]) for x in inputf if (not x.startswith("#") and int(x.split()[0])<cutoff)]
rates = [float(x.split()[-1]) for x in inputf if (not x.startswith("#") and int(x.split()[0])<cutoff)]

plt.title("HK4me3 marks - Aggregated over all piRNA Clusters")
plt.xlabel("Position")
plt.ylabel("HK4me3 ChiP signal")
plt.plot(positions, rates)
plt.savefig(inputfile.replace(".txt", ".svg"), format="svg", dpi=600)


