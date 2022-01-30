from multiprocessing import Process
import os

def weighted(bed2):
    with open(bed2) as bed2_in, open(bed2.replace(".bed2", "_weighted.bed2"), "w") as bed2_out:
        for line in bed2_in:
           readcount = int(line.split("\t")[3])/int(line.split("\t")[4])
           line = line.strip("\n")
           bed2_out.write(f"{line}\t{readcount}\n")

bed2s = [x for x in os.listdir() if x.endswith("mapped.bed2")]
processes = list()

for bed2 in bed2s:
    p = Process(target=weighted, args=(bed2,))
    processes.append(p)
    p.start()

for p in processes:
    p.join()

