import os
import polars as pl
import pandas as pd
"""
try:
    os.mkdir("origin_noted")
    os.mkdir("origin_noted/pdens0.03")
except:
    pass
gtfs = sorted([x for x in os.listdir("pdens0.03") if x.endswith("gtf")])

for gtf in gtfs:
    gtf_name = "_".join(gtf.split("_")[1:]).replace(".gtf", "")
    out = pl.read_csv("pdens0.03/"+gtf, sep="\t", has_header=False, new_columns=["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attribute"])
    source_replacement = pl.Series([gtf_name for x in range(0,len(out))])
    out["attribute"] = out["attribute"].apply(lambda x: ";".join(x.split(";")[:2]))
    out.replace("source", source_replacement)
    out.to_csv(f"./origin_noted/pdens0.03/pdens003_{gtf_name}_OGnoted.gtf", sep="\t", has_header=False)
"""
## Don't... all the information is lost with merging anyways.

## Bedtools intersect example:
## How to merge the cluster files to "one" cluster, if overlapping?

# Firt concat them, then fix the annotation (02-Merge/fix_merged_annotation.py) and then you are good to go.
#cat pdens003_proTRACluster_E1_r1.gtf pdens003_proTRACluster_E3_r2.gtf pdens003_proTRACluster_E4_r3.gtf | sort -k1,1 -k4,4n | bedtools merge -c 6,7,8 -o distinct -s -i stdin | sort -k1,1V -k2,2n > /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08-piRNAClusters/02-Merge/proTRACluster_E_merged.bed
#cat pdens003_proTRACluster_P1_r1.gtf pdens003_proTRACluster_P2_r2.gtf pdens003_proTRACluster_P3_r3.gtf | sort -k1,1 -k4,4n | bedtools merge -c 6,7,8 -o distinct -s -i stdin | sort -k1,1V -k2,2n > /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08-piRNAClusters/02-Merge/proTRACluster_P_merged.bed
#cat pdens003_proTRACluster_WT_r1.gtf pdens003_proTRACluster_WT_r2.gtf pdens003_proTRACluster_WT_r3.gtf | sort -k1,1 -k4,4n | bedtools merge -c 6,7,8 -o distinct -s -i stdin | sort -k1,1V -k2,2n > /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08-piRNAClusters/02-Merge/proTRACluster_WT_merged.bed
#cat pdens003_proTRACluster_E1_r1.gtf pdens003_proTRACluster_E3_r2.gtf pdens003_proTRACluster_E4_r3.gtf pdens003_proTRACluster_P1_r1.gtf pdens003_proTRACluster_P2_r2.gtf pdens003_proTRACluster_P3_r3.gtf pdens003_proTRACluster_WT_r1.gtf pdens003_proTRACluster_WT_r2.gtf pdens003_proTRACluster_WT_r3.gtf | sort -k1,1 -k4,4n | bedtools merge -c 6,7,8 -o distinct -s -i stdin | sort -k1,1V -k2,2n > /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08-piRNAClusters/02-Merge/proTRACluster_uni_merged.bed
