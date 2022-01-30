import os
import polars as pl
import pandas as pd

type_info = {"E": "Ecoli_treated_Smed", "P": "PolyIC_treated_Smed", "WT": "WildType_Smed", "uni": "Combined_of_WTEP"}
## Goal: ["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attribute"]
for bed in [x for x in sorted(os.listdir("merge_bed")) if x.endswith(".bed")]:
    bed_df = pl.read_csv("./merge_bed/"+bed, sep="\t", has_header=False, new_columns=["chr", "start", "end",
                                                                       "score", "strand", "phase"])
    bed_df["source"] = pl.Series(["proTRAC" for x in range(0,len(bed_df))])
    sample_type = bed.split("_")[1]
    bed_df["feature"] = pl.Series([f"piRNA_cluster" for x in range(0,len(bed_df))])
    strand = bed_df.to_series(4).to_list()
    type_entry = type_info.get(sample_type)
    bed_df["attribute"] = pl.Series([f'CLUSTER §piRNA_Cluster_{x+1}§; '
                                     f'DIRECTIONALITY §mono_{"plus" if bed_df["strand"][x]=="+" else ("minus" if bed_df["strand"][x]=="-" else "bi")}§; '
                                     f'ORIGIN §{type_entry}§'
                                     for x in range(0,len(bed_df))])
    bed_df_ordered = bed_df[["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attribute"]]
    #bed_pdf = bed_df_ordered.to_pandas()
    #bed_pdf = bed_pdf.sort_values(["chr", "start"], ascending=(False, False))
    #bed_pdf.to_csv(bed.replace(".bed", "_fixed.bed"), sep="\t", has_header=False)
    bed_df_ordered.to_csv(bed.replace(".bed", ".gtf"), sep="\t", has_header=False)
    with open(bed.replace(".bed", ".gtf"), "r") as bed_read:
        bed_fix = bed_read.read().replace("§", "\"")
    with open(bed.replace(".bed", ".gtf"), "w") as bed_write:
        bed_write.write(bed_fix)
