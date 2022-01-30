import os
import shutil
import subprocess
import sys
import polars as pl

if __name__ == "__main__": 
	if not sys.argv[1]:
		print("Missing Folder to merge...")
	inputfolder=sys.argv[1].strip("/")
	if inputfolder not in [x for x in os.listdir() if os.path.isdir(x)]:
		print("Not a valid Folder...")

	print(f"Folder selected: {inputfolder}")

	outputfile = f"../../{inputfolder}.merged.gtf"
	os.chdir(inputfolder)
	
	gtf_files = [x for x in os.listdir() if x.endswith("gtf")]

	cat1=subprocess.Popen(["cat",]+gtf_files, stdout=subprocess.PIPE)
	fileout=open(f"{inputfolder}.cat_sort.gtf", "wb")
	sort1=subprocess.run(["sort", "-k1,1", "-k4,4n"], stdin=cat1.stdout, stdout=fileout)
	fileout.close()
	bedtools=subprocess.Popen(["bedtools", "merge", "-c", "6,7,8", "-o", "distinct", "-s", "-i", f"{inputfolder}.cat_sort.gtf",], stdout=subprocess.PIPE)
	fileout=open(outputfile, "wb")
	sort2=subprocess.run(["sort", "-k1,1V", "-k2,2n",], stdin=bedtools.stdout, stdout=fileout)
	fileout.close()
		
	type_info = {"E": "Ecoli_treated_Smed", "P": "PolyIC_treated_Smed", "WT": "WildType_Smed", "uni": "Combined_of_WTEP"}
	## Goal: ["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attribute"
	bed_df = pl.read_csv(outputfile, sep="\t", has_header=False, new_columns=["chr", "start", "end","score", "strand", "phase"])
	bed_df["source"] = pl.Series(["proTRAC" for x in range(0,len(bed_df))])
	sample_type = "uni"
	bed_df["feature"] = pl.Series([f"piRNA_cluster" for x in range(0,len(bed_df))])
	strand = bed_df.to_series(4).to_list()
	type_entry = type_info.get(sample_type)
	bed_df["attribute"] = pl.Series([f'CLUSTER §piRNA_Cluster_{x+1}§; '
									 f'DIRECTIONALITY §mono_{"plus" if bed_df["strand"][x]=="+" else ("minus" if bed_df["strand"][x]=="-" else "bi")}§; '
									 f'ORIGIN §{type_entry}§'
									 for x in range(0,len(bed_df))])
	bed_df_ordered = bed_df[["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attribute"]]
	bed_df_ordered["start"]=bed_df_ordered["start"].apply(lambda x: x+1)
	bed_df_ordered["end"]=bed_df_ordered["end"].apply(lambda x: x+1)
	bed_df_ordered.to_csv(outputfile.replace(".bed", ".gtf"), sep="\t", has_header=False)
	with open(outputfile, "r") as bed_read:
		bed_fix = bed_read.read().replace("§", "\"")
	with open(outputfile, "w") as bed_write:
		bed_write.write(bed_fix)
	print(f"Merged file saved here: {outputfile}")
		

