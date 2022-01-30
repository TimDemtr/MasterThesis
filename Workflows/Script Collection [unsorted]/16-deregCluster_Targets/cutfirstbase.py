### Little script to get rid of the first base... seed region defined as 2-8, but we want to align only in the seed region / with 3 mismatches in total

import os

for folder in [x for x in os.listdir() if os.path.isdir(x)]:
	with open(f"{folder}/combined.insert") as inp, open(f"{folder}/combined_1stcut.insert", "w") as outp: 
		for line in inp:
			outp.write(line[1:])


## Command used afterwards for mapping and filtering: 
"""
for d in */; 
do bowtie -r -p 8 -n 0 -l 7 /Work2/Planarian_genome_2021/Planarian_genome/Bowtie-Index_Smed_g4/Bowtie-Index_Smed_g4 $d/combined_1stcut.insert -S $d/combined_1stcut_mapped.sam 2> $d/combined_mapped.log;
samtools view -bS -@ 7 $d/combined_1stcut_mapped.sam | samtools sort -@ 7 - > $d/combined_1stcut_mapped.bam;
rm $d/combined_1stcut_mapped.sam; bedtools bamtobed -i $d/combined_1stcut_mapped.bam > $d/combined_1stcut_mapped.bed;
piPipes_insertBed_to_bed2 $d/combined_1stcut.insert $d/combined_1stcut_mapped.bed > $d/combined_1stcut_mapped.bed2; 
bedtools intersect -a $d/combined_1stcut_mapped.bed2 -b $d/${d%*/}.gtf -v > $d/combined_1stcut_mapped_originfiltered.bed2; 
awk '{print $0"\t"$4/$5}' $d/combined_1stcut_mapped_originfiltered.bed2 > $d/combined_1stcut_mapped_originfiltered_weigthedCounts.bed2;
done
"""
