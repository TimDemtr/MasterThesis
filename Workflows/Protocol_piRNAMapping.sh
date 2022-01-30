#/bin/bash/
cd /data1/Andreas/Planaria/Glyco_Pathway/piRNA_Pops

###1)Raw Data
cd 01-Raw
mkdir fastqc

#for file in *.fq.gz; do
#fastqc $file -o ./fastqc/
#done

###2)Adapter-Trimming
mkdir ../02-AdapterTrimming

#for file in X1*; do
#output=../02-AdapterTrimming/${file%.fq.gz}_cutadapt.fq.gz
#Log=../02-AdapterTrimming/${file%.fq.gz}.Log
#cutadapt -j 6 -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 40 $file 1> $output 2> $Log
#done

#Epidermal Samples have UMIs => different Adapter Trimming necessary
for file in Epi_pi*; do
output=../02-AdapterTrimming/${file%.fq.gz}_cutadapt.fq
Log=../02-AdapterTrimming/${file%.fq.gz}.Log
cutadapt -j 6 -a TGGAATTCTCGGGTGCCAAGG -m 42 -M 64 --discard-untrimmed $file 1> $output 2> $Log
done

cd ../02-AdapterTrimming
mkdir fastqc


for file in *.fq; do
fastqc $file -o ./fastqc/
done


###3)UMI-Filtering
# Epidermal Samples have UMIs -> need to be filtered
mkdir ../../03-UMI
mkdir ../../03-UMI/piRNASeq

for file in *_PolyIC_pi_*.fq; do
single=../../03-UMI/piRNASeq/${file%_cutadapt.fq}_single.fq
dup=../../03-UMI/piRNASeq/${file%_cutadapt.fq}_dups.fq
improperumi=../../03-UMI/piRNASeq/${file%_cutadapt.fq}_impropUMI.fq
log=../../03-UMI/piRNASeq/${file%_cutadapt.fq}.log
umitools reformat_sra_fastq -i $file -o $single -d $dup --reads-with-improper-umi $improperumi -e 1 -v -5 NNNCGANNNTACNNN,NNNATCNNNAGTNNN -3 NNNGTCNNN > $log
done
# -e sets allowed errors in UMIs (1); UMI sequences are added via -5 & -3 arguments

gzip *.fq
cd ../03-UMI/piRNASeq
mkdir fastqc

for file in *.fq; do
fastqc $file -o ./fastqc &
done

###6)Mapping
#generate Insert files; needed for bed2 file generation later on
mkdir ../../06-Mapping/piRNA

for file in *single.fq; do
piPipes_fastq_to_insert $file ../06-Mapping/${file%.fq}.insert
done

gzip *.fq &&
gzip *.sam

#build bowtie index of planarian genome
cd /data1/Andreas/Genome/Planarian_genome

bowtie-build dd_Smes_g4_genome.fasta ./Bowtie-Index_Smed_g4

#bowtie mapping using insert files vs genome + reformatting to bam file
cd /data1/Andreas/Planaria/Glyco_Pathway/piRNA_Pops/06-Mapping

for file in *.insert; do
bowtie -r -p 7 -v 1 -a --best --strata /Work2/Planarian_genome/Bowtie-Index_Smed_g4/Bowtie-Index_Smed_g4 $file -S ${file%.insert}_mapped.sam 2>> ${file%.insert}.log
samtools view -bS -@ 7 ${file%.insert}_mapped.sam | samtools sort -@ 7 - > ${file%.insert}_mapped.bam
rm ${file%.insert}_mapped.sam
done

###7)IGV files for Coverage
# generated for visualization purposes; can be skipped
mkdir ../07-IGV

for file in *.bam; do
samtools index -b $file ${file%.bam}.bai
done

for file in *.bam; do
bamCoverage -p 8 -b $file -o ../07-IGV/${file%.bam}.bigwig
done

for file in *.bam; do
bamCoverage -p 8 -b $file --normalizeUsing RPKM --filterRNAstrand reverse -o /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/07-IGV/piRNA/RPKM_normalized/${file%.bam}.plus.bw
bamCoverage -p 8 -b $file --normalizeUsing RPKM --filterRNAstrand forward -o /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/07-IGV/piRNA/RPKM_normalized/${file%.bam}.minus.bw
done

#Generating Bedgraph files for strandedness from BAM-Files (have to be sorted)
#-ibam: Input=BAM-File; -bg=report depth; -strand=define strand specifity
for file in *.bam; do
output_forward=../07-IGV/${file%.bam}_forward.bedgraph
output_reverse=../07-IGV/${file%.bam}_reverse.bedgraph
bedtools genomecov -ibam $file -bg -strand + > $output_forward
bedtools genomecov -ibam $file -bg -strand - > $output_reverse
done

#create negative coverage for reverse bedgraph (IGV visualizazion purposes)
cd ../07-IGV
for file in *reverse.bedgraph; do
awk '{print $1"\t"$2"\t"$3"\t-"$4;}' $file > ${file%.bedgraph}_IGV.bedgraph
rm $file
done


###8)Bed-File Generation
# creating bed files from mapped bams; bed is then transformed to bed2 format (including read sequence)
mkdir ../08-Bed
cd ../06-Mapping

for file in *.bam; do
bedtools bamtobed -i $file > ../../08-Bed/${file%bam}bed
piPipes_insertBed_to_bed2 ${file%_mapped.bam}.insert ../../08-piRNAClusters/${file%bam}bed > ../../08-piRNAClusters/${file%bam}bed2
done

## I will use proTRAC with adjusted parameters... 
# proTRAC takes an input map file created with sRNAmapper.pl...
# > perl sRNAmapper.pl -input piRNAs.fasta -genome genome.fasta -alignments best
# However, I want to keep it streamlined with my BAM file... So I am building the map output myself based on the BAM/BED(2) files..

# From proTRAC documentation: https://www.smallrnagroup.uni-mainz.de/software/proTRAC_documentation_v2.4.4.pdf
##With the above command, sRNAmapper will create a map file that is named piRNAs.fasta.map (you can give it a 
##different name using the option -output othername.map) and that will have the following format: 
##Chr1  238  TGTTACGGCTAGCTCAGTACGGC  23  TGTTACGGCTAGCTCAGTACGAA  2  + 
##Chr1  291  TACGCCAGCTCGACTCGCCTGTGCA  23  TACGCCAGCTCGACTCGCCTGTGCA  0  -

mapping_folder = "/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/06-Mapping/piRNASeq/"
cd ../../08-piRNAClusters/

for file in *.bed; do 
bedtools getfasta -fi /Work2/Planarian_genome_2021/Planarian_genome/dd_Smes_g4_genome.fasta -bed $file -tab -fo ../getfasta/${file%.bed}_getfasta.bed -s &
done


# proTRAC2 with build MAP -> See script bed2bedinsert_to_map.py -> It also accounts for the read counts in the insert file. 
for file in ../map/*.map; do
echo $file
perl proTRAC_2.4.4.pl -map $file -genome /Work2/Planarian_genome_2021/Planarian_genome/dd_Smes_g4_genome.fasta -geneset /Work2/Planarian_genome_2021/Planarian_genome/smes_v2_hconf_FullAnnotation.gff3 -pimin 25 -pimax 35 
done

## Tried out many different settings, yielding 50-1000 clusters, depending on the settings... pdens 0.03 with a clustersize of 1000 (ref Iana) did yield around 500, going with that.
## Reanalyzed Sebastian's WT as well - Good here: Other worms (Beth), never analyzed with proTRAC (big difference between CLIP and FULL piRNA population...
## Epidermal piRNA give a completly new light on the whole problem, as we can account only two loci - should have used proTRAC to compare his cluster "prediction" pipeline... 
## But compared only to 270 piRNA clusters of Iana (pdens 0.09 - Really loose...) 

## In the end, proTRAC does similiar stuff than what Sebastian's pipeline did, but basically better - rolling window, density measurment, etc. 
### Guess, it was more a "wanting to create a software" than actually needing a new one. "Das Rad muss nicht neu erfunden werden".

## Also, using SAM as input, yielded around the same number of clusters - so this validates my map files.. Only problem: Sometimes the left most position is +1, sometimes not. So there CAN BE a slight offset of the left most start of a cluster. But that should not matter. 
### Look into this in the future. Additionally, I want to write a pipeline utilizing Histone marks, expression levels of RNASeq etc. automatically to actually validate all the found clusters. 

######################################################################
######################################################################

## Multiple GTFs from different replicates. What to do? 
## -> Create a GTF with all overlapping regions of a replicate. 
## -> Create a GTF with ALL found clusters -> DESeq Input later to analyse differential expression of the found clusters. 
## -> Intersect piRNA Clusters with SMESG genes -> we only want to look at the genes that are "annotated/predicted"
## -> HMM analysis of the piRNA clusters found -> Look for HMM domains related to our up/down regulated genes :) 


cd /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08-piRNAClusters/01-Cluster_GTF_Collection
for file in *.bed2; do
bedtools intersect -split -wo -f 0.5 -a $file -b /data1/Andreas/Genome/Planarian_genome/smes_v2_hconf_genes_SMESG.bed > ../11-Intersect/${file%_weightCounts.bed2}_GeneIntersect.bed
bedtools intersect -split -wo -s -f 1 -a $file -b /data1/Andreas/Genome/Planarian_genome/smes_v2_hconf_genes_SMESG.bed > ../11-Intersect/${file%_weightCounts.bed2}_StrictSenseIntersect.bed
bedtools intersect -split -wa -c -s -F 1 -a /data1/Andreas/Genome/Planarian_genome/smes_v2_hconf_genes_SMESG.bed -b $file > ../11-Intersect/${file%_weightCounts.bed2}_StrictSenseCounts
done


######################################################################
######################################################################
#bed2-Tips:
# to retrieve unique mappers
#awk '$5==1' $BED2

# to retrieve multi-mappers
#awk '$5>1' $BED2

# to count unique-mappers reads
#awk '{if ($5==1) c+=$4 }END{print c}' $BED2

# to count unique + multi-mappers reads, with multi-mappers been partitioned to the number of times
# it can be mapped
#awk '{c+=$4/$5}END{print c}' $BED2
######################################################################
######################################################################

###9)Length Distribution & Cutoff
cd ../08-Bed
mkdir ../09-Length_Cutoff

for file in *.bed2; do
awk '{print length($7);}' $file | sort -n | uniq -c > ../09-Length_Cutoff/${file%.bed2}_Length.txt
done

for file in *.bed2; do
awk  'length($7) > 25' $file > ../09-Length_Cutoff/${file%.bed2}_min25nt.bed2
done

gzip *.bed*

###10) Calculate weighted counts
#create weighted counts by ReadCount/MappingPositions
cd ../09-Length_Cutoff
mkdir ../10-WeightCounts

for file in *25nt.bed2; do
awk '{print $0"\t"$4/$5}' $file > ../10-WeightCounts/${file%_mapped_min25nt.bed2}_weightCounts.bed2
done

gzip *.bed2

###11) Intersect vs Gene annotations
mkdir ../11-Intersect
cd ../10-WeightCounts


#LINE= /data1/Andreas/Genome/Planarian_genome/Repeat-Maskers/RepeatMaskerLINE.gtf
#DNA= /data1/Andreas/Genome/Planarian_genome/Repeat-Maskers/RepeatMaskerDNA.gtf
#LTR= /data1/Andreas/Genome/Planarian_genome/Repeat-Maskers/RepeatMaskerLTR.gtf
#Unk= /data1/Andreas/Genome/Planarian_genome/Repeat-Maskers/RepeatMaskerUnk.gtf
#No_annot=
#Cluster=
#Exon= /data1/Andreas/Genome/Planarian_genome/smes_v2_hconf_genes_SMESG.gtf

for file in *.bed2; do
bedtools intersect -split -wo -f 0.5 -a $file -b /data1/Andreas/Genome/Planarian_genome/smes_v2_hconf_genes_SMESG.bed > ../11-Intersect/${file%_weightCounts.bed2}_GeneIntersect.bed
bedtools intersect -split -wo -s -f 1 -a $file -b /data1/Andreas/Genome/Planarian_genome/smes_v2_hconf_genes_SMESG.bed > ../11-Intersect/${file%_weightCounts.bed2}_StrictSenseIntersect.bed
bedtools intersect -split -wa -c -s -F 1 -a /data1/Andreas/Genome/Planarian_genome/smes_v2_hconf_genes_SMESG.bed -b $file > ../11-Intersect/${file%_weightCounts.bed2}_StrictSenseCounts
done

#-wa/wo decides which original entries are reported; -c counts overlaps of B vs A, -s requires same strandedness; -F sets fraction of B that has to overlap feature in A


##########
#Intersect bed fields
#1-3: Read Coordinates
#4: Read Count
#5: No. of Mapping positions per Read
#6: Strand
#7: Sequence of Read
#8: Weighted Count
#9-11: Gene Coordinates
#12: ID
#13:
#14: Gene Strand
#15: irrelevant (source of gene prediction)
#16: Region Type (exon/CDS/gene/etc)
#17: 
#18: Transcript IDs (is split into $17-$20)
#19: No of bp overlaps

##RPM normalization
#filter intersect files for gene entries only
cd ../11-Intersect
mkdir Sense
mkdir AntiSense

for file in *GeneIntersect.bed; do
awk '{if($5==1 && $11=="exon" && $6==$15) print $0}' $file > Sense/${file%Intersect.bed}MappingUniq.bed
done

for file in *GeneIntersect.bed; do
awk '{if($5==1 && $11=="exon" && $6!=$15) print $0}' $file > AntiSense/${file%Intersect.bed}MappingUniq.bed
done

#Use file of uniq piRNA mappers to genes for RPM calculation
cd Sense
for file in *Uniq.bed; do
awk '{a[$18]+=$4}END{for (b in a) {printf("%s\t%d\n",b,a[b])}}' $file > ${file%.bed}_Counts.csv 
done

#$18 = $20 when working on Intersects vs gtf file

cd ../AntiSense
for file in *Uniq.bed; do
awk '{a[$18]+=$4}END{for (b in a) {printf("%s\t%d\n",b,a[b])}}' $file > ${file%.bed}_Counts.csv 
done







