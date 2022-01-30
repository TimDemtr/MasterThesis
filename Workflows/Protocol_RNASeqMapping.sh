#/bin/bash/
cd /data1/Andreas/Planaria/Glyco_Pathway/piRNA_Pops

###1)Raw Data
#cd 01-Raw
#mkdir fastqc

#for file in *.fq.gz; do
#fastqc $file -o ./fastqc/
#done

###2)Adapter-Trimming
mkdir ../02-AdapterTrimming


#Epidermal Samples have UMIs => different Adapter Trimming necessary
#36, because 9 UMI (9+9 -> read through = empty adapter), 18 UMI + 18 minimal length to map properly (source Andi)
for file in *; do
output=../../02-AdapterTrimming/RNASeq/${file%.fq.gz}_cutadapt.fq
Log=../../02-AdapterTrimming/RNASeq/${file%.fq.gz}.Log
cutadapt -j 6 -a AGATCGGAAGAGCACACGTCT -m 36 $file 1> $output 2> $Log
done

cd ../02-AdapterTrimming/RNASeq/
mkdir fastqc


for file in *.fq; do
fastqc $file -o ./fastqc/
done


###3)UMI-Filtering
# Epidermal Samples have UMIs -> need to be filtered
mkdir ../../03-UMI/RNASeq/


for file in *.fq; do
out=../../03-UMI/RNASeq/${file%_cutadapt.fq}_umitools_out.fq
log=../../03-UMI/RNASeq/${file%_cutadapt.fq}.log
improper=../../03-UMI/RNASeq/${file%_cutadapt.fq}_improp.fq
umi_tools extract --extract-method=regex -p "^(?P<umi_1>.{5}((ATC)|(GGG)|(TCA))){s<=1}" -I $file -S $out --filtered-out=$improper -L $log &
done
# In the above: Maybe use REGEX to only look for UMIs with the UMI Locators... so ^(?P<umi_1>.{5}((ATC)|(GGG)|(TCA))){s<=1} [allowing for one mismatch) 

gzip Epi*.fq
cd ../03-UMI
mkdir fastqc

for file in *umitools_out.fq; do
fastqc $file -o ./fastqc/ &
done

###4)rRNA filtering
#mapping against rRNA databases using SortMeRNA
mkdir ../../04-rRNA_Filter/
mkdir ../../04-rRNA_Filter/RNASeq/

for file in *out.fq; do
aligned=../../04-rRNA_Filter/RNASeq/${file%_single.fq}_rRNAs
unaligned=../../04-rRNA_Filter/RNASeq/${file%_single.fq}_unaligned
sortmerna --ref /Work2/Planarian_genome/rRNA_silva_LSU-SSU/planarian_rRNA_full_sequences.fa,/Work2/Planarian_genome/rRNA_silva_LSU-SSU/planarian_rRNA_full_sequences.idx --reads $file --fastx -a 7 --aligned $aligned --other $unaligned -v --log
done
gzip *.fq

cd ../04-rRNA_Filter
mkdir fastqc

for file in *.fq; do
fastqc $file -o ./fastqc
done


###5)tRNA filtering
build index for planarian tRNAs
cd /data1/Andreas/Genome/Planarian_genome
mkdir tRNA_BowtieIndex
bowtie-build -f insilico_Smed_tRNA_sequences.fa ./tRNA_BowtieIndex/Smed_tRNA

#map vs tRNAs; unmapped will be used further
#cd /data1/Andreas/Planaria/Glyco_Pathway/piRNA_Pops/04-rRNA_Filter
mkdir ../../05-tRNA_Filter
mkdir ../../05-tRNA_Filter/RNASeq

for file in *unaligned.fq; do
tRNA=../../05-tRNA_Filter/RNASeq/${file%_unaligned.fq}_tRNA_hits.sam
unaligned=../../05-tRNA_Filter/RNASeq/${file%_unaligned.fq}_RNASeq.fq
log=../../05-tRNA_Filter/RNASeq/${file%_unaligned.fq}.log
bowtie -q -p 7 -v 1 -a /Work2/Planarian_genome/tRNA/Bowtie_indices/Smed_tRNA $file -S $tRNA --un $unaligned 2>> $log
done
#gzip *.fq

cd ../../05-tRNA_Filter/RNASeq/
mkdir fastqc

for file in *.fq; do
fastqc $file -o ./fastqc
done



#Check final preprocessed output
#for file in *.log; do
#	TotalReads=`grep 'reads processed:' $file | awk '{print $0}'`
#	tRNA=`grep 'at least' $file | awk '{print $0}'`
#	RNA=`grep 'failed' $file | awk '{print $0}'`
#echo -e "${file%_.log}:\t"${TotalReads} >> tRNA_filtering_output.txt
#echo -e "${file%_.log}:\t"${tRNA} >> tRNA_filtering_output.txt
#echo -e "${file%_.log}:\t"${RNA} >> tRNA_filtering_output.txt
#done

###6)Mapping
#generate Insert files; needed for bed2 file generation later on
mkdir ../../06-Mapping/RNASeq

gzip *fq_RNASeq.fq

for file in *fq_RNASeq.fq.gz; do
# Loose Settings
echo "Start Loose of $file"
STAR --runThreadN 8 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --seedSearchStartLmax 30 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0  --outFilterMatchNmin 20 --outSAMattributes All --genomeDir /Work2/Planarian_genome/STAR_index/ --readFilesIn $file --readFilesCommand zcat --outFileNamePrefix "../../06-Mapping/RNASeq/loose_parameters/${file%fq_RNASeq.fq.gz}" --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
# Default Settings & Output of unmapped into file
echo "Start Default of $file"
STAR --runThreadN 8 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outSAMattributes All --genomeDir /Work2/Planarian_genome/STAR_index/ --readFilesIn $file --readFilesCommand zcat --outFileNamePrefix "../../06-Mapping/RNASeq/default_parameters/${file%fq_RNASeq.fq.gz}" --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
done

#a) Map vs TE Elements

##build STAR index for TE fasta from repBase
#cd /Work2/Planarian_genome/Dresden/planarian_genome/RepBase/
#mkdir STAR_INDEX_TE
#STAR --runMode genomeGenerate --genomeDir './STAR_INDEX_TE' --genomeFastaFiles consensus_planarian_TE_from_repBase.fa --runThreadN 8 --genomeSAindexNbases 9

#Perform Alignment of unmapped reads on TE Index
#cd to either /Work2/Tim_D_Sequencing/(...)/06-Mapping/RNASeq/loose_parameters|default_parameters/
for file in *Unmapped.out.mate1; do
STAR --runThreadN 4 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outSAMattributes All --genomeDir /Work2/Planarian_genome/Dresden/planarian_genome/RepBase/STAR_INDEX_TE/ --readFilesIn $file --outFileNamePrefix "./TE_STAR_ALIGN/${file%.out.mate1}" --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
done

cd /Work2/Tim_D_Sequencing/
# Loose Settings with --seedSearchStartLmax 50 
for file in *fq_RNASeq.fq.gz; do
echo "Start Loosen_but_tigther of $file"
STAR --runThreadN 8 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0  --outFilterMatchNmin 20 --outSAMattributes All --genomeDir /Work2/Planarian_genome/STAR_index/ --readFilesIn $file --readFilesCommand zcat --outFileNamePrefix "../../06-Mapping/RNASeq/loose_parameters_50split/${file%fq_RNASeq.fq.gz}" --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
done

# Default_Sam
for file in *fq_RNASeq.fq.gz; do

echo "Start Loosen_but_tigther of $file"
STAR --runThreadN 8 --outSAMattributes All --genomeDir /Work2/Planarian_genome/STAR_index/ --readFilesIn $file --readFilesCommand zcat --outFileNamePrefix "../../06-Mapping/RNASeq/Default_Sam/${file%fq_RNASeq.fq.gz}" --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
done

##-----------------------------------------------------------

## Continuing with loose_parameters_50split until "Default Sam" is ready. 

cd /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/06-Mapping/RNASeq/loose_parameters_50split/
# build index files of 
for file in *.bam; do
samtools index -b $file ${file%.bam}.bai
done
for file in *.bam; do
bamCoverage -p 8 -b $file -o ../../../07-IGV/RNASeq/Default_Sam/${file%.bam}.bigwig
done

for file in *.bam; do
bamCoverage -p 8 -b $file --normalizeUsing RPKM --filterRNAstrand reverse -o /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/07-IGV/RNASeq/RPKM_normalized/${file%.bam}.reverse.bw
bamCoverage -p 8 -b $file --normalizeUsing RPKM --filterRNAstrand forward -o /Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/07-IGV/RNASeq/RPKM_normalized/${file%.bam}.forward.bw
done