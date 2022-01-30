library(DESeq2)
library(Rsubread)
# Make a GTF from the GFF3 
library(rtracklayer)
library(EnhancedVolcano)
#?import.gff3#
#test_path <- "/Work2/Planarian_genome_2021/Planarian_genome/"
#test_gff3 <- file.path(test_path, "smes_v2_hconf_FullAnnotation.gff3")
#test <- import(test_gff3)
#export(test, "smes_v2_hconf_FullAnnotation.gtf", "gtf")
?plotMA
# use featureCounts from Rsubread to count exons... 
setwd('/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08.2-DE/piRNASeq')

bamfiles <- list.files("./BAM_files/", pattern='*.bam', full.names = TRUE)
?featureCounts
#for (mapped in c("0-antisense_mapped", "0-sense_mapped")) {
#for (mapped in c("0-sense_mapped")) {
#  if (mapped=="0-antisense_mapped"){
#    mappid = "antisense"}
#  else {
#  mappid = "sense"
#  }
mapped="0-sense_mapped"
mappid = "sense"
  ###/pdens0.3_swsize3000_clsize3000.merged.clusterXhconfXRNASeq_## -> This was the final one used!
  hconfannotation_file="/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/12-FINAL_Relevant_Cluster_Files/pdens0.3_swsize3000_clsize3000.merged.clusterXhconfXRNASeq_sense_100.gtf"
  #annotation_file=paste("/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08.2-DE/piRNASeq/1.2-cluster_gtfs_final/",mapped,"/pdens0.3_swsize3000_clsize3000.merged.clusterXhconfXRNASeq_",mappid,"_100.gtf", sep = "")
  count <- featureCounts(bamfiles, 
                        isGTFAnnotationFile=TRUE,
                        GTF.featureType = "piRNA_cluster",
                        GTF.attrType = "CLUSTER",
                        annot.ext = hconfannotation_file,
                        allowMultiOverlap=TRUE, 
                        largestOverlap=TRUE,
                        nthreads=8,
                        strandSpecific = 1,
                        countMultiMappingReads=TRUE,
                        isPairedEnd=FALSE,
                        )
  head(count$annotation)
  head(count$annotation$Length)
  head(count$counts)
  RNASeqOnlyannotation_file="/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/12-FINAL_Relevant_Cluster_Files/pdens0.3_swsize3000_clsize3000.merged.clusterXRNASeq_sense_100.gtf"
  countxRNASeq <- featureCounts(bamfiles, 
                         isGTFAnnotationFile=TRUE,
                         GTF.featureType = "piRNA_cluster",
                         GTF.attrType = "CLUSTER",
                         annot.ext = RNASeqOnlyannotation_file,
                         allowMultiOverlap=TRUE, 
                         largestOverlap=TRUE,
                         nthreads=8,
                         strandSpecific = 1,
                         countMultiMappingReads=TRUE,
                         isPairedEnd=FALSE,
  )
  head(countxRNASeq$annotation)
  head(countxRNASeq$annotation$Length)
  head(countxRNASeq$counts)

  dir.create("0-featurecounts")
  dir.create(paste("0-featurecounts/",mappid, sep = ""))

  write.table(x = data.frame(count$annotation$GeneID, count$counts, stringsAsFactors = FALSE), 
            file = paste('0-featurecounts/',mappid,'/unique-counts.txt', sep = ""),
            quote = FALSE,sep = "\t", row.names = FALSE)
  
  write.table(x = data.frame(countxRNASeq$annotation$GeneID, countxRNASeq$counts, stringsAsFactors = FALSE), 
              file = paste('0-featurecounts/',mappid,'/unique-counts_RNASeqintersect_only.txt', sep = ""),
              quote = FALSE,sep = "\t", row.names = FALSE)

  counts <- read.csv(paste('0-featurecounts/',mappid,'/unique-counts.txt', sep = ""), header=T, sep = "\t")

  rownames(counts) <- counts$count.annotation.GeneID
  counts <- counts[,-1]
  colnames(counts) <- c("ecoli_r1", "ecoli_r2", "ecoli_r3","polyic_r1","polyic_r2", "polyic_r3", "wt_r1", "wt_r2", "wt_r3")

  counts <- as.matrix(counts)

  ecoli_set = counts[,c(1:3, 7:9)]
  polyic_set = counts[,c(4:6,7:9)]


# Differential gene expression analysis using DESeq2. 
# For more information read manual https://bioconductor.org/packages/release/bioc/html/DESeq2.html

## Ecoli piRNA
# Generate a dataframe with a set of conditions that you are going to use in your analysis. The sequence of the conditions should be the same as the columns in your dataframe
  colData <- c("ecoli", "ecoli", "ecoli", "wt", "wt","wt")
  colData <- data.frame("condition" = colData)
dds <- DESeqDataSetFromMatrix(countData = ecoli_set,
                              colData = colData,
                              design = ~ condition)
colData(dds)$condition=relevel(colData(dds)$condition,
                               ref="wt")
?DESeq
?results
dds
dds <- DESeq(dds)
res <- results(dds)
ecolires <- res
ecolires
summary(res)
dir.create("3-Analysis")
dir.create(paste("./3-Analysis/",mappid, sep = ""))

# Analyse the results with MA plot and dispersion plot. Write the plot into pdf file.
dir.create("2-DE")
dir.create(paste("2-DE/",mappid, sep = ""))
pdf(file=paste("./2-DE/",mappid,"/ecoli_pi.pdf", sep = ""), width = 10, height = 5)
plotMA(dds,ylim=c(-2,2))
dev.off()

pdf(file=paste("./2-DE/",mappid,"/ecoli_pi_dispersion.pdf", sep = ""), width = 8, height = 5)
plotDispEsts(dds)
dev.off()

# Filter significant genes with padj less than 0,1 and write data into the txt file

sum(res$padj < 0.1, na.rm=TRUE)
de_significant=res[which(res$padj < 0.1), ]
head(de_significant[order(de_significant$padj), ])
de <- as.data.frame(de_significant)

write.table(de, file=paste("./3-Analysis/",mappid,"/de_unique_ecoli.txt", sep = ""), row.names = T, col.names = T, sep="\t")

## PolyIC
colData <- c("polyic", "polyic", "polyic", "wt", "wt","wt")
colData <- data.frame("condition" = colData)
dds <- DESeqDataSetFromMatrix(countData = polyic_set,
                              colData = colData,
                              design = ~ condition)
colData(dds)$condition=relevel(colData(dds)$condition,
                               ref="wt")
dds
dds <- DESeq(dds)
res <- results(dds)
polyicres <- res
polyicres
summary(res)

# Analyse the results with MA plot and dispersion plot. Write the plot into pdf file.
pdf(file=paste("./2-DE/",mappid,"/polyic_pi.pdf", sep = ""), width = 10, height = 5)
plotMA(dds,ylim=c(-2,2),)
dev.off()

pdf(file=paste("./2-DE/",mappid,"/polyic_pi_dispersion.pdf", sep = ""), width = 8, height = 5)
plotDispEsts(dds)
dev.off()

sum(res$padj < 0.1, na.rm=TRUE)
de_significant=res[which(res$padj < 0.1), ]

head(de_significant[order(de_significant$padj), ])
de <- as.data.frame(de_significant)

write.table(de, file=paste("./3-Analysis/",mappid,"/de_unique_polyIC.txt", sep = ""), row.names = T, col.names = T, sep="\t")

#}
selectiontoplot = c("")
### Volcano plots 
pdf(file="./2-DE/ecoli_volcano.pdf", width = 6, height = 8)
EnhancedVolcano(ecolires,
                lab = rownames(ecolires),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.15,
                title = 'E. coli treated vs Untreateed',
                selectLab = selectiontoplot,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                xlim=c(-2,2),
                ylim=c(0,30),
                col=c('#797979', '#797979', '#4878d0', '#d65f5f'),
                colConnectors = '#701f57')
dev.off()

pdf(file="./2-DE/polyic_volcano.pdf", width = 6, height =8)
EnhancedVolcano(polyicres,
                lab = rownames(polyicres),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.15,
                title = 'Poly(I:C) treated vs Untreateed',
                selectLab = selectiontoplot,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                xlim=c(-2,2),
                ylim=c(0,30),
                col=c('#797979', '#797979', '#4878d0', '#d65f5f'),
                colConnectors = '#701f57')
dev.off()

