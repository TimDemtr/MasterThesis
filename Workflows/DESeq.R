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

# use featureCounts from Rsubread to count exons... 


selectiontoplot =  c('SMESG000000502.1', 'SMESG000024584.1', 'SMESG000027263.1', 'SMESG000069927.1', 'SMESG000069927.1',
                     'SMESG000067336.1','SMESG000027263.1', )

setwd('/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08.2-DE/RNASeq')

bamfiles <- list.files("./BAM_files/", pattern='*out.bam', full.names = TRUE)
#bamfiles <- bamfiles[1]
annotation.file = "/Work2/Planarian_genome_2021/Planarian_genome/smes_v2_hconf_FullAnnotation.gtf"

count <- featureCounts(bamfiles, 
                       isGTFAnnotationFile=TRUE,
                       GTF.featureType = "gene",
                       GTF.attrType = "ID",
                       annot.ext = annotation.file,
                       allowMultiOverlap=TRUE, 
                       largestOverlap=TRUE,
                       nthreads=8,
                       strandSpecific = 2,
                       countMultiMappingReads=TRUE,
                       isPairedEnd=FALSE,
)

head(count$annotation)
head(count$annotation$Length)
head(count$counts)



dir.create("0-featurecounts")

write.table(x = data.frame(count$annotation$GeneID, count$counts, stringsAsFactors = FALSE), 
            file = '0-featurecounts/unique-counts.txt',
            quote = FALSE,sep = "\t", row.names = FALSE)

counts <- read.csv('0-featurecounts/unique-counts.txt', header=T, sep = "\t")

rownames(counts) <- counts$count.annotation.GeneID
counts <- counts[,-1]
colnames(counts) <- c("ecoli_r1", "ecoli_r2", "ecoli_r3", "wt_r1", "wt_r2", "wt_r3","polyic_r1","polyic_r2", "polyic_r3")

counts <- as.matrix(counts)

ecoli_set = counts[,c(1:3, 4:6)]
polyic_set = counts[,c(7:9, 4:6)]

# Differential gene expression analysis using DESeq2. 
# For more information read manual https://bioconductor.org/packages/release/bioc/html/DESeq2.html

# Differential gene expression analysis using DESeq2. 
# For more information read manual https://bioconductor.org/packages/release/bioc/html/DESeq2.html

## Ecoli RNASeq
# Generate a dataframe with a set of conditions that you are going to use in your analysis. The sequence of the conditions should be the same as the columns in your dataframe
colData <- c("ecoli", "ecoli", "ecoli", "wt", "wt","wt")
colData <- data.frame("condition" = colData)
dds <- DESeqDataSetFromMatrix(countData = ecoli_set,
                              colData = colData,
                              design = ~ condition)
#mcols(dds)$basepairs = count$annotation$Length
#fpkmvalues = fpkm(dds)
#write.table(x = fpkmvalues, 
#            file = '0-featurecounts/ecoli_fpkms.txt',
#            quote = FALSE,sep = "\t", row.names = TRUE)

colData(dds)$condition=relevel(colData(dds)$condition,
                               ref="wt")
dds
dds <- DESeq(dds)
res <- results(dds)
ecolires <- res
summary(res)
dir.create("3-Analysis")

# Analyse the results with MA plot and dispersion plot. Write the plot into pdf file.
dir.create("2-DE")
pdf(file="./2-DE/ecoli.pdf", width = 10, height = 5)
plotMA(dds,)#
dev.off()

pdf(file="./2-DE/ecoli_dispersion.pdf", width = 8, height = 5)
plotDispEsts(dds,ylim=c(1e-6,7e0))
dev.off()


# Filter significant genes with padj less than 0,1 and write data into the txt file
sum(res$padj < 0.1, na.rm=TRUE)
de_significant=res[which(res$padj < 0.1), ]
head(de_significant[order(de_significant$padj), ])
de <- as.data.frame(de_significant)
alle <- as.data.frame(res)

write.table(de, file="./3-Analysis/de_unique_ecoli.txt", row.names = T, col.names = T, sep="\t")
write.table(alle, file="./3-Analysis/de_alle_ecoli.txt", row.names = T, col.names = T, sep="\t")

## PolyIC
colData <- c("polyic", "polyic", "polyic", "wt", "wt","wt")
colData <- data.frame("condition" = colData)
dds <- DESeqDataSetFromMatrix(countData = polyic_set,
                              colData = colData,
                              design = ~ condition)
#mcols(dds)$basepairs = count$annotation$Length
#fpkmvalues = fpkm(dds)
#write.table(x = fpkmvalues, 
#            file = '0-featurecounts/polyic_fpkms.txt',
#            quote = FALSE,sep = "\t", row.names = TRUE)
colData(dds)$condition=relevel(colData(dds)$condition,
                              ref="wt")
dds
dds <- DESeq(dds)
res <- results(dds)
res
polyicres <- res
summary(res)

# Analyse the results with MA plot and dispersion plot. Write the plot into pdf file.
dir.create("2-DE")
pdf(file="./2-DE/polyic.pdf", width = 10, height = 5)
plotMA(dds,)#
dev.off()

pdf(file="./2-DE/polyic_dispersion.pdf", width = 8, height = 5)
plotDispEsts(dds,ylim=c(1e-6,7e0))
dev.off()



# Filter significant genes with padj less than 0,1 and write data into the txt file

sum(res$padj < 0.1, na.rm=TRUE)
de_significant=res[which(res$padj < 0.1), ]
head(de_significant[order(de_significant$padj), ])
de <- as.data.frame(de_significant)
alle <- as.data.frame(res)
write.table(de, file="./3-Analysis/de_unique_polyIC.txt", row.names = T, col.names = T, sep="\t")
write.table(alle, file="./3-Analysis/de_alle_polyIC.txt", row.names = T, col.names = T, sep="\t")

### Volcano plots 
selectiontoplot = c('SMESG000000502.1', 'SMESG000049844.1', 'SMESG000063594.1', 'SMESG000063621.1', 'SMESG000067250.1', 'SMESG000070162.1', 'SMESG000027308.1', 'SMESG000052389.1', 'SMESG000010558.1', 'SMESG000025566.1', 'SMESG000049783.1', 'SMESG000060403.1', 'SMESG000005043.1',
                    'SMESG000001806.1', 'SMESG000002673.1', 'SMESG000003379.1',)
                  # for extended list see polyIC
selectiontoplot = c('SMESG000049844.1', 'SMESG000063621.1', 'SMESG000027308.1', 'SMESG000010558.1', 'SMESG000049783.1', 'SMESG000060403.1')

pdf(file="./2-DE/ecoli_volcano.pdf", width =6, height = 8)
EnhancedVolcano(ecolires,
                lab = rownames(ecolires),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1.5,
                title = 'E. coli treated vs Untreateed',
                selectLab = selectiontoplot,
                labCol = 'black',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                xlim=c(-10,10),
                ylim=c(0,35),
                labSize =3,
                col=c('#797979', '#797979', '#4878d0', '#d65f5f'),
                colConnectors = '#701f57')
dev.off()


selectiontoplot = c('SMESG000000502.1', 'SMESG000049844.1', 'SMESG000063621.1', 'SMESG000067250.1', 'SMESG000070162.1', 'SMESG000037135.1', 'SMESG000065681.1', 'SMESG000010558.1', 'SMESG000024968.1', 'SMESG000012317.1', 'SMESG000012332.1', 'SMESG000016883.1', 'SMESG000035518.1', 'SMESG000057992.1', 'SMESG000014755.1', 'SMESG000014759.1', 'SMESG000052389.1',)
                    #'SMESG000001806.1', 'SMESG000002673.1', 'SMESG000003379.1', 'SMESG000003977.1', 'SMESG000007869.1', 'SMESG000015520.1', 'SMESG000016079.1', 'SMESG000017099.1', 'SMESG000018333.1', 'SMESG000020535.1', 'SMESG000026745.1', 'SMESG000031033.1', 'SMESG000037102.1', 'SMESG000037593.1', 'SMESG000042447.1', 'SMESG000042470.1', 'SMESG000052688.1', 'SMESG000060163.1', 'SMESG000061113.1', 'SMESG000065562.1', 'SMESG000066537.1', 'SMESG000066540.1', 'SMESG000066644.1', 'SMESG000074545.1', 'SMESG000076644.1', 'SMESG000006596.1', 'SMESG000023166.1', 'SMESG000033348.1', 'SMESG000033907.1', 'SMESG000052204.1')

selectiontoplot = c('SMESG000000502.1', 'SMESG000049844.1', 'SMESG000063621.1',  'SMESG000037135.1', 'SMESG000010558.1', 'SMESG000024968.1', 'SMESG000012332.1', 'SMESG000014755.1', 'SMESG000014759.1', 'SMESG000052389.1')
pdf(file="./2-DE/polyic_volcano.pdf", width = 6, height = 8)
EnhancedVolcano(polyicres,
                lab = rownames(polyicres),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1.5,
                title = 'Poly(I:C) treated vs Untreateed',
                selectLab = selectiontoplot,
                labCol = 'black',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                xlim=c(-10,10),
                ylim=c(0,35),
                labSize =3,
                col=c('#797979', '#797979', '#4878d0', '#d65f5f'),
                colConnectors = '#701f57',)
dev.off()


