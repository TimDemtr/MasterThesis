library(topGO)
library(plotrix)
library(Biostrings)
library(Rgraphviz)
options(max.print = 999999999)
fisher.mode="parentChild"
#fisher.mode="classic"
########################GO-Analysis
setwd("/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08.2-DE/piRNASeq/4-GOAnalysis")
dir.create("./GO_out_graphs")
dir.create("./GO_out_graphs/GO_Networks")
dir.create("./GO_out_graphs/barplot_out/percentageordered", recursive = TRUE)
dir.create("./GO_out_graphs/barplot_out/fisherordered", recursive = TRUE)
dir.create("./GO_out_graphs/TSVs")
dir.create("./GO_out_graphs/TSVs/allRes")
dir.create("./GO_out_graphs/TSVs/GO_Stats")
dir.create("./GO_out_graphs/TSVs/sigGenes")
dir.create("./GO_out_graphs/TSVs/GOtoGene")
dir.create("./GO_out_graphs/histogramm")
dir.create("./GO_out_graphs/TSVs/Plotted_Terms")
#Read in GO Annotation
#geneID2GO <- read.csv('/Work2/Planarian_genome_2021/SMEST2G_2GO/smes_v2_hconf_toGO.tsv', header=F, sep = '\t', stringsAsFactors = T)
geneID2GO <- read.csv('/Work2/Planarian_genome_2021/SMEST2G_2GO/smes_v2_hconf_only1SMEST_toGO.tsv', header=F, sep = '\t', stringsAsFactors = T) # only SMESTXXXXX1.1...
#geneID2GO <- read.csv('smes_v2_hconf_toGO.tsv', header=F, sep = '\t', stringsAsFactors = T)
#geneID2GO <- geneID2GO[,c(3,5)]
colnames(geneID2GO)=c("Transcript.ID","GO-ID","Description","GO-Term")

#create Gene and GO Universe 
geneID2GO$`GO-Term` <- as.character(geneID2GO$`GO-Term`)                                                    # has to be done as character
geneID2GO$Transcript.ID <- as.character(geneID2GO$Transcript.ID)
geneID2GO$`GO-ID` <- as.character(geneID2GO$`GO-ID`)
geneUniverse <- geneID2GO$Transcript.ID
geneUniverse <- unique(geneUniverse)
gene.to.GO <- split(geneID2GO$`GO-ID`, geneID2GO$Transcript.ID)                                             # creating a list of all associated GO-Terms for each gene
str(head(gene.to.GO))

#GO-ID to GO-Term Mapping
GOs <- as.data.frame(Term(GOTERM))
GOs[2] <- rownames(GOs)
colnames(GOs)=c("GO-Term","GO-ID")
GOs$`GO-Term` <- as.character(GOs$`GO-Term`)
rownames(GOs)=NULL
write.table(GOs, file= "GOdescriptions.txt", col.names = F, row.names = F, quote = F, sep='\t')
#output_path="/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08.2-DE/piRNASeq/4-GOAnalysis/GO_out_graphs/"
#file_paths = "/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08.2-DE/piRNASeq/4-GOAnalysis/DE_unique_Transcripts/"
#all_files <- setdiff(list.files(file_paths), list.dirs(file_paths,recursive = FALSE, full.names = FALSE))
#all_files

output_path="./GO_out_graphs/"
tsv_path="./GO_out_graphs/TSVs/"
file_paths = "./DE_unique_Transcripts/"
all_files <- setdiff(list.files(file_paths), list.dirs(file_paths,recursive = FALSE, full.names = FALSE))
all_files

#Input-Data: Gene-Symbols from all reps

for (inputfile in all_files) { 

  all_genes <- read.table(paste(file_paths,inputfile,sep=""), header=F, sep = "\t")
  #AC_Transcripts_cutExon
  
  genesOfInterest <- as.character(all_genes$V1)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  str(geneList)
  
  GOdata_all <- new("topGOdata", description = paste("GO-Analysis of Genes in",inputfile), ontology = "BP",
                    allGenes = geneList, 
                    annot = annFUN.gene2GO, gene2GO = gene.to.GO)
  
  GOdata_all
  sig <- as.data.frame(sigGenes(GOdata_all)) # Liste der signifikanten Gene
  write.table(sig, file= paste(tsv_path,"sigGenes/",inputfile,".sigGenes.txt", sep=""), col.names = F, row.names = F, quote = F, sep='\t')
  
  
  #resultFisher <- runTest(GOdata_all, algorithm = "classic", statistic = "fisher")
  resultFisher <- runTest(GOdata_all, algorithm = fisher.mode, statistic = "fisher")
  #elim.ks <- runTest(GOdata_all, algorithm = "elim", statistic = "ks")
  #weight01.t <- runTest(GOdata_all, algorithm = "weight01", statistic = "t")
  
  allRes <- GenTable(GOdata_all, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(resultFisher@score)) #KS=elim.ks, weight=weight01.t, 
  ?GenTable
  #write.table(allRes, file='./GO-Terms_p_less0.01.csv', col.names = T, row.names = F, quote = F, sep='\t')
  
  #pdf-Output of GO-Analysis
  pdf(file = paste(output_path,"GO_Networks/",inputfile,".",fisher.mode,".pdf",sep=""), width = 10, height = 10, )
  nodes=showSigOfNodes(GOdata_all, score(resultFisher), useInfo ='all') #rectangles show most significant ones. significance from dark red (most) to bright yellow (least)
  dev.off()
  pvalGO <- score(resultFisher)
  pdf(file = paste(output_path,"histogramm/",inputfile,".",fisher.mode,".hist.pdf", sep=""), width = 10, height = 10, )
  hist(pvalGO[which(pvalGO<1)],100, xlab = "p-values", main = "p-value Frequency of GO-Terms with p-value < 1")
  dev.off()
  pdf(file = paste(output_path,"histogramm/",inputfile,".",fisher.mode,".histincl1.pdf", sep=""), width = 10, height = 10, )
  hist(pvalGO,100, xlab = "p-values", main = "p-value Frequency of GO-Terms with p-value < 1")
  dev.off()
  geneData(resultFisher)
  
  
  ###################################################################
  #Statistics for all annotated GO-Terms
  go <- usedGO(GOdata_all) # Liste der verfügbaren GO-Terms
  GO_Stats <- termStat(GOdata_all, go)#sel.terms) # gibt die Statistik für die eingegebenen GO-Terms aus
  
  #Mapping GO-Terms into GO-Stats
  termcol <- c()
  for(i in 1:length(rownames(GO_Stats))){
    t=GOs$`GO-Term`[which(GOs$`GO-ID`==rownames(GO_Stats)[i])]
    termcol <- append(termcol,t)
  }
  GO_Stats[4] <- termcol
  colnames(GO_Stats)[4]=c("GO-Term")
  
  
  #Calculating Percentage of present genes
  pct <- GO_Stats$Significant/GO_Stats$Annotated
  GO_Stats[5] <- pct
  colnames(GO_Stats)[5]=c("Percentage")
  
  #Mapping fisher p-values to GO-IDs
  length(rownames(GO_Stats))
  completeRes <- GenTable(GOdata_all, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(rownames(GO_Stats)))
  fishercol <- c()
  for(i in 1:length(rownames(GO_Stats))){
    f=as.numeric(completeRes$classicFisher[which(completeRes$GO.ID==rownames(GO_Stats)[i])])
    fishercol <- append(fishercol,f)
  }
  GO_Stats[6] <- fishercol
  colnames(GO_Stats)[6]=c("classicFisher")
  head(GO_Stats)
  head(allRes)
  GO_Stats <- GO_Stats[order(GO_Stats$classicFisher, decreasing = F),]
  ################GO-Plots
  #Top Genes sorted by fisher.mode fisher
  GOStatPlot <- GO_Stats[which(GO_Stats$Annotated>=10),]
  GOStatPlot <- GOStatPlot[which(GOStatPlot$classicFisher < 0.05),]
  GOStatPlot <- GOStatPlot[order(GOStatPlot$classicFisher, decreasing = F),]
  write.table(GO_Stats, file= paste(tsv_path,"GO_Stats/",inputfile,".",fisher.mode,".top35.GO_Stats.tsv", sep=""), col.names = T, row.names = T, quote = F, sep='\t')
  GOStatPlot <- GOStatPlot[which(GOStatPlot$Expected>=0.1),]
  Plotting <- GOStatPlot[c(1:35),]
  #Plotting <- GOStatPlot[c(1:15),]
  Plotting <- Plotting[order(Plotting$classicFisher, decreasing=T),]
  head(GOStatPlot)
  head(Plotting)
  write.table(Plotting, file= paste(tsv_path,"Plotted_Terms/",inputfile,".",fisher.mode,".plotted.tsv", sep=""), col.names = T, row.names = T, quote = F, sep='\t')
  pdf(file = paste(output_path,"barplot_out/fisherordered/",inputfile,".",fisher.mode,".top35.pdf", sep=""), width = 10, height = 10, )
  par(mar=c(5,32,21,5))
  barplot(Plotting$Percentage,
          names.arg = paste(row.names(Plotting),Plotting$`GO-Term`, sep=" : "), horiz = TRUE, #über col könnte noch eine Farbcodierung eingefügt werden
          xlab="Percentage associated Genes", cex.names = 0.7, 
          xlim=c(0,1), las=1, ylim=c(0,70), width=3) 
  mtext(paste("Top 50 GO-Terms by % of deregulated Genes\n in",inputfile), at=getFigCtr()[1], line=18, font=2)
  dev.off()
  
  
  #Top Genes sorted by percentage
  GOStatPlot <- GO_Stats[which(GO_Stats$Annotated>=10),]
  GOStatPlot <- GOStatPlot[which(GOStatPlot$classicFisher < 0.05),]
  GOStatPlot <- GOStatPlot[order(GOStatPlot$classicFisher, decreasing = F),]
  #GOStatPlot <- GOStatPlot[order(GOStatPlot$Percentage, decreasing = T),]
  write.table(GO_Stats, file= paste(tsv_path,"GO_Stats/",inputfile,".",fisher.mode,".top35.GO_Stats.tsv", sep=""), col.names = NA, row.names = T, quote = F, sep='\t')
  GOStatPlot <- GOStatPlot[which(GOStatPlot$Expected>=0.1),]
  Plotting <- GOStatPlot[c(1:35),]
  #Plotting <- GOStatPlot[c(1:15),]
  Plotting <- Plotting[order(Plotting$Percentage, decreasing=F),]
  
  head(GOStatPlot)
  head(Plotting)
  write.table(Plotting, file= paste(tsv_path,"Plotted_Terms/",inputfile,".",fisher.mode,".plotted_percentage.tsv", sep=""), col.names = NA, row.names = T, quote = F, sep='\t')
  pdf(file = paste(output_path,"barplot_out/percentageordered/",inputfile,".",fisher.mode,".top35_percentage.pdf",sep=""), width = 10, height = 10, )
  par(mar=c(5,32,21,5))
  barplot(Plotting$Percentage,
          names.arg = paste(row.names(Plotting),Plotting$`GO-Term`, sep=" : "), horiz = TRUE, #über col könnte noch eine Farbcodierung eingefügt werden
          xlab="Percentage associated Genes", cex.names = 0.7, 
          xlim=c(0,1), las=1, ylim=c(0,70), width=3) 
  mtext(paste("Top 50 GO-Terms by % of deregulated Genes\n in",inputfile), at=getFigCtr()[1], line=18, font=2)
  dev.off()
  
  
  
  allRes <- allRes[order(allRes$classicFisher, decreasing=F),]
  write.table(allRes, file= paste(tsv_path,"allRes/",inputfile,".",fisher.mode,".top35.allRes.tsv", sep=""), col.names = T, row.names = T, quote = F, sep='\t')
  ?genesInTerm
  genesInTermOutput <- genesInTerm(GOdata_all)
  capture.output((genesInTermOutput),  file= paste(tsv_path,"GOtoGene/",inputfile,".GOtoGene.txt", sep=""))

  
}
#go <- usedGO(GOdata_all)
#usedGO(GOdata_all)
#length(go)
#sample(go, 10)
#samplesel.terms <- sample(go, 10)
#ann.genes <- genesInTerm(GOdata_all, sel.terms)
#str(ann.genes)
#num.ann.genes <- countGenesInTerm(GOdata_all)
#?countGenesInTerm
#str(num.ann.genes)
#?par
#?barplot
#par(mar=c(5,23.5,17,5))
#head(allRes)
#sigGenes(GOdata_all,)
#?showGroupDensity
#goID <- "GO:0033209"
#print(showGroupDensity(GOdata_all, goID, ranks = TRUE,rm.one = F))
#dev.off()