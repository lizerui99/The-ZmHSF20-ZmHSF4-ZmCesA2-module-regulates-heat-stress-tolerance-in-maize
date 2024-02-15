rm(list=ls())

library(WGCNA)
library(tidyverse)
# setwd('C:/Users/Public/Desktop/Training/wgcna')
options(stringsAsFactors = FALSE) # Multi-threading is allowed
enableWGCNAThreads()

RNAseq_voom <- read.delim('d45wgcna.tab',header=T,row.names = 1,sep="\t",na.strings = "-")



datExprOri=t(RNAseq_voom)    # Transposition, one sample per row, one gene per column

#======================= Filtering of samples Filtering of samples with more than 10% number of missing genes===========
nGenes = ncol(datExprOri)      # Number of genes
NumberMissingBySample=apply(is.na(data.frame(datExprOri)),1, sum)   # of NA's calculated per line
KeepSample= NumberMissingBySample<0.1*nGenes  # Determine which samples have a missing rate<10%，
table(KeepSample)  #Statistics on the number of samples filtered

datExpr=datExprOri[KeepSample,]  # Expression Matrix Filtering Samples
SampleName=rownames(datExprOri)[KeepSample] # Sample Name Variable Filtration Sample

#=====================Filtering for genes with zero variance, and genes with more than 10% missing samples=====
nSamples = nrow(datExpr)  # Number of statistical samples
KeepGenes= variancedatExpr>0 & no.missingdatExpr<0.1*nSamples  # Retain genes with variance not equal to 0 and with less than 10% deletions
table(KeepGenes) # Filtering statistics

datExpr=datExpr[, KeepGenes] 
GeneName=colnames(datExpr) 


# Draw sample trees based on hierarchical clustering
tree=hclust(dist(datExpr),method = 'average')
pdf(file='sample clustering.pdf',w=10,h=7)
plot(tree,xlab="", sub="", cex = 0.7,main="Sample Clustering") 
dev.off()


#================Setting the network build parameter selection range======
powers = c(1:30)  
cor <- WGCNA::cor 
sft = pickSoftThreshold(datExpr, RsquaredCut = 0.9,powerVector = powers, verbose = 5, networkType = 'unsigned')  
sft$powerEstimate


sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")
abline(h=0.85,col="blue")

# Plot the average connectivity of each candidate beta
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#===============Automatic construction of WGCNA models==================
#========The expression matrix is converted to an adjacency matrix, and then the adjacency matrix is converted to a topology matrix, identifying the module=============
net = blockwiseModules(
  datExpr,
  power = 9,#sft$powerEstimate,         
  maxBlockSize = 15000,              
  TOMType = "unsigned",                 
  deepSplit = 2, minModuleSize = 30,    
  mergeCutHeight = 0.25,                
  numericLabels = TRUE,                       
  pamRespectsDendro = FALSE,  
  saveTOMs = TRUE,                     
  saveTOMFileBase = "FPKM-TOM",
  loadTOMs = FALSE,
  verbose = 3                            
)

# Get module color
moduleColors = labels2colors(net$colors)
# Get the feature values for each module
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Output the module in which each gene is located and the KME value with that module
#file.remove('All_Gene_KME.txt')
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=as.data.frame(MEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F)
}

#===============模块绘图==================
# log10（FPKM+1）
expColor=t(numbers2colors(scale((datExpr)),colors=blueWhiteRed(100),naColor="grey"))
colnames(expColor)=rownames(datExpr)
# Draw a tree diagram of the genes, a module diagram, and an expression graph for each sample
pdf("wgcna.dendroColors.pdf",height = 7,width = 9)

plotDendroAndColors(net$dendrograms[[1]], 
                    colors=cbind(moduleColors[net$blockGenes[[1]]],expColor),
                    c("Module",colnames(expColor)),
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    cex.rowText=0.5)
dev.off()

# Plotting the adjacency matrix between two-by-two modules
pdf("wgcna.adjacency.heatmap.pdf",height = 10,width = 9)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",plotDendrograms = F,
                      marDendro = c(4,4,2,4))
dev.off()


# Plot heatmaps of expression values and bar charts of eigenvalues for all modules
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=MEs[,paste("ME",module,sep="")]
  pdf(paste("wgcna.", module, ".express.barplot.pdf", sep=""),height = 7,width = 9)
  par(mfrow=c(2,1),mar=c(0.3,5.5,3,2))
  plotMat(t(scale(datExpr[,moduleColors==module])),
          rlabels=F,main=module,cex.main=2,clabels=F)
  
  par(mar=c(5,4.2,0,0.7))
  barplot(ME,col=module,main="",cex.main=2,ylab="eigengene expression",xlab="sample")
  dev.off()
}
# Save the results in cytoscape's input file format

# Read in the TOM matrix saved at build time
load('FPKM-TOM-block.1.RData') 
TOM=as.matrix(TOM)

# or recalculate the TOM matrix
#TOM = TOMsimilarityFromExpr(datExpr, power =sft$powerEstimate,TOMType = "unsigned"); 

for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", module, ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", module, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.1,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
}

probes = colnames(datExpr)
dimnames(TOM) <- list(probes, probes)
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("edges-", ".edges.txt", sep=""),
                               nodeFile = paste("nodes-", ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0,
                               nodeNames = probes, nodeAttr = moduleColors)

adjacency = adjacency(datExpr,power=9)
kIM <- intramodularConnectivity(adjacency, moduleColors)
kIM$module <- moduleColors
write.csv(kIM,'kim.csv')

