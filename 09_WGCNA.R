# WGCNA
suppressPackageStartupMessages({
library(tidyverse)
library(WGCNA)
library(cluster)
library(ggplot2)
library(ggpubr)
})
enableWGCNAThreads()

dir.create("wgcna_output")

# Load tables 
load("Expression_Input.RData")

pd <- data.frame(row.names=colnames(rpkm), Genotype = as.factor(do.call(rbind,strsplit(colnames(rpkm),"_"))[,1]), 
	             Region = as.factor(do.call(rbind,strsplit(colnames(rpkm),"_"))[,2])) 

pd <- pd %>%
	  rownames_to_column("ID") %>%
	  unite("Sample",Genotype:Region,remove=F) %>%
	  mutate(Sample = as.factor(Sample),Genotype = as.factor(Genotype),Region = as.factor(Region)) %>%
	  column_to_rownames("ID")

datTraits <- model.matrix(~ .-1, data = pd) %>%
				as.data.frame()

# Filter the expression by condition
filter=apply(rpkm, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5) | all(x[4:6] >= 0.5) | all(x[7:9] >= 0.5) | all(x[10:12] >= 0.5)))
rpkm_filt <- rpkm[filter,]

datExpr <- t(log2(rpkm_filt + 1))

## Powers analysis
powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr,powerVector=powers,verbose = 5, blockSize= 16000, networkType = "signed") 
pdf("wgcna_output/SoftThresholdingPower_signed.pdf")
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

######################################################################################################################
############################################ MODULE construction #####################################################
############################################    signed network   #####################################################
PWR=22

net <- blockwiseModules(datExpr,
corType="bicor",
maxBlockSize = 16000,
networkType="signed",
minCoreKME = 0.4, 
minKMEtoStay = 0.5,
power=PWR, 
checkMissingData = TRUE,
minModuleSize=50,
nThreads=15,
TOMType = "signed",
TOMDenom = "mean",
deepSplit=4,
verbose=1,
mergeCutHeight=0.10,
reassignThreshold = 1e-10,
numericLabels=TRUE,
saveIndividualTOMs=TRUE,
saveConsensusTOMs=TRUE)
save(net,file="wgcna_output/Network_Data.RData")

moduleLabelsAutomatic=net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
MEsAutomatic=net$MEs
unique(moduleColorsAutomatic)
table(moduleColorsAutomatic)
write.table(moduleColorsAutomatic, "wgcna_output/Net_Colors.txt",sep="\t",quote=F)


KMEs<-signedKME(datExpr, net$MEs,corFnc = "cor", corOptions = "use = 'p'")
kme=data.frame(rownames(rpkm_filt), moduleColorsAutomatic, KMEs)
colnames(kme)[1]="Symbol"
rownames(kme)=NULL
write.table(kme,"wgcna_output/KMEs.txt",sep="\t",quote=F)

# Heatmap
moduleColorsIEGG=moduleColorsAutomatic
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsIEGG)$eigengenes
MEsIEGG = MEs0
MEsIEGG$MEgrey=NULL
modTraitCor= cor(MEsIEGG, datTraits,method="spearman")
write.table(modTraitCor,"wgcna_output/modTraitCor.txt",sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitP,"wgcna_output/modTraitP.txt",sep="\t",quote=F)
textMatrix = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor) 
par(mar = c(6, 8.5, 3, 3))
pdf("wgcna_output/Heatmap_DatTraits.pdf",width=7,height=7)
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), yLabels = names(MEsIEGG), ySymbols = names(MEsIEGG), colorLabels =FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix, setStdMargins = FALSE, cex.text = 0.4, zlim = c(-1,1), main = paste("Module Association"))
dev.off()

#EigenGeneNet
MET=orderMEs(MEs0)
pdf("wgcna_output/EigengeneNetworks.pdf")
plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.5,xLabelsAngle=90)
dev.off()

#GGplotInput
MEs0$Rows=colnames(rpkm_filt)
MEs0$Sample = pd$Sample
write.table(MEs0, "wgcna_output/Matrix_module_correlation.txt",sep="\t",quote=F)

# Adjacency matrix
Adj = adjacency(datExpr, power = PWR,type="signed",corFnc = "bicor")
moduleOutput <- data.frame(rownames(rpkm_filt))
moduleOutput[,2]<- as.factor(moduleColorsAutomatic)
intraCon <- intramodularConnectivity(Adj, moduleColorsAutomatic)
moduleOutput[,3]<-intraCon$kWithin
colnames(moduleOutput) <- c("Gene", "ModuleColor", "kWithin")
moduleOutput <- droplevels(moduleOutput[!(moduleOutput$ModuleColor=="grey"),])
moduleOutput$ModuleName <- plyr::mapvalues(moduleOutput$ModuleColor, as.character(sort(unique(moduleOutput$ModuleColor))), paste("Adhd","M",1:length(levels(moduleOutput$ModuleColor)),sep=""))
moduleOutput <- moduleOutput[order(moduleOutput$ModuleColor),]
write.table(moduleOutput, "wgcna_output/ModuleOutput.txt", sep="\t", quote=F,row.names=FALSE)

# TO connectivity (you need the table as single gene list column in the directory)
TOM = TOMsimilarityFromExpr(datExpr, power= PWR,corType = "bicor",networkType="signed",TOMType="signed",TOMDenom = "mean",nThreads = 15,verbose = 5, indent = 0)
colnames(TOM)=rownames(TOM)=colnames(datExpr)
save(TOM,file="wgcna_output/TOM_SIGNED.RData")
Connectivity=apply(TOM,1,sum)
save(Connectivity,file="wgcna_output/Connectivity.RData")

# NetDendro
x=cor(datTraits,datExpr,method="spearman")
x[x > 0.2]=6
x[x < -0.2]=30
x[ x > -0.2 & x < 0.2] = 27
listLabels2 <- t(labels2colors(x))
rownames(listLabels2) = colnames(x) 
colnames(listLabels2) = rownames(x)
plotColors=as.data.frame(cbind(moduleColorsAutomatic,listLabels2));
colnames(plotColors)=c("Module",names(datTraits)) 
geneTree = hclust(1-as.dist(TOM), method="average")
pdf("wgcna_output/NetworkDendrogram_ADHD.pdf",width=10,height=4)
plotDendroAndColors(geneTree, plotColors,colorHeight = 0.1, dendroLabels=FALSE,hang = 0.03, addGuide = TRUE, guideHang = 0.05 )
dev.off()

# CytoScape output
dir.create("wgcna_output/Cyto")
for(module in unique(moduleColorsAutomatic)){
inModule <- is.finite(match(moduleColorsAutomatic, module))
modTOM <- TOM[inModule, inModule]
cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste("wgcna_output/Cyto/CytoEdge",paste(module,collapse="-"),".txt",sep=""), nodeFile=paste("wgcna_output/Cyto/CytoNode",paste(module,collapse="-"),".txt",sep=""), weighted = TRUE, threshold = 0, nodeAttr = moduleColorsAutomatic[inModule], nodeNames = names(datExpr)[inModule])
}

# WGCNA barplot output
library(ggplot2)
library(reshape2)
library(RColorBrewer)
df=melt(MEs0)
dir.create("wgcna_output/Plot")
setwd("wgcna_output/Plot")
df$Rows=factor(df$Rows, levels = colnames(rpkm_filt))
df=df[!df$variable == "MEgrey", ]
# Function to make the barplot and saving as pdf according to the module name
cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
cols <- sample(cl, 12)
doPlot = function(sel_name) 
{
    df = subset(df, variable == sel_name)
    PLOT=ggplot(data=df, aes(x=Rows, y=value)) +
 geom_bar(aes(fill=Sample),stat="identity",position = "dodge")+
 scale_y_continuous(limits=c(-1,+1))+
 theme_bw()+
 theme(strip.text.x = element_text(size=12, face="bold"),
 strip.background = element_rect(colour="black", fill="#CCCCFF"))+
 scale_fill_manual(values = cols)+
 theme(axis.title.x = element_blank(),
            axis.text.x  = element_text(face="bold", size=6,angle = 45, hjust = 1))+
 theme(axis.title.y = element_blank(),
            axis.text.y  = element_text(face="bold", size=6))
    print(PLOT)
    ggsave(sprintf("%s.pdf", sel_name))
 }

lapply(unique(df$variable), doPlot)





