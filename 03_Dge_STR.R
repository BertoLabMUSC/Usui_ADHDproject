# Load libraries
suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(tidyverse)
library(preprocessCore)
library(future.apply)
library(DESeq2)
library(pheatmap)
})

dir.create("str_dge")

load("futcounts/Expression_Input_STR.RData")

pd <- data.frame(row.names=colnames(exp_str), Genotype = as.factor(do.call(rbind,strsplit(colnames(exp_str),"_"))[,1]))

# Filter the expression by condition
filter=apply(rpkm_str, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5)))
count_filt <- exp_str[filter,]
rpkm_filt <- rpkm_str[filter,]

logCPM <- log2(rpkm_filt+1)
p <- normalize.quantiles(as.matrix(logCPM))
rownames(p) <- rownames(logCPM)
colnames(p) <- colnames(logCPM)

pdf("str_dge/PCA_STR.pdf",width=6,height=6,useDingbats=FALSE)
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x,Genotype=pd$Genotype,ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Genotype",palette=c("black","red"), 
          shape = "Genotype", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic() 
dev.off()


# DGE 
dds <- DESeqDataSetFromMatrix(countData = count_filt,colData = pd,design = ~ Genotype)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=FALSE)

STR_FullTab <- as.data.frame(results(dds,contrast=c("Genotype","ADHD","Control"),cooksCutoff=FALSE,independentFiltering = FALSE)) %>%
          rownames_to_column("Gene")


STR_DGE <- STR_FullTab %>%
                  mutate(Abs = abs(log2FoldChange)) %>%
                  filter(padj < 0.05 & Abs > 0.3) %>%
                  arrange(desc(Abs))


save(STR_FullTab,STR_DGE,pd, file = "str_dge/STR_Dge_Data.RData")


openxlsx::write.xlsx(STR_FullTab, 
                     file = "str_dge/STR_DGE_FullStats.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(STR_DGE, 
                     file = "str_dge/STR_DGE.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")


# Input for viz Plot
df <- STR_FullTab %>% 
        mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% 
        mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(padj) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
mat <- p[rownames(p)%in% top_labelled$Gene,] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("str_dge/Boxplots_TopGenes_STR.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
 palette = c("red", "black")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=3) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("str_dge/Vulcano_Plot_STR.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "log2FoldChange", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"))+
      theme(legend.position="none")+
      ylim(0,30) + xlim(-5,+5)
dev.off()

# heatmap
mat <- p[rownames(p)%in% STR_DGE$Gene,]
anno <- pd
Genotype        <- c("red", "black")
names(Genotype) <- c("Control", "ADHD")
anno_colors <- list(Genotype = Genotype)
pdf("str_dge/Heatmap_STR.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()










