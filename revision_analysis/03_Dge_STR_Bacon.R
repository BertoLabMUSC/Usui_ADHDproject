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
library(sva)
library(qqman)
library(bacon)
library(UpSetR)
})

dir.create("str_dge_bacon")


# BACON analsis CTX
dge <- readxl::read_excel("original_dge/STR_DGE_FullStats.xlsx")

beta <- dge$log2FoldChange
se <- dge$lfcSE

bc <- bacon(effectsize = dge$log2FoldChange, standarderror = dge$lfcSE, teststatistics = dge$stat)


beta_corrected <- es(bc)
se_corrected <- se(bc)
pval_corrected <- pval(bc)

STR_FullTab_Bacon <- data.frame(
  Gene = dge$Gene,
  log2FoldChange_original = beta,
  lfcSE_original = se,
  pvalue_original = dge$pvalue,
  padj_original = dge$padj,
  log2FoldChange_corrected = beta_corrected,
  lfcSE_corrected = se_corrected,
  pvalue_corrected = pval_corrected,
  padj_corrected = p.adjust(pval_corrected,"BH")
)


STR_DGE_Original <- STR_FullTab_Bacon %>%
                  mutate(Abs = abs(log2FoldChange_original)) %>%
                  filter(padj_original < 0.05 & Abs > 0.3) %>%
                  arrange(desc(Abs))

STR_DGE_Bacon <- STR_FullTab_Bacon %>%
                  mutate(Abs = abs(log2FoldChange_corrected)) %>%
                  filter(padj_corrected < 0.05 & Abs > 0.3) %>%
                  arrange(desc(Abs))


openxlsx::write.xlsx(STR_FullTab_Bacon, 
                     file = "str_dge_bacon/STR_DGE_FullStats_BACON.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(STR_DGE_Original, 
                     file = "str_dge_bacon/STR_DGE_Original.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(STR_DGE_Bacon, 
                     file = "str_dge_bacon/STR_DGE_Bacon.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")



# Viz and stats
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



# Input for viz Plot
df <- STR_FullTab_Bacon %>% 
        mutate(LOG = -log10(padj_corrected), ABS = abs(log2FoldChange_corrected)) %>% 
        mutate(Threshold = if_else(padj_corrected < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(log2FoldChange_corrected > 0.3 & padj_corrected < 0.05 ~ "UpReg", 
                                     log2FoldChange_corrected < -0.3 & padj_corrected < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(padj_corrected) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
mat <- p[rownames(p)%in% top_labelled$Gene,] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("str_dge_bacon/Boxplots_TopGenes_STR_Bacon.pdf",width=6,height=5,useDingbats=FALSE)
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

pdf("str_dge_bacon/Vulcano_Plot_STR_Bacon.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "log2FoldChange_corrected", 
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
mat <- p[rownames(p)%in% STR_DGE_Bacon$Gene,]
anno <- pd
Genotype        <- c("red", "black")
names(Genotype) <- c("Control", "ADHD")
anno_colors <- list(Genotype = Genotype)
pdf("str_dge_bacon/Heatmap_CTX_Bacon.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()

# QQ plots
dir.create("str_dge_bacon/qqplots_comparative")

pdf("str_dge_bacon/qqplots_comparative/STR_Bacon_Scatter.pdf",width=5,height=5)
ggscatter(STR_FullTab_Bacon, 
        x = "log2FoldChange_original", 
        y = "log2FoldChange_corrected",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
   ) +
 ylim(-10,10)+
 xlim(-10,10)
dev.off()

pdf("str_dge_bacon/qqplots_comparative/STR_Bacon_QQ_pvals.pdf",width=5,height=3)
plot(bc, type="qq")
dev.off()


pdf("str_dge_bacon/qqplots_comparative/STR_Bacon_QQplots.pdf",width=8,height=5)
par(mfrow = c(1, 2))
qqnorm(beta / se, main="QQ Plot Before Correction", pch=20)
qqline(beta / se, col="red")
qqnorm(beta_corrected / se_corrected, main="QQ Plot After Correction", pch=20)
qqline(beta_corrected / se_corrected, col="blue")
dev.off()


# Upset Plot
STR_DGE_Original <- readxl::read_excel("str_dge_bacon/STR_DGE_Original.xlsx")
STR_DGE_Bacon <- readxl::read_excel("str_dge_bacon/STR_DGE_Bacon.xlsx")

tab1 <- STR_DGE_Original %>%
          mutate(Class = "Original") %>%
          select(Gene, Class)

tab2 <- STR_DGE_Bacon %>%
          mutate(Class = "Corrected") %>%
          select(Gene, Class)

tmp_upset <- rbind(tab1,tab2)

l <- split(as.character(tmp_upset$Gene),tmp_upset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))

pdf("str_dge_bacon/Upset_Plot_Intersection.pdf", width = 6, height = 4)
upset(fromList(l),,nsets = 4, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(Original = "#89E651", Corrected = "#DB8BE4"), 
    alpha = 0.5))))
dev.off()







