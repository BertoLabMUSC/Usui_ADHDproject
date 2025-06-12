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
library(qqman)
library(bacon)

})

dir.create("revision_plots")

# QQ plots for CTX
dge <- readxl::read_excel("ctx_dge/CTX_DGE_FullStats.xlsx")

dge <- dge %>% 
        mutate(qq = -2*log(padj))

pdf("revision_plots/CTX_QQplot.pdf",width=4,height=4)
qq(dge$padj)
dev.off()

median(qchisq(dge$padj, df = 2, lower.tail = FALSE)) / qchisq(0.5, df = 2)
#0.4386712

median(qchisq(dge$pvalue, df = 2, lower.tail = FALSE)) / qchisq(0.5, df = 2)
#1.438351


# QQ plots for STR
dge <- readxl::read_excel("str_dge/STR_DGE_FullStats.xlsx")

dge <- dge %>% 
        mutate(qq = -2*log(padj))

pdf("revision_plots/STR_QQplot.pdf",width=4,height=4)
qq(dge$padj)
dev.off()


median(qchisq(dge$padj, df = 2, lower.tail = FALSE)) / qchisq(0.5, df = 2)
#0.5378166


median(qchisq(dge$pvalue, df = 2, lower.tail = FALSE)) / qchisq(0.5, df = 2)
#1.537688



# BACON analsis STR
dge <- readxl::read_excel("str_dge/STR_DGE_FullStats.xlsx")

beta <- dge$log2FoldChange
se <- dge$lfcSE

bc <- bacon(effectsize = beta, standarderror = se)


beta_corrected <- es(bc)
se_corrected <- se(bc)
pval_corrected <- pval(bc)

corrected_dge <- data.frame(
  gene = dge$Gene,
  log2FoldChange_original = beta,
  lfcSE_original = se,
  pvalue_original = dge$pvalue,
  padj_original = dge$padj,
  log2FoldChange_corrected = beta_corrected,
  lfcSE_corrected = se_corrected,
  pvalue_corrected = pval_corrected,
  padj_corrected = p.adjust(pval_corrected,"BH")
)


original <- corrected_dge %>%
                  mutate(Abs = abs(log2FoldChange_original)) %>%
                  filter(padj_original < 0.05 & Abs > 0.3) %>%
                  arrange(desc(Abs))

corrected <- corrected_dge %>%
                  mutate(Abs = abs(log2FoldChange_corrected)) %>%
                  filter(padj_corrected < 0.05 & Abs > 0.3) %>%
                  arrange(desc(Abs))


openxlsx::write.xlsx(corrected_dge, 
                     file = "revision_plots/STR_DGE_FullStats_BACON.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(original, 
                     file = "revision_plots/STR_DGE_BACON_Original.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(corrected, 
                     file = "revision_plots/STR_DGE_BACON_Corrected.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")



pdf("revision_plots/STR_Bacon_Scatter.pdf",width=5,height=5)
ggscatter(corrected_dge, 
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


pdf("revision_plots/STR_Bacon_QQ_pvals.pdf",width=5,height=3)
plot(bc, type="qq")
dev.off()

pdf("revision_plots/STR_Bacon_QQplots.pdf",width=8,height=5)
par(mfrow = c(1, 2))
qqnorm(beta / se, main="QQ Plot Before Correction", pch=20)
qqline(beta / se, col="red")
qqnorm(beta_corrected / se_corrected, main="QQ Plot After Correction", pch=20)
qqline(beta_corrected / se_corrected, col="blue")
dev.off()


# BACON analsis CTX
dge <- readxl::read_excel("ctx_dge/CTX_DGE_FullStats.xlsx")

beta <- dge$log2FoldChange
se <- dge$lfcSE

bc <- bacon(effectsize = beta, standarderror = se)


beta_corrected <- es(bc)
se_corrected <- se(bc)
pval_corrected <- pval(bc)

corrected_dge <- data.frame(
  gene = dge$Gene,
  log2FoldChange_original = beta,
  lfcSE_original = se,
  pvalue_original = dge$pvalue,
  padj_original = dge$padj,
  log2FoldChange_corrected = beta_corrected,
  lfcSE_corrected = se_corrected,
  pvalue_corrected = pval_corrected,
  padj_corrected = p.adjust(pval_corrected,"BH")
)


original <- corrected_dge %>%
                  mutate(Abs = abs(log2FoldChange_original)) %>%
                  filter(padj_original < 0.05 & Abs > 0.3) %>%
                  arrange(desc(Abs))

corrected <- corrected_dge %>%
                  mutate(Abs = abs(log2FoldChange_corrected)) %>%
                  filter(padj_corrected < 0.05 & Abs > 0.3) %>%
                  arrange(desc(Abs))


openxlsx::write.xlsx(corrected_dge, 
                     file = "revision_plots/CTX_DGE_FullStats_BACON.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(original, 
                     file = "revision_plots/CTX_DGE_BACON_Original.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(corrected, 
                     file = "revision_plots/CTX_DGE_BACON_Corrected.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")



pdf("revision_plots/CTX_Bacon_Scatter.pdf",width=5,height=5)
ggscatter(corrected_dge, 
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

pdf("revision_plots/CTX_Bacon_QQ_pvals.pdf",width=5,height=3)
plot(bc, type="qq")
dev.off()

pdf("revision_plots/CTX_Bacon_Scatter_OnlyOriginal.pdf",width=5,height=5)
ggscatter(original, 
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


pdf("revision_plots/CTX_Bacon_QQplots.pdf",width=8,height=5)
par(mfrow = c(1, 2))
qqnorm(beta / se, main="QQ Plot Before Correction", pch=20)
qqline(beta / se, col="red")
qqnorm(beta_corrected / se_corrected, main="QQ Plot After Correction", pch=20)
qqline(beta_corrected / se_corrected, col="blue")
dev.off()





