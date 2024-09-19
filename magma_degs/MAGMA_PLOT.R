library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(WGCNA)
library(tidyverse)
files=list.files(pattern="STATISTICS")
tmp=as.data.frame(lapply(files,read.table,sep="\t",header=T))

df <- tmp[-grep("z",tmp$Sample),]
#df <- melt(tmp)
df$log <- -log10(df$P)

df$log[df$log < 1.3] <- NA


df <- df %>% separate(VARIABLE, c("Cluster","Class"), "_",remove=FALSE)

pdf("MAGMA_BUBBLE.pdf",width=5,height=3.5)
ggscatter(df, 
                        x = "Sample",
                        y = "VARIABLE",
                        size="log",
                        color="BETA",
                        alpha = 0.8,
                        xlab = "",ylab = "",) +
                        theme_minimal() + 
                        gradient_color(c("red", "darkred"))+
                        rotate_x_text(angle = 45)+
        coord_flip()+
        scale_size(range = c(2, 8))

dev.off()


