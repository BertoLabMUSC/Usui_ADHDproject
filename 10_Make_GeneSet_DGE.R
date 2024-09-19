suppressPackageStartupMessages({
library(tidyverse)
})

ctx <- read.table("ctx_dge/CTX_Dge_MouseID.txt",header=T)
str <- read.table("str_dge/STR_Dge_MouseID.txt",header=T)

df <- rbind(ctx,str)
df$Direction <- as.factor(df$Direction)
GeneSets <- split(df,df$Direction)

save(GeneSets, file = "utils/geneset/GeneSets_ADHD_DGE.RData")