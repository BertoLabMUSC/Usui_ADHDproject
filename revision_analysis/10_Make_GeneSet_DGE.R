suppressPackageStartupMessages({
library(tidyverse)
})

ctx <- read.table("ctx_dge_bacon/CTX_Dge_Bacon_MouseID.txt",header=T)
str <- read.table("str_dge_bacon/STR_Dge_Bacon_MouseID.txt",header=T)

df <- rbind(ctx,str)
df$Direction <- as.factor(df$Direction)
GeneSets <- split(df,df$Direction)

save(GeneSets, file = "utils/geneset/GeneSets_ADHD_DGE_Bacon.RData")