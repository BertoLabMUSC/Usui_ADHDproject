suppressPackageStartupMessages({
library(biomaRt)
library(tidyverse)
})

# Load Data
mod <- read.table("wgcna_output/ModuleOutput.txt",header=T)

# Convert to human 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mod$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb <- merge(mod,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

df <- signComb %>%
                dplyr::select(HGNC.symbol,ModuleColor,kWithin, ModuleName) %>%
                dplyr::rename(Gene = HGNC.symbol) %>%
                arrange(ModuleColor)

write.table(df,"wgcna_output/ModuleOutput_HumanID.txt",sep="\t",quote=F,row.names=F)
