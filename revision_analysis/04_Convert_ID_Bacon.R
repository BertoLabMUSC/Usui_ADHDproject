suppressPackageStartupMessages({
library(biomaRt)
library(tidyverse)
})

# Load Data
CTX_DGE <- readxl::read_excel("ctx_dge_bacon/CTX_DGE_Bacon.xlsx")

# Convert to human 
human = useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = CTX_DGE$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb <- merge(CTX_DGE,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

df <- signComb %>%
                mutate(Direction = case_when(log2FoldChange_corrected > 0 ~ "CTX_Upreg", log2FoldChange_corrected < 0  ~ "CTX_Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp <- data.frame(Gene = df$Gene, Direction = rep("CTX_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df,"ctx_dge_bacon/CTX_Dge_Bacon_HumanID.txt",sep="\t",quote=F,row.names=F)


df <- CTX_DGE %>%
                mutate(Direction = case_when(log2FoldChange_corrected > 0 ~ "CTX_Upreg", log2FoldChange_corrected < 0  ~ "CTX_Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp <- data.frame(Gene = df$Gene, Direction = rep("CTX_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df,"ctx_dge_bacon/CTX_Dge_Bacon_MouseID.txt",sep="\t",quote=F,row.names=F)


# STR
STR_DGE <- readxl::read_excel("str_dge_bacon/str_DGE_Bacon.xlsx")


MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = STR_DGE$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb <- merge(STR_DGE,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

df <- signComb %>%
                mutate(Direction = case_when(log2FoldChange_corrected > 0 ~ "STR_Upreg", log2FoldChange_corrected < 0  ~ "STR_Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp <- data.frame(Gene = df$Gene, Direction = rep("STR_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df,"str_dge_bacon/STR_Dge_Bacon_HumanID.txt",sep="\t",quote=F,row.names=F)


df <- STR_DGE %>%
                mutate(Direction = case_when(log2FoldChange_corrected > 0 ~ "STR_Upreg", log2FoldChange_corrected < 0  ~ "STR_Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp <- data.frame(Gene = df$Gene, Direction = rep("STR_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df,"str_dge_bacon/STR_Dge_Bacon_MouseID.txt",sep="\t",quote=F,row.names=F)
