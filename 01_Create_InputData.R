library(tidyverse)
library(data.table)
library(future.apply)
plan(multicore)

setwd("futcounts")

files = list.files(pattern = 'Primary_Gene_Assigned_exon.txt')
samples <- gsub("Primary_Gene_Assigned_exon.txt", "", files )
myfiles = lapply(files, read.table,sep="\t",fill=TRUE,header=T)

tmp <- Reduce(dplyr::full_join, myfiles) %>%
		 as.data.frame()

colnames(tmp) <-  sub(".UNIQUE.SORTED.bam", "", colnames(tmp))
colnames(tmp) <-  sub("X", "", colnames(tmp))

exp <- tmp %>%
			select(-Chr,-Start,-End,-Strand,-Length,-Status) %>%
			group_by(Geneid) %>%
			summarise_all(sum) %>% 
			na.omit() %>%
			column_to_rownames("Geneid")

Length <- tmp %>%
			select(Geneid,Length) %>%
			group_by(Geneid) %>%
			summarise_all(max) %>% 
			na.omit() %>%			
			column_to_rownames("Geneid")

cpm <- future_apply(exp, 2, function(x) x/sum(as.numeric(x)) * 10^6)

# Calculate the RPKM
rpkm <- future_apply(exp, 2, function(x) 10^9 * x / as.numeric(Length$Length) / sum(as.numeric(x)))

# Calculate the TPM
rpk <- future_apply(exp, 2, function(x) x/(as.numeric(Length$Length)/1000))
tpm <- future_apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

save(exp,cpm,rpkm,tpm, file = "Expression_Input.RData")

# Filter by region
exp_ctx <- exp[,grep("CT_",colnames(exp))] 
exp_str <- exp[,grep("STR_",colnames(exp))] 

rpkm_ctx <- rpkm[,grep("CT_",colnames(rpkm))] 
rpkm_str <- rpkm[,grep("STR_",colnames(rpkm))] 

save(exp_ctx,rpkm_ctx, file = "Expression_Input_CTX.RData")
save(exp_str,rpkm_str, file = "Expression_Input_STR.RData")





