suppressPackageStartupMessages(library(GOstats))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(multtest))

dogo <- function(names,universe,species="human", goP = 0.05, 
	cond=FALSE, ontology = "BP"){
    if(species=="human"){
		golib="org.Hs.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Hs.egSYMBOL2EG
  } else  if (species == "mouse") {
		golib="org.Mm.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Mm.egSYMBOL2EG
	}
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA))
 Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

# How to run
#load("NoStress_dge_output/NoStress_Dge_Data.RData")
mod=read.table("ctx_dge/CTX_Dge_HumanID.txt",sep="\t",header=T)
mod$Gene <- as.character(mod$Gene)
mod$Direction <- as.factor(mod$Direction)

uni <- read.table("utils/gencode.v24.annotation.protcod.length.genebody.txt")
uni_genes <- rownames(uni)

listBP=list()
dfBP=list()
listMF=list()
dfMF=list()
listCC=list()
dfCC=list()
for (net in levels(mod$Direction)) {
listBP[[net]]=mod[mod$Direction == net,]
listMF[[net]]=mod[mod$Direction == net,]
listCC[[net]]=mod[mod$Direction == net,]
dfBP[[net]]=dogo(listBP[[net]]$Gene,uni_genes,species="human", goP = 1,cond=FALSE, ontology = "BP")
dfBP[[net]]$padj=p.adjust(dfBP[[net]]$Pvalue,method="BH")
dfBP[[net]]=dfBP[[net]][dfBP[[net]]$padj < 0.2,]
dfMF[[net]]=dogo(listMF[[net]]$Gene,uni_genes,species="human", goP = 1,cond=FALSE, ontology = "MF")
dfMF[[net]]$padj=p.adjust(dfMF[[net]]$Pvalue,method="BH")
dfMF[[net]]=dfMF[[net]][dfMF[[net]]$padj < 0.2,]
dfCC[[net]]=dogo(listCC[[net]]$Gene,uni_genes,species="human", goP = 1,cond=FALSE, ontology = "CC")
dfCC[[net]]$padj=p.adjust(dfCC[[net]]$Pvalue,method="BH")
dfCC[[net]]=dfCC[[net]][dfCC[[net]]$padj < 0.2,]
}

dir.create("ctx_dge/GO_CTX/")
dir.create("ctx_dge/GO_CTX/GO_BP/")
dir.create("ctx_dge/GO_CTX/GO_MF/")
dir.create("ctx_dge/GO_CTX/GO_CC/")

net=levels(mod$Direction)
for (i in 1:length(dfBP)) write.table(dfBP[[i]], file = paste("ctx_dge/GO_CTX/GO_BP/BP_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)
for (i in 1:length(dfMF)) write.table(dfMF[[i]], file = paste("ctx_dge/GO_CTX/GO_MF/MF_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)
for (i in 1:length(dfCC)) write.table(dfCC[[i]], file = paste("ctx_dge/GO_CTX/GO_CC/CC_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)

# Create Plots
setwd("ctx_dge/GO_CTX/")

library(ggpubr)
library(tidyverse)

files = list.files(path="GO_BP",pattern = '*.txt',full.names = TRUE)
names <- gsub( "GO_BP/|BP_|.txt", "", files )
GeneSets = lapply(files, read.table,header=T,sep="\t")
names(GeneSets) <- names


filt <- vector("list", length = length(GeneSets))
names(filt) <- names(GeneSets)
class <- names(GeneSets)
for (i in 1:length(GeneSets))
    {
      filt[[i]] <- GeneSets[[i]] %>% 
              filter(Count > 5 & 
              Pvalue < 0.05) %>%
              select(Term,Pvalue,OddsRatio) %>% 
              mutate(log = -log10(Pvalue)) %>%
              as.data.frame()
  }

for (i in 1:length(GeneSets))
    {
  filt[[i]]$Class <- factor(rep(class[i],nrow(filt[[i]])))
  }

df <- do.call(rbind,filt)


top_labelled <- tbl_df(df) %>% 
                  group_by(Class) %>% 
                  top_n(n = 5, wt = abs(log))


pdf("BP_CTX_enrichment.pdf",width=8,height=3,useDingbats=FALSE)
ggscatter(top_labelled, x = "log", y = "OddsRatio",
   color = "Class", palette = c("royalblue", "royalblue", "royalblue"),size = 2,
   label = "Term", repel = TRUE,font.label = c(8, "plain"))+
xlab("-log10(p-value)")+ 
ylab("Odds Ratio") +
facet_wrap(~Class,ncol=3,nrow=1,scales="free")+
theme_classic()+
theme(legend.position="none")

dev.off()


files = list.files(path="GO_MF",pattern = '*.txt',full.names = TRUE)
names <- gsub( "GO_MF/|MF_|.txt", "", files )
GeneSets = lapply(files, read.table,header=T,sep="\t")
names(GeneSets) <- names


filt <- vector("list", length = length(GeneSets))
names(filt) <- names(GeneSets)
class <- names(GeneSets)
for (i in 1:length(GeneSets))
    {
      filt[[i]] <- GeneSets[[i]] %>% 
              filter(Count > 5 & 
              Pvalue < 0.05) %>%
              select(Term,Pvalue,OddsRatio) %>% 
              mutate(log = -log10(Pvalue)) %>%
              as.data.frame()
  }

for (i in 1:length(GeneSets))
    {
  filt[[i]]$Class <- factor(rep(class[i],nrow(filt[[i]])))
  }

df <- do.call(rbind,filt)


top_labelled <- tbl_df(df) %>% 
                  group_by(Class) %>% 
                  top_n(n = 5, wt = abs(log))


pdf("MF_CTX_enrichment.pdf",width=8,height=3,useDingbats=FALSE)
ggscatter(top_labelled, x = "log", y = "OddsRatio",
   color = "Class", palette = c("royalblue", "royalblue", "royalblue"),size = 2,
   label = "Term", repel = TRUE,font.label = c(8, "plain"))+
xlab("-log10(p-value)")+ 
ylab("Odds Ratio") +
facet_wrap(~Class,ncol=3,nrow=1,scales="free")+
theme_classic()+
theme(legend.position="none")

dev.off()


files = list.files(path="GO_CC",pattern = '*.txt',full.names = TRUE)
names <- gsub( "GO_CC/|CC_|.txt", "", files )
GeneSets = lapply(files, read.table,header=T,sep="\t")
names(GeneSets) <- names


filt <- vector("list", length = length(GeneSets))
names(filt) <- names(GeneSets)
class <- names(GeneSets)
for (i in 1:length(GeneSets))
    {
      filt[[i]] <- GeneSets[[i]] %>% 
              filter(Count > 5 & 
              Pvalue < 0.05) %>%
              select(Term,Pvalue,OddsRatio) %>% 
              mutate(log = -log10(Pvalue)) %>%
              as.data.frame()
  }

for (i in 1:length(GeneSets))
    {
  filt[[i]]$Class <- factor(rep(class[i],nrow(filt[[i]])))
  }

df <- do.call(rbind,filt)


top_labelled <- tbl_df(df) %>% 
                  group_by(Class) %>% 
                  top_n(n = 5, wt = abs(log))


pdf("CC_CTX_enrichment.pdf",width=8,height=3,useDingbats=FALSE)
ggscatter(top_labelled, x = "log", y = "OddsRatio",
   color = "Class", palette = c("royalblue", "royalblue", "royalblue"),size = 2,
   label = "Term", repel = TRUE,font.label = c(8, "plain"))+
xlab("-log10(p-value)")+ 
ylab("Odds Ratio") +
facet_wrap(~Class,ncol=3,nrow=1,scales="free")+
theme_classic()+
theme(legend.position="none")
dev.off()

