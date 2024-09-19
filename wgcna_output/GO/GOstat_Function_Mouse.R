suppressPackageStartupMessages(library(GOstats))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(multtest))
suppressPackageStartupMessages(library(tidyverse))

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
wgcna = list.files(pattern = 'ModuleOutput*')
mod=read.table(wgcna,sep="\t",header=T)

mod <- mod %>%
       mutate(ModuleColor = as.factor(ModuleColor), ModuleName = as.factor(ModuleName))


uni=read.table("Mouse_Exp_Genes.txt",header=T)

listBP=list()
dfBP=list()
listMF=list()
dfMF=list()
listCC=list()
dfCC=list()
for (net in levels(mod$ModuleName)) {
listBP[[net]]=mod[mod$ModuleName == net,]
listMF[[net]]=mod[mod$ModuleName == net,]
listCC[[net]]=mod[mod$ModuleName == net,]
dfBP[[net]]=dogo(listBP[[net]]$Gene,uni$hgnc_symbol,species="mouse", goP = 0.05,cond=FALSE, ontology = "BP")
dfBP[[net]]$adj=p.adjust(dfBP[[net]]$Pvalue,"BH")
dfBP[[net]]=dfBP[[net]][dfBP[[net]]$adj < 0.05,]
dfMF[[net]]=dogo(listMF[[net]]$Gene,uni$hgnc_symbol,species="mouse", goP = 0.05,cond=FALSE, ontology = "MF")
dfMF[[net]]$adj=p.adjust(dfMF[[net]]$Pvalue,"BH")
dfMF[[net]]=dfMF[[net]][dfMF[[net]]$adj < 0.05,]
dfCC[[net]]=dogo(listCC[[net]]$Gene,uni$hgnc_symbol,species="mouse", goP = 0.05,cond=FALSE, ontology = "CC")
dfCC[[net]]$adj=p.adjust(dfCC[[net]]$Pvalue,"BH")
dfCC[[net]]=dfCC[[net]][dfCC[[net]]$adj < 0.05,]
}

dir.create("GO_BP/")
dir.create("GO_MF/")
dir.create("GO_CC/")

net=levels(mod$ModuleColor)
for (i in 1:length(dfBP)) write.table(dfBP[[i]], file = paste("GO_BP/MODULE_GO_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)
for (i in 1:length(dfMF)) write.table(dfMF[[i]], file = paste("GO_MF/MODULE_GO_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)
for (i in 1:length(dfCC)) write.table(dfCC[[i]], file = paste("GO_CC/MODULE_GO_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)



