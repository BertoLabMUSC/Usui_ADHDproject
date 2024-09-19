# WGCNA

# WGCNA
suppressPackageStartupMessages({
library(tidyverse)
library(WGCNA)
library(cluster)
library(ggplot2)
library(ggpubr)
library(igraph)
})
enableWGCNAThreads()


# Module Viz
dir.create("wgcna_output/modules_visualizations/")

load("futcounts/Expression_Input.RData")


mod <- read.table("wgcna_output/ModuleOutput.txt",header=T)



# Hubs 
hubs <- mod %>% 
            group_by(ModuleColor) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

hub_list <- split(hubs, hubs$ModuleColor)


tmp <- list()
adjMat <- list()
modColors <- list()
valueList <- list()
colorList <- list()
g1 <- list()
layoutFR <- list()

for(i in 1:length(hub_list)){
tmp[[i]] <- log2(exp[rownames(exp) %in% hub_list[[i]]$Gene,]+1)
adjMat[[i]] <- bicor(t(tmp[[i]]))
adjMat[[i]][abs(adjMat[[i]])<0.5] <- 0
adjMat[[i]] <- abs(adjMat[[i]])
adjMat[[i]] <- adjMat[[i]][match(hub_list[[i]]$Gene,rownames(adjMat[[i]])), match(hub_list[[i]]$Gene,colnames(adjMat[[i]]))]
modColors[[i]] <- as.matrix(t(hub_list[[i]]$ModuleColor))
valueList[[i]] <- lapply(1:ncol(modColors[[i]]), function(x) as.numeric(!is.na(modColors[[i]][,x])))
colorList[[i]] <- lapply(1:ncol(modColors[[i]]), function(x) modColors[[i]][,x])
g1[[i]] <- graph.adjacency(as.matrix(adjMat[[i]]),mode="undirected",weighted=TRUE,diag=FALSE)
g1[[i]] <- delete.vertices(g1[[i]],which(degree(g1[[i]])<1))
layoutFR[[i]] <- layout_with_lgl(g1[[i]],maxiter = 500)

ggsave(paste0("wgcna_output/modules_visualizations/",names(hub_list)[[i]], ".pdf"),
plot.igraph(g1[[i]],
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hub_list[[i]]$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR[[i]],
            edge.color=adjustcolor("grey", alpha.f = .3)),
height = 5, width = 5)
}


