# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776, brain expressed = 15585, WGCNA list = 6029)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir wgcna_output/enrichments_wgcna/

cp utils/geneset/*.RData wgcna_output/enrichments_wgcna/
cp utils/Enrichment_wgcna.r wgcna_output/enrichments_wgcna/
cp wgcna_output/ModuleOutput.txt wgcna_output/enrichments_wgcna/
cp wgcna_output/ModuleOutput_HumanID.txt wgcna_output/enrichments_wgcna/

cd wgcna_output/enrichments_wgcna/

mkdir STATS/

# scRNA healty
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l Allen_MultiReg_CellMarkers_GeneSet.RData -p -b 15585 -o STATS/Allen_Markers_wgcna -W 10 -H 8
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_scMouse.RData -p -b 15585 -o STATS/scMouse_wgcna -W 10 -H 6

# scRNA Disorders
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l ALZ_SingleCell_DEGs.RData -p -b 15585 -o STATS/ALZ_SingleCell_wgcna -W 10 -H 8
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l ASD_SingleCell_DEGs.RData -p -b 15585 -o STATS/ASD_SingleCell_wgcna -W 10 -H 8

# Neuropsy
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l PsychENCODE_DEGs.RData -p -b 15585 -o STATS/PSY_DEGS_wgcna -W 10 -H 6
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l PsychEncode_Modules.RData -p -b 15585 -o STATS/PSY_MODS_wgcna -W 10 -H 6
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l ASD_SFARI.RData -p -b 15585 -o STATS/ASD_Sfari_wgcna -W 10 -H 6

# DGE
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_ADHD_DGE.RData -p -b 15585 -o STATS/ADHD_DGE_wgnca -W 10 -H 8

rm *.RData

