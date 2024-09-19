# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776, brain expressed = 15585, WGCNA list = 6029)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir str_dge/enrichments_str/

cp utils/geneset/*.RData str_dge/enrichments_str/
cp utils/Enrichment.r str_dge/enrichments_str/
cp str_dge/STR_Dge_HumanID.txt str_dge/enrichments_str/
cp str_dge/STR_Dge_MouseID.txt str_dge/enrichments_str/

cd str_dge/enrichments_str/

mkdir STATS/

# scRNA healty
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l Allen_MultiReg_CellMarkers_GeneSet.RData -p -b 15585 -o STATS/Allen_Markers -W 4 -H 4
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l Allen_MultiReg_CellMarkers_GeneSet_Simple.RData -p -b 19776 -o STATS/Allen_Markers_Simple -W 4 -H 3
Rscript Enrichment.r -g STR_Dge_MouseID.txt -l GeneSets_scMouse.RData -p -b 15585 -o STATS/scMouse -W 4 -H 5

# scRNA Disorders
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l ALZ_SingleCell_DEGs.RData -p -b 15585 -o STATS/ALZ_SingleCell -W 4 -H 5
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l ASD_SingleCell_DEGs.RData -p -b 15585 -o STATS/ASD_SingleCell -W 4 -H 5
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l GeneSets_BroadTrans_ASD_SingleCell.RData -p -b 15585 -o STATS/ASD_BroadTrans_SingleCell -W 4 -H 5
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l GeneSets_BroadTranscript_ASD.RData -p -b 15585 -o STATS/ASD_BroadTranscript -W 4 -H 5

# Neuropsy
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l PsychENCODE_DEGs.RData -p -b 15585 -o STATS/PSY_DEGS -W 4 -H 3
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l PsychEncode_Modules.RData -p -b 15585 -o STATS/PSY_MODS -W 4 -H 5
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l ASD_SFARI.RData -p -b 15585 -o STATS/ASD_Sfari -W 4 -H 2
Rscript Enrichment.r -g STR_Dge_HumanID.txt -l GeneSets_Disorders.RData -p -b 15585 -o STATS/Disorders -W 4 -H 4

rm *.RData

