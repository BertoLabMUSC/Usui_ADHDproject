ADHD project 
==========================

This repository contains analysis code for the ADHD mouse model gene expression project carried out by researchers at the [Shimada Lab, MUSC](https://www.med.osaka-u.ac.jp/eng/introduction/research-5/anatomy/neuroscienceandcell) and [Berto Lab, MUSC](https://bertolab.org/).

## Cite this

If you use anything in this repository please cite the following publication:

URL: 


## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`futcounts`](futcounts/) | Input data of the initial processing and quality check. | 01_Create_InputData.R|
| [`ctx_dge`](ctx_dge/) | Output of the Differential expression analysis for the CTX. | 02_Dge_CTX.R \ 05_Enrichments_CTX.sh \ 07_GO_CTX.R|
| [`str_dge`](str_dge/) | Output of the Differential expression analysis for the STR. | 02_Dge_STR.R \ 05_Enrichments_STR.sh \ 07_GO_STR.R|
| [`wgcna_output`](wgcna_output/) | Output data of the coexpression analysis. | 09_WGCNA.R \ 11_Module_Visiualization.R|
| [`utils`](utils/) | Utility functions and data. | 08_Enrichment_Modules.sh |
