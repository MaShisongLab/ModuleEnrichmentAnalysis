# ModuleEnrichmentAnalysis

**An algorithm for performing GO and MP term enrichment analysis in mouse**

**ModuleEnrichmentAnalysis** is an algorithm to perform enrichment analysis for mouse gene sets. The GO term annotations were obtained from Ensembl BioMart, while the MP term annotations were obtained from MGI ([Mouse Genome Informatics](https://www.informatics.jax.org/downloads/reports/index.html)). 

## Table of Contents
- [Install](#Install)
- [Usage](#Usage)
- [References](#References)

## Install
This algorithm works in R. Copy the file `ModuleEnrichmentAnalysis.R` and the subfolder `data` that contains the annotation files to your R working folder and start using it.


## Usage

The algorithm takes a list of clusters (saved in the file `clusters.txt`) and a list of background genes (saved in the file `back_ground_genes.txt`) as input. The file `clusters.txt` should contain two columns separated by tab, with first column being gene id, and the second column being cluster id. The file `back_ground_genes.txt` should contain a list of back ground genes used for the analysis. Only Ensembl gene IDs are supported.

The functions are:

<b>go = go_enrichment_analysis (`cluster_file = 'clusters.txt'`, `bk_gene_file = 'back_ground_genes.txt'`, `go_annotation_file = 'data/Mouse.GO.annotation.txt'`, `go_name_file = 'data/GO.names.txt'`, `gene_symbl_file = 'data/Mouse.gene.symbl.txt'`)</b>

`cluster_file` - the file saving the clusters information
</br>`back_ground_genes.txt` - the file with a list of back ground genes
</br>`go_annotation_file` - the file containing GO annotations for mouse genes
</br>`go_name_file` - the file containing GO names
</br>`gene_symbl_file` - the file containing gene symbls for mouse genes

<b>mp = mp_enrichment_analysis (`cluster_file = 'clusters.txt'`, `bk_gene_file = 'back_ground_genes.txt'`, `mp_annotation_file = 'data/Mouse.MP.annotation.txt'`, `mp_name_file = 'data/MP.names.txt'`, `gene_symbl_file = 'data/Mouse.gene.symbl.txt'`)</b>

`cluster_file` - the file saving the clusters information
</br>`back_ground_genes.txt` - the file with a list of back ground genes
</br>`mp_annotation_file` - the file containing MP annotations for mouse genes
</br>`mp_name_file` - the file containing MP names
</br>`gene_symbl_file` - the file containing gene symbls for mouse genes



```R
# R code
# make sure the file 'ModuleEnrichmentAnalysis.R' and the 'data' folder is in your R working folder.
source('ModuleEnrichmentAnalysis.R')
go = go_enrichment_analysis('clusters.txt')
mp = mp_enrichment_analysis('clusters.txt')
write.csv(go,'go.results.csv')
write.csv(mp,'mp.results.csv')
```
The results are save to `go.results.csv` and `mp.results.csv`.

## References

To be added.

