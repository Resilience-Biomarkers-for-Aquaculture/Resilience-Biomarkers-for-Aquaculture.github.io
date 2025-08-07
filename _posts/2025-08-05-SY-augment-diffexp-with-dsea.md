---
layout: post
title: Creating a GMT file to use with GSEA
author: Steve Yost
tags: NextFlow, differential abundance, GSEA
---

# Motivation
Based on successful nf-core/differentialabundance runs, we wished to augment them
with gene ontology information in order to examine relationships with biological functions
of the differentially expressed genes. nf-core/differentialabundance supports this, and
requires one or more gene matrix files. The present task was to create a GMT file for
C. Virginica, because none existed, and re-run nf-core/differentialabundance.

# Results
Here is the [differential abundance summary report](https://steveyost-seqera.s3.us-east-1.amazonaws.com/Cvirg_Pmarinus/deseq2/study4_genelength_D28_results/report/study.html) including GSEA results in the *Gene set analysis* section. The link will download the HTML file for local display. 
The GMT file used, and its precursor files, are in our Git repository [here](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/tree/main/data), all prefixed with *4706222.C_virginica*, based on the GOA file name from which they were derived.

# Method
Because apparently no GMT file existed for C. Virginica, I used [ChatGPT to give me the following recommendations](https://chatgpt.com/share/686d4d43-9b44-800d-97ee-02da85a7c7a5):

1. Download C. virginica GO GAF from UniProt or NCBI
1. Parse GAF to map genes → GO terms
1. Write out a standard tab-delimited GMT
I found a [.goa file for C. Virginica](
https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/4706222.C_virginica.goa).
and added the file to our Git data directory.

I used the following scripts (created mostly by ChatGPT) in the following order:
1. [goa_to_gmt.py](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/scripts/goa_to_gmt.py) to create a base GMT file from the downloaded GOA file.
Noting that *all* evidence codes were `IEA`, meaning Inferred from Electronic Annotation, I had it filter only `ND` (No biological data available) codes.
ChatGPT noted:
    > You can proceed with IEA-only gene sets, with the understanding that:
    > * The annotations are computational predictions, not experimentally validated.
    > * Enrichment results will reflect putative function, useful for hypothesis generation but not definitive.
    >
    > This is a common reality when working with non-model species like Crassostrea virginica.
1. [add_go_names_to_gmt.py](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/scripts/add_go_names_to_gmt.py) to add gene ontology terms. Output file was [4706222.C_virginica_named.gmt](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/4706222.C_virginica_named.gmt)
1. [filter_gmt_by_geneset_size.py](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/scripts/filter_gmt_by_geneset_size.py) to filter out gene sets < 10 (not enough signal) or > 500 (too broad), producing [4706222.C_virginica_named_filtered.gmt](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/4706222.C_virginica_named_filtered.gmt)
1. [filter_gmt_by_go_depth.py](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/scripts/filter_gmt_by_go_depth.py) to produce 4706222.C_virginica_named_filtered_depth3.gmt.
1. After running GSEA within an nf-core/differentialabundance pipeline, the GSEA section only lists the GO terms, like `GO:0031032`, not the function, so I used [enrich_gmt_with_go_namespaces.py](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/scripts/enrich_gmt_with_go_namespaces.py) to have the GO IDs include encoded gene function, producing 4706222.C_virginica_named_filtered_depth3_enrichednames.gmt. I used this GMT file in the differentialabundance pipeline to produce the results linked in the *Results* section above.



 


