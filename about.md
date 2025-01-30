---
layout: page
title: Resilience Biomarkers for Aquaculture
subtitle:
---

Jump to:  
- [Project Summary](#project-summary)
- [Omics Datasets](#omics-datasets)
- [Biomarker Database](#biomarker-database)
- [People](#people)

### Project Summary

Economic challenges imposed by climate change and disease on the aquaculture industry
necessitate advances for improved animal welfare and resiliency. Biomarkers associated with
environmental and disease resilience traits can be leveraged in breeding and management
strategies. However, their discovery has been limited in part by the complexity of molecular
systems and the cost of genomics tools used to understand them. Advances in computational
approaches including machine learning algorithms, together with the wealth of genomic data that
has amassed, enable powerful meta-analyses for improved biomarker discovery in aquaculture
species.

This project aims to advance the discovery and characterization of biomarkers
through mining publicly available shellfish genomic datasets from resilient populations.

**The objectives are to:**

1) Develop standardized open-access, user-friendly, reproducible bioinformatics pipelines for resilience biomarker discovery through systematic reanalysis, data integration and meta-analysis  

2) Build a user-friendly open-access comprehensive database of candidate resilience biomarkers that is widely available for use by the aquaculture community.

The resulting database will enable improved molecular tool development for more efficient
phenotype selection and health monitoring, implementation of selection methods that use a systems biology approach for simultaneous
improvement of multiple traits, and ultimately increased animal fitness and improved animal welfare.

View more project details [here](https://github.com/Resilience-Biomarkers-for-Aquaculture/Resilience-Biomarkers-for-Aquaculture.github.io/blob/master/docs/ProjectSummaryandNarrative.pdf)

![](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Resilience-Biomarkers-for-Aquaculture.github.io/master/img/fig1.png)
a) Graphical summary of Objective 1. b) Example pipelines and software for systematic reanalysis and data integration. c) Example of resilience biomarker database (modeled after Tamborero et al. 201846).

### Omics datasets processed

**Species**|**Data Type**|**stress class**|**stressor**|**Phenotype**|**Phenotpe Summary**|**Reference**|**DOI**|**SRA links**|**meta data file**|**counts table**|**DEG data**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
_C. gigas_|T|environment|thermal |thermotolerance different among lines|resilience| Arredondo-Espinoza _et al._ 2023|[https://doi.org/10.1016/j.cbd.2023.101089](https://doi.org/10.1016/j.cbd.2023.101089)|[https://www.ncbi.nlm.nih.gov/biosample?LinkName=bioproject_biosample_all&from_uid=516210](https://www.ncbi.nlm.nih.gov/biosample?LinkName=bioproject_biosample_all&from_uid=516210)|[SraRunTable.csv](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/refs/heads/main/data/SraRunTable.csv) | [salmon.merged.gene_counts.tsv](https://gannet.fish.washington.edu/emma.strand/rnaseq/Cgigas_ArredondoEspinoza2023/salmon.merged.gene_counts.tsv)|
_C. virginica_|T|disease|Perkinsus|infection tolerance|resilience|Proestou _et al._ 2023|[https://doi.org/10.3389/fgene.2023.1054558](https://doi.org/10.3389/fgene.2023.1054558)|[https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=5&WebEnv=MCID_679bbe17db915954f5c764fa&o=acc_s%3Aa](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=5&WebEnv=MCID_679bbe17db915954f5c764fa&o=acc_s%3Aa)|[SraRunTable.csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/SraRunTable.csv) | [salmon.merged.gene_counts.tsv](https://gannet.fish.washington.edu/emma.strand/rnaseq/Cvir_Prkns_rnaseq_dataset1/salmon.merged.gene_counts.tsv)|
_C. virginica_|E,T|disease|Perkinsus|infected|sensitivity|Johnson _et al._ 2020|[https://doi.org/10.3389/fmars.2020.00598](https://doi.org/10.3389/fmars.2020.00598)|[https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP246310&o=acc_s%3Aa](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP246310&o=acc_s%3Aa)|[SraRunTable(1).csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/SraRunTable%20(1).csv) |[salmon.merged.gene_counts.tsv](https://gannet.fish.washington.edu/emma.strand/rnaseq/Cvir_Prkns_rnaseq_dataset2/salmon.merged.gene_counts.tsv) |
_C. gigas & virginica_|T|disease|Perkinsus|infection tolerance|resilience|Chan _et al._ 2021|[10.3389/fgene.2021.795706](https://doi.org/10.3389/fgene.2021.795706)|[https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN11031730&o=acc_s%3Aa](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN11031730&o=acc_s%3Aa)|[SraRunTable(2).csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/SraRunTable%20(2).csv) | [salmon.merged.gene_counts.tsv](https://gannet.fish.washington.edu/emma.strand/rnaseq/Cvir_Prkns_rnaseq_dataset3/salmon.merged.gene_counts.tsv) *Cgigas data [here](https://gannet.fish.washington.edu/emma.strand/rnaseq/Cvir_Prkns_rnaseq_dataset3_Cgigas/salmon.merged.gene_counts.tsv)|
_C. virginica_|T|disease|Perkinsus|infection |sensitivity|Sullivan and Proestou 2021|[https://doi.org/10.1016/j.aquaculture.2021.736831](https://doi.org/10.1016/j.aquaculture.2021.736831)|[https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP301630&o=acc_s%3Aa](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP301630&o=acc_s%3Aa)|[SraRunTable(3).csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/SraRunTable%20(3).csv) |[salmon.merged.gene_counts.tsv](https://gannet.fish.washington.edu/emma.strand/rnaseq/Cvir_Prkns_rnaseq_dataset4/salmon.merged.gene_counts.tsv) |
_C. virginica_|T|disease|Perkinsus|infection |sensitivity|Proestou and Sullivan 2020|[https://doi.org/10.1016/j.fsi.2019.12.001](https://doi.org/10.1016/j.fsi.2019.12.001)|[https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=4&WebEnv=MCID_679bbe17db915954f5c764fa&o=acc_s%3Aa](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=4&WebEnv=MCID_679bbe17db915954f5c764fa&o=acc_s%3Aa)|[SraRunTable(4).csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/SraRunTable%20(4).csv) |[salmon.merged.gene_counts.tsv](https://gannet.fish.washington.edu/emma.strand/rnaseq/Cvir_Prkns_rnaseq_dataset5/salmon.merged.gene_counts.tsv) |
_C. virginica_|E|environment|OA ||resilience| Roberts _et al._ _unpublished_|[https://github.com/sr320/ceasmallr](https://github.com/sr320/ceasmallr)|[https://gannet.fish.washington.edu/seashell/bu-github/ceasmallr/data/](https://gannet.fish.washington.edu/seashell/bu-github/ceasmallr/data/)|[L18_larvae_meta.csv](https://github.com/sr320/ceasmallr/blob/main/data/L18_larvae_meta.csv) | [bismark.cov](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/bismark/methylation_calls/methylation_coverage/)|

### Biomarker Database

### People
Project Director:  
Shelly Wanamaker, PhD (she/her)  
Research Scientist  
[shelly.wanamaker@gmgi.org](mailto:shelly.wanamaker@gmgi.org)  
[https://github.com/shellywanamaker](https://github.com/shellywanamaker)

Project Personnel:  
Emma Strand, PhD (she/her)  
Postdoctoral Scientist  
[emma.strand@gmgi.org](mailto:emma.strand@gmgi.org)  
[https://github.com/emmastrand](https://github.com/emmastrand)

Steve Yost (he/him)  
Bioinformatics Software Engineer  
[steve.yost+gmgi@gmail.com](mailto:steve.yost@gmail.com)

 **Office and mailing address:**  
Gloucester Marine Genomics Institute  
417 Main Street  
Gloucester, MA 01930 USA  

**Funding:**  
This project is supported by the USDA Agriculture and Food Research Initiative Animal Health and Production and Animal Products program under the Animal Breeding, Genetics, and Genomics section award number 2024-67015-41794 to Dr. Shelly A. Wanamaker.
