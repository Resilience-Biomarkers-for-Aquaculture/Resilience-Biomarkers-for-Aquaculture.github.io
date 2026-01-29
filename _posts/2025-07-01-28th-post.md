---
layout: post
title: Run nf-core differentialabundance pipeline on _C.virg_ data
author: Shelly Wanamaker
tags: RNAseq differentialabundance nf-core virginica perkinsus
---

### Main Goals: 
- Understand how Diff. Abundance pipeline works
- Understand the differences in immune responses that lead to tolerant and sensitive phenotypes
	- compare fed vs. injected animals from study 4 (Steve)
	- compare breeds in dosed animals from study 1 (Shelly)	

For Study 1 [Prosteu et al. 2023](https://doi.org/10.3389/fgene.2023.1054558), RNAseq analysis included:

1. Preliminary RNAseq analysis: PCA of VST transformed data
	- Separation of treatment groups 	
2. WGCNA to identify modules related to traits
3. DESeq on dose 10^8 samples only
	- threshold: transcripts must have >=10 counts across all samples
	- comparison of family variable only 
	- GLM shrinkage applied to Log2 fold change estimates
	- BH FDR correction

Question to answer: 	

- Does breed impact expression under dose-independent infection?
- Dosed vs. control 
- Ideally we'd want to test: expression ~ Breed + treatment + Breed:treatment 
- Since testing for interactive effects is not straight forward in this pipeline (), we'll select only dosed animals and compare breeds
- another option to try that we may be able to extract the effects from breed:
	id,variable,reference,target,blocking
	treatment_control_dosed_blockrep,treatment,Control,Dose,breed

### Prepare input files for pipeline
**GTF**
The GTF had some errors that Emma [found previously](https://github.com/Resilience-Biomarkers-for-Aquaculture/Resilience-Biomarkers-for-Aquaculture.github.io/blob/master/_posts/2024-12-11-ES-differential_abundance_workflow.md#script). One problem was empty gene_ids and one problem was omission of 'transcript' in the 3rd field [(see github issue #25)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/25). A workaround for transcript not being in the 3rd field is changing the feature type that the pipeline uses from transcript to gene using the `--features_gtf_feature_type` parameter. 

Correcting the GTF for empty gene_ids: 

```
# check the gene_id fields in entire file
zcat GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz | awk -F"\t" '{print $9}'| awk -F";" '{print $1}' | awk '{print $1, substr($2,1,4)}' | sort | uniq -c | less
      5  
     26 gene_id ""
      4 gene_id "ATP
     12 gene_id "COX
      5 gene_id "CYT
1530292 gene_id "LOC
      4 gene_id "ND1
      4 gene_id "ND2
      4 gene_id "ND3
      8 gene_id "ND4
      4 gene_id "ND5
      4 gene_id "ND6
   1175 gene_id "Trn
```
only 26 lines have empty gene_ids   
   
```
#check the lines that don't have 'gene_id'
zcat GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz | grep -v 'gene_id'

#gtf-version 2.2
#!genome-build C_virginica-3.0
#!genome-build-accession NCBI_Assembly:GCF_002022765.2
#!annotation-source NCBI Crassostrea virginica Annotation Release 100
###
```
it's only the header that doesn't have gene_id in the line.
For lines without gene_id, create one using the position
```
# make gene ID the position
   
zcat GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz | awk -F'\t' '{if ($9 ~ /gene_id ""/) {$9 = gensub(/gene_id ""/, "gene_id \"LOC_" $4 "-" $5 "\"", 1, $9);}OFS="\t"; print}' \

#check example that previously had 'gene_id ""'
| grep 'pos:15923..15925'

NC_007175.2	RefSeq	exon	15892	15957	.	+	.	gene_id "LOC_15892-15957"; transcript_id "unknown_transcript_1"; anticodon "(pos:15923..15925)"; gbkey "tRNA"; product "tRNA-Ala"; exon_number "1";


#that code works so save new GTF file
cat /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gtf |  awk -F'\t' '{if ($9 ~ /gene_id ""/) {$9 = gensub(/gene_id ""/, "gene_id \"LOC_" $4 "-" $5 "\"", 1, $9);}OFS="\t"; print}'  > /gscratch/srlab/strigg/GENOMES/GCF_002022765.2_C_virginica-3.0_genomic_noEmptyGeneIDs.gtf

#check new GTF doesn't have empty gene_id fields
grep 'gene_id ""' /gscratch/srlab/strigg/GENOMES/GCF_002022765.2_C_virginica-3.0_genomic_noEmptyGeneIDs.gtf

```

**Sampleheet**
We're modifying the samplesheet used in RNAseq pipeline because we need to add the meta data that the contrasts file will use to the samplesheet data

I subsetted the samples and excluded the controls for study 1. 

[samplesheet.csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differential_abundance_sheets/study1_samplesheet.csv)


Here is Steve's samplesheet for example: [rnaseq_difabundance_study4_s3_samplesheet.csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differential_abundance_sheets/rnaseq_difabundance_study4_s3_samplesheet.csv). There are repeated columns in this file (Collection_Interval_Days, Collection_Date, Treatment)

**Contrasts file**

[contrasts.csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differential_abundance_sheets/study1_contrasts.csv)

```
id,variable,reference,target
breed_dose_120_84,BREED,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017084
breed_dose_90_84,BREED,ABC_VIMS_Family_2017090,ABC_VIMS_Family_2017084
breed_dose_90_89,BREED,ABC_VIMS_Family_2017090,ABC_VIMS_Family_2017089
breed_dose_120_89,BREED,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017089
breed_dose_120_90,BREED,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017090
breed_dose_89_84,BREED,ABC_VIMS_Family_2017089,ABC_VIMS_Family_2017084
```

Here are some examples:

- Emma's  

```
id,variable,reference,target
condition_test_experiment,condition,test,experiment
```
-Steve's
```
id,variable,reference,target,blocking
D6_treatment_control_fed,treatment,Control,Fed,batch
D6_treatment_control_inj,treatment,Control,Injected,batch
```

**Gene counts matrix**
Steve used merged gene counts from all studies
https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/rnaseq_gene_counts/merged_gene_counts.tsv

Emma's
salmon.merged.gene_counts_length_scaled.tsv

### Running the pipeline
- Our matrix is merged_gene_counts so we should be using `gene` for the feature type. So we have to set `features_gtf_feature_type` because the default is transcript 


```
screen -S diffabund
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=24:00:00
mamba activate nextflow

#download gene counts matrix
wget https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/rnaseq_gene_counts/salmon.merged.gene_counts_length_scaled_dataset1.tsv

#download sample contrasts
wget https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/differential_abundance_sheets/study1_contrasts.csv

#download samplesheet
wget https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/differential_abundance_sheets/study1_samplesheet.csv

nextflow run nf-core/differentialabundance \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
-resume \
--study_name dataset_1 \
--input study1_samplesheet.csv \
--contrasts study1_contrasts.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250701_diffabund \
--matrix /gscratch/scrubbed/strigg/analyses/20250701_diffabund/salmon.merged.gene_counts_length_scaled_dataset1.tsv \
--feature_type gene \
--features_gtf_feature_type gene \
--gtf /gscratch/srlab/strigg/GENOMES/GCF_002022765.2_C_virginica-3.0_genomic_noEmptyGeneIDs.gtf 
```

- pipeline results report here: [20250701_diffabund/report](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250701_diffabund/report/dataset_1.html#Differential_analysis)


## attempt with thresholding
```
nextflow run nf-core/differentialabundance \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
-resume \
--study_name dataset_1 \
--input study1_samplesheet.csv \
--contrasts study1_contrasts.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250701_diffabund \
--matrix /gscratch/scrubbed/strigg/analyses/20250701_diffabund/salmon.merged.gene_counts_length_scaled_dataset1.tsv \
--feature_type gene \
--features_gtf_feature_type gene \
--filtering_min_abundance 10 \
--filtering_min_proportion 0.75 \
--gtf /gscratch/srlab/strigg/GENOMES/GCF_002022765.2_C_virginica-3.0_genomic_noEmptyGeneIDs.gtf 


```

- pipeline results report here: [20250701_diffabund/read_thresh10_minSamp075/](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250701_diffabund/read_thresh10_minSamp075/report/dataset_1.html#Results)

- complete results here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250701_diffabund/read_thresh10_minSamp075/](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250701_diffabund/read_thresh10_minSamp075/)


## Attempt combining resistant and sensitive families as one  

84 and 89 = resistant; 90 and 120 = sensitive
made a trait column in the samplesheet that indicates resistant or sensitive

contrasts file: 

```
id,variable,reference,target
trait_dose,trait,sensitive,resistant
```

run pipeline:

```
nextflow run nf-core/differentialabundance \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
-resume \
--study_name dataset_1 \
--input study1_samplesheet.csv \
--contrasts study1_contrasts.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250702_diffabund \
--matrix /gscratch/scrubbed/strigg/analyses/20250701_diffabund/salmon.merged.gene_counts_length_scaled_dataset1.tsv \
--feature_type gene \
--features_gtf_feature_type gene \
--filtering_min_abundance 10 \
--filtering_min_proportion 0.75 \
--gtf /gscratch/srlab/strigg/GENOMES/GCF_002022765.2_C_virginica-3.0_genomic_noEmptyGeneIDs.gtf 
```
- pipeline results report here: [20250701_diffabund/read_thresh10_minSamp075/](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250702_diffabund/report/dataset_1.html)

- complete results here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250702_diffabund/](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250702_diffabund/)



## Pipeline Errors

```

 Error in checkListIsSubset(val, samples[[var]], "contrast levels", "sample metadata variable") :
    Not all contrast levels (ABC_VIMS_Family_2017190) are available in the sample metadata variabl
e (ABC_VIMS_Family_2017084,ABC_VIMS_Family_2017089,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017084
,ABC_VIMS_Family_2017084,ABC_VIMS_Family_2017089,ABC_VIMS_Family_2017089,ABC_VIMS_Family_2017090,A
BC_VIMS_Family_2017090,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017089,ABC_VIMS_Family_2017084,ABC
_VIMS_Family_2017089,ABC_VIMS_Family_2017090,ABC_VIMS_Family_2017090,ABC_VIMS_Family_2017120,ABC_V
IMS_Family_2017084,ABC_VIMS_Family_2017084,ABC_VIMS_Family_2017089,ABC_VIMS_Family_2017090,ABC_VIM
S_Family_2017120,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017084,ABC_VIMS_Family_2017089,ABC_VIMS_
Family_2017089,ABC_VIMS_Family_2017090,ABC_VIMS_Family_2017084,ABC_VIMS_Family_2017090,ABC_VIMS_Fa
mily_2017120,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017084,ABC_VIMS_Family_2017084,ABC_VIMS_Fami
ly_2017089,ABC_VIMS_Family_2017090,ABC_VIMS_Family_2017090,ABC_VIMS_Family_2017120,ABC_VIMS_Family
_2017089,ABC_VI
  Calls: validate_inputs -> read_contrasts -> checkListIsSubset
  Execution halted

```
This error came from an error in my contrast file which had ABC_VIMS_Family_2017190 instead of ABC_VIMS_Family_2017090. I corrected the contrasts file and reran the pipeline with the resume option

I got another error of file name collision and I think it's because my IDs in my contrast file were all the same so i made them distinct and ran it again.

old contrasts file:
id,variable,reference,target
breed_dose,BREED,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017084
breed_dose,BREED,ABC_VIMS_Family_2017190,ABC_VIMS_Family_2017084
breed_dose,BREED,ABC_VIMS_Family_2017190,ABC_VIMS_Family_2017089
breed_dose,BREED,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017089
breed_dose,BREED,ABC_VIMS_Family_2017120,ABC_VIMS_Family_2017190
breed_dose,BREED,ABC_VIMS_Family_2017089,ABC_VIMS_Family_2017084

-------------------------------------------------------------------------
#checking fields in GTF file
zcat GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz | awk -F'\t' '{print NF}' | sort | uniq -c
      5 1
1531542 9

#checking fields in GFF file
zcat GCF_002022765.2_C_virginica-3.0_genomic.gff.gz | awk -F'\t' '{print NF}' | sort | uniq -c


zcat GCF_963853765.1_xbMagGiga1.1_genomic.gtf.gz | awk -F'\t' '{print $3,$4,$5,$9}' | awk -F';' '{print $1}' | awk '{if(($1~/gene/)||($1~/transcript/)) print $0}' | awk '{if (NR % 2 == 1) {v1 = $2; v2 = $3} else {if ($2 != v1 || $3 != v2) {print "Mismatch at line " NR ": Expected", v1, v2, "but got", $2, $3}}}END{if (NR % 2 == 1) {print "Warning: Last line has no pair"} else {print "All pairs matched"}}'



-----------------

Dataset 4: 2021 Proestou

Condition: Fed/Inj/Control
Batch: sampling 
Time: days


expression ~ condition + time + condition:time + 1|batch

different lines in contrast file for each different comparison. For example:

- For day 6: 
	- expression ~ condition + batch
	
id,variable,reference,target,blocking
D6_treatment_control_fed_blockrep,treatment,Control,Fed,replicate;batch
D6_treatment_control_inj_blockrep,treatment,Control,Injected,replicate;batch


id,variable,reference,target,blocking
D28_treatment_control_fed_block,treatment,Control,Fed,batch
D28_treatment_control_inj_block,treatment,Control,Injected,batch
D28_treatment_control_fed_block,treatment,Fed,Injected, batch

id,variable,reference,target,blocking

