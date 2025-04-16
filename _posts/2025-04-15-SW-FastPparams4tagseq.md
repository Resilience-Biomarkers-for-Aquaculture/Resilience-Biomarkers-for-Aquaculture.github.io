---
layout: post
title: Reprocess TAG-seq data with different FastP parameters 
author: Shelly Wanamaker
tags: TAG-seq RNA-seq fastp adapter
---

### This post is related to github issues [# 26](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/26) and [# 28](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/28)

The [multiQC report](https://gannet.fish.washington.edu/emma.strand/rnaseq/Cvir_Prkns_rnaseq_dataset2/pipeline_info/multiqc_report.html#fastqc_trimmed_fastqc_per_sequence_gc_content) for Johnson et al 2020 dataset 2 from [Emma's RNAseq analysis](https://resilience-biomarkers-for-aquaculture.github.io/ES-rnaseq_Cvig_Prkns/) showed a lot of adapter contamination as well as overrepresentation of polyA sequence. We realized this is because it is TAG-seq data, which should be trimmed of adapters and polyA sequence. Here is an example of how others do this: [Sam Gurr's pipeline for trimming TAGseq data](https://github.com/SamGurr/Pgenerosa_OA_TagSeq/blob/main/TagSeq/HPC_work/Geoduck_TagSeq_Bioinf.md#Trimming-polyA-tail)

The [Johnson et al. 2020](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2020.00598/full#h3) samples were run at the UT Austin genomics facility and here is their TAG-seq [library protocol](https://github.com/z0on/tag-based_RNAseq/blob/master/TagSeq_sample_prep_june2019.docx) and [preprocessing bioinformatics pipeline](https://cloud.wikis.utexas.edu/wiki/spaces/bioiteam/pages/47720861/TAG-Seq+Preprocessing). These resources are also on their Github page [https://github.com/z0on/tag-based_RNAseq](https://github.com/z0on/tag-based_RNAseq).

To test out different FastP parameters, I ran RNAseq many times on a subset of the data (50K reads/sample). First, I had to redownload the data and did that with the FetchNGS pipeline:

```
nextflow run \
nf-core/fetchngs \
-resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250407_FetchNGS/ids.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250407_FetchNGS \
--download_method sratools

```

Next I ran the following iterations of the RNAseq pipeline:

### USING `--trim_poly_x` and `--adapter_sequence` 

Adapter sequences were copied from the [fastP readme] (https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters)

```

nextflow run nf-core/rnaseq -resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/samplesheet_rnaseq_dataset2.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250407_RNAseq \
--gtf /gscratch/srlab/elstrand/genomes/C_virginica/mod_GCF_002022765.2_C_virginica-3.0_genomic.gtf \
--gff /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz \
--fasta /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz \
--trimmer fastp \
--extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10 --reads_to_process 50000 --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' \
--aligner star_salmon \
--skip_pseudo_alignment \
--multiqc_title Cvir_disease_rnaseq_dataset2 \
--deseq2_vst

#copy data over to gannet
rsync --progress --verbose --archive --exclude adapter_fasta/ --exclude star_salmon/ --exclude fastqc/ --exclude work/ 20250407_RNAseq shellytrigg@gannet.fish.washington.edu:/volume2/web/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/
```
multiqc output here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/multiqc/star_salmon/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/multiqc/star_salmon/multiqc_report.html)



### USING `--adapter_fasta`

Adapter fasta file used: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/polyA.fa](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/polyA.fa)

```
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina TruSeq Adapter Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>polyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

The sequences in the adapter fasta file were copied from the [FastP Github readme](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters)

```
nextflow run nf-core/rnaseq -resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/samplesheet_rnaseq_dataset2.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/adapter_fasta \
--gtf /gscratch/srlab/elstrand/genomes/C_virginica/mod_GCF_002022765.2_C_virginica-3.0_genomic.gtf \
--gff /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz \
--fasta /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz \
--trimmer fastp \
--extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10 --reads_to_process 50000 --adapter_fasta /mmfs1/gscratch/scrubbed/strigg/analyses/20250407_RNAseq/adapter_fasta/polyA.fa' \
--aligner star_salmon \
--skip_pseudo_alignment \
--multiqc_title Cvir_disease_rnaseq_dataset2 \
--deseq2_vst
```

multiqc report here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/multiqc/star_salmon/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/multiqc/star_salmon/multiqc_report.html)

### USING `--adapter_fasta`  with UT Austin sequence

Adapter fasta file used: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/UT_tag_adapt/polyA.fa](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/UT_tag_adapt/polyA.fa)

```
>Illumina TruSeq Adapter Read 1
AGATCGGAAG
>polyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

The adapter sequence was copied from the UT Austin TAGseq perl script here: [https://github.com/z0on/tag-based_RNAseq/blob/master/tagseq_trim_launch.pl](https://github.com/z0on/tag-based_RNAseq/blob/master/tagseq_trim_launch.pl)


```
nextflow run nf-core/rnaseq -resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/samplesheet_rnaseq_dataset2.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/adapter_fasta/UT_tag_adapt \
--gtf /gscratch/srlab/elstrand/genomes/C_virginica/mod_GCF_002022765.2_C_virginica-3.0_genomic.gtf \
--gff /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz \
--fasta /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz \
--trimmer fastp \
--extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10 --reads_to_process 50000 --adapter_fasta /mmfs1/gscratch/scrubbed/strigg/analyses/20250407_RNAseq/adapter_fasta/UT_tag_adapt/polyA.fa' \
--aligner star_salmon \
--skip_pseudo_alignment \
--multiqc_title Cvir_disease_rnaseq_dataset2 \
--deseq2_vst
```

multiqc report here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/UT_tag_adapt/multiqc/star_salmon/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/UT_tag_adapt/multiqc/star_salmon/multiqc_report.html)

### USING `--adapter_fasta` with UT Austin 3ILL and 5ILL oligos

Adapter fasta file used: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/3ILL-30_5ILL/polyA.fa](hhttps://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/3ILL-30_5ILL/polyA.fa) and looks like the following:

```
>polyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

>3ILL-polyT
TTAGATCGGAAGAGCACACGT

>5ILL_RC
AGATCGGAAGAGCGTCGTGTAG

```

The adapter sequences were copied from the UT Austin [TAGseq oligo spreadsheet](https://github.com/z0on/tag-based_RNAseq/blob/master/tagseq_oligos_order_192samples.xls) on [github](https://github.com/z0on/tag-based_RNAseq).

```
nextflow run nf-core/rnaseq -resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/samplesheet_rnaseq_dataset2.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/adapter_fasta/3ILL-30_5ILL \
--gtf /gscratch/srlab/elstrand/genomes/C_virginica/mod_GCF_002022765.2_C_virginica-3.0_genomic.gtf \
--gff /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz \
--fasta /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz \
--trimmer fastp \
--extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10 --reads_to_process 50000 --adapter_fasta /mmfs1/gscratch/scrubbed/strigg/analyses/20250407_RNAseq/adapter_fasta/3ILL-30_5ILL/polyA.fa' \
--aligner star_salmon \
--skip_pseudo_alignment \
--multiqc_title Cvir_disease_rnaseq_dataset2 \
--deseq2_vst
```

multiqc report is here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/3ILL-30_5ILL/multiqc/star_salmon/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/3ILL-30_5ILL/multiqc/star_salmon/multiqc_report.html)

### USING `adapter_fasta` with UT Austin 3ILL and 5ILL and a couple other seqs 
Adapter fasta file used: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/3ILL_5ILL_sw/polyA.fa](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/3ILL_5ILL_sw/polyA.fa) and looks like the following:

```
>3ILL-polyT_RC
TTAGATCGGAAGAGCACACGT
>5ILL_RC
AGATCGGAAGAGCGTCGTGTAG
>5ILL_switch
ACCCCATGGGGCTACACGACGCTCTTCCGATCT
>IC2-P7_RC
TCGTATGCCGTCTTCTGCTTG
>polyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

```

These sequences were also taken from the UT Austin [TAGseq oligo spreadsheet](https://github.com/z0on/tag-based_RNAseq/blob/master/tagseq_oligos_order_192samples.xls) on [github](https://github.com/z0on/tag-based_RNAseq).


```
nextflow run nf-core/rnaseq -resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/samplesheet_rnaseq_dataset2.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250407_RNAseq/adapter_fasta/3ILL_5ILL_sw \
--gtf /gscratch/srlab/elstrand/genomes/C_virginica/mod_GCF_002022765.2_C_virginica-3.0_genomic.gtf \
--gff /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz \
--fasta /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz \
--trimmer fastp \
--extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10 --reads_to_process 50000 --adapter_fasta /mmfs1/gscratch/scrubbed/strigg/analyses/20250407_RNAseq/adapter_fasta/3ILL_5ILL_sw/polyA.fa' \
--aligner star_salmon \
--skip_pseudo_alignment \
--multiqc_title Cvir_disease_rnaseq_dataset2 \
--deseq2_vst

```
multiqc report is here:[https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/3ILL_5ILL_sw/multiqc/star_salmon/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250407_RNAseq/adapter_fasta/3ILL_5ILL_sw/multiqc/star_salmon/multiqc_report.html)

## Summary of Results

| Method| Adapter Content   |  GC Content |
:--|:-------------------------:|:-------------------------:
default FastP | <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.43.01%20PM.png?raw=true" width="600" height="300"> | <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.47.25%20PM.png?raw=true" width="600" height="300"> 
`--trim_poly_x` and `--adapter_sequence`| <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.43.27%20PM.png?raw=true" width="600" height="300"> |  <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.47.40%20PM.png?raw=true" width="600" height="300"> 
`--adapter_fasta` from [FastP](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters)| <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.45.20%20PM.png?raw=true" width="600" height="300"> |  <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.47.53%20PM.png?raw=true" width="600" height="300">
`--adapter_fasta`  with UT Austin pipeline sequence| <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.45.52%20PM.png?raw=true" width="600" height="300"> |  <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.48.04%20PM.png?raw=true" width="600" height="300">
`--adapter_fasta` with UT Austin 3ILL and 5ILL oligo seqs |<img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.46.03%20PM.png?raw=true" width="600" height="300">  | <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.48.12%20PM.png?raw=true" width="600" height="300">  
`adapter_fasta` with UT Austin 3ILL and 5ILL and a couple other seqs |  <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.46.11%20PM.png?raw=true" width="600" height="300"> |  <img src="https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/images/CvirgDataSet2FastP_screenshots/Screen%20Shot%202025-04-11%20at%203.48.21%20PM.png?raw=true" width="600" height="300">

### Rerunning the pipeline on the whole dataset

I went with the last option and reran the RNAseq analysis on the whole dataset

```
nextflow run nf-core/rnaseq -resume \
-c /gscratch/srlab/nextflow/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250415_RNAseq/samplesheet_rnaseq_dataset2.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250415_RNAseq \
--gtf /gscratch/srlab/elstrand/genomes/C_virginica/mod_GCF_002022765.2_C_virginica-3.0_genomic.gtf \
--gff /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz \
--fasta /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz \
--trimmer fastp \
--extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10 --adapter_fasta /mmfs1/gscratch/scrubbed/strigg/analyses/20250415_RNAseq/polyA.fa' \
--aligner star_salmon \
--skip_pseudo_alignment \
--multiqc_title Cvir_disease_rnaseq_dataset2 \
--deseq2_vst

```
multiqc here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250415_RNAseq/multiqc/star_salmon/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250415_RNAseq/multiqc/star_salmon/multiqc_report.html#fastqc_trimmed_fastqc_adapter_content)

It seems there is still an overrepresented adapter sequence `TGTGCTCTTCCGATCTAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT `


####Attempt 2

I created a new adapter fasta containing the primers in both orientations to try to account for this and included polyT sequence. I also removed the switch primer sequences because I don't think 

```
>3ILL-polyT_RC
TTAGATCGGAAGAGCACACGT
>5ILL_RC
AGATCGGAAGAGCGTCGTGTAG
>3ILL-polyT
ACGTGTGCTCTTCCGATCTAA
>5ILL
CTACACGACGCTCTTCCGATCT
>polyT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>polyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

And reran the pipeline as follows:

```
nextflow run nf-core/rnaseq -resume \
-c /gscratch/srlab/nextflow/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250415_RNAseq_2/samplesheet_rnaseq_dataset2.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250415_RNAseq_2 \
--gtf /gscratch/srlab/elstrand/genomes/C_virginica/mod_GCF_002022765.2_C_virginica-3.0_genomic.gtf \
--gff /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz \
--fasta /gscratch/srlab/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz \
--trimmer fastp --extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10 --adapter_fasta /mmfs1/gscratch/scrubbed/strigg/analyses/20250415_RNAseq_2/polyA.fa' \
--aligner star_salmon \
--skip_pseudo_alignment \
--multiqc_title Cvir_disease_rnaseq_dataset2 \
--deseq2_vst
```
### Output

multiqc report is here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250415_RNAseq_2/mutliqc/star_salmon/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250415_RNAseq_2/mutliqc/star_salmon/multiqc_report.html)

counts matrices are here: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250415_RNAseq_2/](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250415_RNAseq_2/)

### Next steps
re-run integrated analysis with newly aligned data