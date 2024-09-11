---
layout: post
title: RNAseq workflow with de novo transcriptome assembly
---

This is a description of the **nf-core/denovotranscript** pipeline for de novo transcriptome assembly of paired-end short reads from bulk RNA-seq. 

![nf-core/transfuse metro map](docs/images/denovotranscript_metro_map.drawio.svg)

1. Read QC of raw reads ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))
3. Read QC of trimmed reads ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Remove rRNA or mitochondrial DNA (optional) ([`SortMeRNA`](https://hpc.nih.gov/apps/sortmeRNA.html))
5. Transcriptome assembly using any combination of the following:

   - [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki) with normalised reads (default=True)
   	- normalization helps with the very non-uniform coverage of RNAseq datasets (e.g. highly vs lowly expressed genes, high-coverage in some regions may imply more sequencing errors complicating assembly down the road)
   	- [normalization will leave poorly covered regions unchanged, but will down-sample reads in high-coverage regions](https://biohpc.cornell.edu/lab/doc/trinity_workshop_part1.pdf)
   - [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki) with non-normalised reads
   - [`rnaSPAdes`](https://ablab.github.io/spades/rna.html) medium filtered transcripts outputted (default=True)
   - [`rnaSPAdes`](https://ablab.github.io/spades/rna.html) soft filtered transcripts outputted
   - [`rnaSPAdes`](https://ablab.github.io/spades/rna.html) hard filtered transcripts outputted
   	- differences between trinity and rnaSPAdes include processing power/time with rnaSPAdes requiring less, and rnaSPAdes can sometimes recover more transcripts (especially for low coverage genes) than trinity however this is dataset dependent [_Bushmanova et al. 2019_](https://academic.oup.com/gigascience/article/8/9/giz100/5559527)
   	- some [do's and don'ts for de novo transcriptome assembly](http://arthropods.eugenes.org/EvidentialGene/evigene/docs/perfect-mrna-assembly-2013jan.txt)

6. Redundancy reduction with [`Evidential Gene tr2aacds`](http://arthropods.eugenes.org/EvidentialGene/). A transcript to gene mapping is produced from Evidential Gene's outputs using [`gawk`](https://www.gnu.org/software/gawk/).
7. Assembly completeness QC ([`BUSCO`](https://busco.ezlab.org/))
8. Other assembly quality metrics ([`rnaQUAST`](https://github.com/ablab/rnaquast))
9. Transcriptome quality assessment with [`TransRate`](https://hibberdlab.com/transrate/), including the use of reads for assembly evaluation. This step is not performed if profile is set to `conda` or `mamba`.
10. Pseudo-alignment and quantification ([`Salmon`](https://combine-lab.github.io/salmon/))
11. HTML report for raw reads, trimmed reads, BUSCO, and Salmon ([`MultiQC`](http://multiqc.info/))


### 1. Set-up samplesheet
prepare a samplesheet with your input data (each row represents a pair of fastq files (paired end)) that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,AEG588A5_S5_L003_R2_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,AEG588A6_S6_L004_R2_001.fastq.gz
```

### 2. Set-up config file with parameter choices  
- run with rnaSPAdes using all 3 filters
- no quantification because will run through nf-core RNAseq pipeline

### 3. Create slurm script

### 4. Assess transcriptomes and decide best to use as reference in RNAseq pipeline with reference

Plan is to start with Day 30 RNAseq datasets from Roberto's 2023 thermotolerance study: [https://doi.org/10.1016/j.cbd.2023.101089](https://doi.org/10.1016/j.cbd.2023.101089)
----
****
