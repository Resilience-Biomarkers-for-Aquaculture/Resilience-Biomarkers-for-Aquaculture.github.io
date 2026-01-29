---
layout: post
title: RNASeq workflow with reference genome
author: Emma Strand
---

## RNASeq workflow construction with reference genome

Author(s): [E. Strand](https://github.com/emmastrand) 

Overview of program options: 

01. Quality Assessment via FastQC and MultiQC           
02. Data pre-processing: Trim reads via FastP or TrimGalore!       
- [FastP](https://github.com/OpenGene/fastp): can be faster and more efficient;  is an all-in-one preprocessor that combines quality control, adapter trimming, read filtering, and base correction;  can automatically detect adapter sequences for single-end and paired-end Illumina data; generates single html file for QC        
- [TrimGalore!](https://github.com/FelixKrueger/TrimGalore): focuses primarily on adapter trimming and quality filtering; requires manual input of adapter sequences or uses a set of predefined Illumina adapter sequences     
03. Option 1: Genome alignment via STAR or HISAT2 OR Option 2: Pseudo-alignment via Salmon or Kallisto    
- [STAR](https://github.com/alexdobin/STAR): potentially higher alignment rates than HISAT2, potentially better at mapping longer transcripts in low-quality genomes, better for draft genomes  
- [HISAT2](https://daehwankimlab.github.io/hisat2/): generally faster than STAR, fewer computational resources, better at handling SNPs when made aware of them  
- [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html): equivalent to Kallisto, salmon can run multithread so potentially faster than Kallisto 
- [Kallisto](https://github.com/pachterlab/kallisto): equivalent to Salmon in quality  
04. Quantification: UMI-tools dedup prior to RSEM or Salmon or no quantification with HISAT2    
- [RSEM](https://github.com/deweylab/RSEM): estimating gene and isoform expression levels from RNA-Seq data; For most targeted RNA-seq analyses in research labs, either Salmon or RSEM could be equally effective choices; salmon is slightly faster than RSEM, with no well-annotated transcriptome then STAR+RSEM might be preferred   
- [UMI-tools](https://github.com/CGATOxford/UMI-tools): dedup = groups PCR duplicates and deduplicates reads to yield one read per group  
05. Post-processing: Sort alignments via SAMtools, deduplicate via picard, transcript assembly and quantification via [StringTie](https://ccb.jhu.edu/software/stringtie/)     
06. Final Quality Assessment via [Preseq](http://smithlabresearch.org/software/preseq/), [Qualimap](http://qualimap.bioinfo.cipf.es/), [dupRadar](https://bioconductor.org/packages/release/bioc/html/dupRadar.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [RSeQC](http://rseqc.sourceforge.net/)

This pipeline utilizes nf-core's [nextflow](https://training.nextflow.io/basic_training/) workflow for RNAseq data: https://nf-co.re/rnaseq/3.15.0

![](https://raw.githubusercontent.com/nf-core/rnaseq/3.15.0//docs/images/nf-core-rnaseq_metro_map_grey_animated.svg)

## Set-up: Samplesheet 

Prior to running nf-core pipeline, build a samplesheet.csv file that contains sample ID and the paths to R1 and R2 files. Recommend to create file on RStudio Interactive on Discovery Cluster using (`create_metadatasheets.R`).

Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Rows with the same sample identifier are considered technical replicates and merged automatically. The strandedness refers to the library preparation and will be automatically inferred if set to `auto`.

```
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,auto
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,auto
```

or 

```
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,forward
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,,reverse
```

- sampleID (required): Unique sample IDs, must start with a letter, and can only contain letters, numbers or underscores (no hyphons!).  
- forwardReads (required): Paths to (forward) reads zipped FastQ files  
- reverseReads (optional): Paths to reverse reads zipped FastQ files, required if the data is paired-end  
- run (optional): If the data was produced by multiple sequencing runs, any string  

*This is an R script, not slurm script. Open RStudio interactive on Discovery Cluster to run this script.* 

Example of R script file from another project:

```
## Creating sample sheet for offshore wind eDNA project 

### Step 1: In terminal 

# cd raw_data 
# ls * > ../metadata/sample_list.txt

### Resume steps below

library(dplyr)
library(stringr)
library(strex)
#library(filesstrings)

### Read in sample sheet 

sample_list <- read.delim2("/work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/rawdata", header=F) %>% 
  dplyr::rename(forwardReads = V1) %>%
  mutate(sampleID = str_after_nth(forwardReads, "data/", 1),
         sampleID = str_before_nth(sampleID, "_R", 1),
         sampleID = ifelse(!grepl('July', sampleID), sub("_.*", "", sampleID), sampleID)
         )

# creating sample ID 
sample_list$sampleID <- gsub("-", "_", sample_list$sampleID)

# keeping only rows with R1
sample_list <- filter(sample_list, grepl("R1", forwardReads, ignore.case = TRUE))

# duplicating column 
sample_list$reverseReads <- sample_list$forwardReads

# replacing R1 with R2 in only one column 
sample_list$reverseReads <- gsub("R1", "R2", sample_list$reverseReads)

# rearranging columns 
sample_list <- sample_list[,c(2,1,3)]

sample_list %>% write.csv("/work/gmgi/Fisheries/eDNA/offshore_wind2023/metadata/samplesheet.csv", 
                          row.names=FALSE, quote = FALSE)

```

## Parameter choices

Data pre-processing: fastp instead of TrimGalore!  
Genome alignment: STAR instead of HISAT2   
Pseudo alignment: Salmon instead of Kallisto   
Quantification: Salmon instead of RSEM 

Alignment options: https://nf-co.re/rnaseq/3.15.0/docs/usage/#alignment-options

## Slurm scripts

Check for updated pipeline: `nextflow pull nf-core/rnaseq`

Module loads and environments TBD depending on the server that we do this on.. UW? Is this right for UW?  
Fill in extra args for fastp, star, and salmon depending on dataset we're using? 

`rnaseq_refgenome.sh`:

```
#!/bin/bash
#SBATCH --error=extract_output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=extract_output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --job-name=rnaseq_pseudoalign
#SBATCH --mem=100GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

module load singularity/3.10.3
module load nextflow/23.10.1

# SET PATHS 
metadata="" 
output_dir=""
genome=""
gtf_file=""
gff_file=""

nextflow run nf-core/rnaseq \
    -resume \
    -profile singularity \
    --input ${metadata}/samplesheet.csv \
    --outdir ${output_dir} \
    --gtf ${gtf_file} \
    --gff ${gff_file} \
    --fasta ${genome} \
    --trimmer fastp \
    --aligner star_salmon \
    --pseudo_aligner salmon \
    --multiqc_title rnaseq_projectitle_refgenome \
    --extra_fastp_args '--cut_mean_quality 30 --trim_front1 10' \
    --extra_star_align_args '' \
    --extra_salmon_quant_args '' \
    --deseq2_vst \
    --multiqc_methods_description

```


To run: `sbatch rnaseq_refgenome.sh`  

Note that the pipeline will generate the following files in the working directory:

```
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```
