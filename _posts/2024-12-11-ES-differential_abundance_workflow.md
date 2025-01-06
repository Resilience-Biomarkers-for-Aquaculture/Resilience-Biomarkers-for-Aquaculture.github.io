---
layout: post
title: Exploring nf-core differential abundance workflow
author: Emma Strand
tags: nextflow Klone differential-abundance
---

We want to test out this nf-core workflow: https://github.com/nf-core/differentialabundance for our RNAseq data output. 

![](https://raw.githubusercontent.com/nf-core/differentialabundance/1.5.0//docs/images/workflow.png)

### Required files 

Observations / Samplesheet input: `--input '[path to samplesheet file]'`

```
sample,fastq_1,fastq_2,condition,replicate,batch
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,control,1,A
CONTROL_REP2,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,control,2,B
```

RNAseq data:

```
--matrix 'salmon.merged.gene_counts.tsv' \
--transcript_length_matrix 'salmon.merged.gene_lengths.tsv'
```

Contrasts information: `--contrasts '[path to contrasts file]'`

```
id,variable,reference,target,blocking
condition_control_treated,condition,control,treated,
condition_control_treated_blockrep,condition,control,treated,replicate;batch
```

Feature annotations `--gtf '[path to gtf file]'`


### Script 

`rnaseq_diffab.sh`

```
#!/bin/bash
#SBATCH --account=srlab
#SBATCH --error=/gscratch/scrubbed/elstrand/Cvir_disease_meta/scripts/output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=/gscratch/scrubbed/elstrand/Cvir_disease_meta/scripts/output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --time=1-20:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=2

## Set paths 
samplesheet=""
contrasts=""
gene_counts=""
transcript_counts=""
output_dir=""
gtf=""

## Run differential abundance workflow 
nextflow run nf-core/differentialabundance -resume \
-profile rnaseq,singularity \
--input ${samplesheet} \
--contrasts ${contrasts} \
--matrix ${gene_counts} \
--transcript_length_matrix ${transcript_counts} \
--outdir ${output_dir} \
--gtf ${gtf}
```

To run: `sbatch rnaseq_diffab.sh`


#### Examples 

Riss ran this on thier rnaseq from sea urchins. They ended up with log2change values that seemed unrealistic so they ended up doing this on their own. I wonder if this was because the output from rnaseq was run with --deseq_vst so if the differential abundance workflow is also normalizing these value it would be a problem..?

```
#!/bin/bash
#SBATCH --error=output_messages/"%x_error.%j" #if job fails, error report is put into this file
#SBATCH --output=output_messages/"%x_output.%j" #final job report
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=nfcore_differentialabundance_controlsovertime_RMK
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

#Change directory to where data files are
cd /work/gmgi/Urchin/nfcoreRNAseq

#Load the 'singularity' container version 3.10.3
module load singularity/3.10.3

#Load the nextflow module, most recent version on Discovery
module load nextflow/24.04.4

#####NEXTFLOW DifferentialAbundance Pipeline w RNAseq Output#########
nextflow -log ./nextflow.log run nf-core/differentialabundance \
--max_cpus 24 \
--input /work/gmgi/Urchin/nfcoreRNAseq/controlsonly_samplesheet.csv \
--contrasts /work/gmgi/Urchin/nfcoreRNAseq/controlsonly_timecontrasts.csv \
--matrix /work/gmgi/Urchin/nfcoreRNAseq/nfcoreOutput/star_salmon/salmon.merged.gene_counts.tsv \
--transcript_length_matrix /work/gmgi/Urchin/nfcoreRNAseq/nfcoreOutput/star_salmon/salmon.merged.gene_lengths.tsv \
--outdir /work/gmgi/Urchin/nfcoreRNAseq/nfcoreDiffAbund_ControlsOverTimeOutput \
--gtf /work/gmgi/Urchin/GenomesAndAnnotationFiles/Spur_5.0.60.gtf.gz \
-profile rnaseq,singularity
```



