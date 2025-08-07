---
layout: post
title: Run methylseq on Mchilensis WGBS data
author: Shelly Wanamaker
tags: methylseq nf-core klone
---

## Purpose
This post is related to [https://github.com/RobertsLab/resources/issues/2251](https://github.com/RobertsLab/resources/issues/2251)

Run Chilean mytilus methylation samples via nextflow methylseq.

Files needed:

- zipped fastqs: [https://gannet.fish.washington.edu/v1_web/Raw_WGBS_Mch/](https://gannet.fish.washington.edu/v1_web/Raw_WGBS_Mch/)
- genome: [https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/project-chilean-mussel/data/Mchi/](https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/project-chilean-mussel/data/Mchi/)
- config file: 
	- [uw\_hyak_srlab.config](https://gannet.fish.washington.edu/metacarcinus/Mchilensis/20250731_methylseq/uw_hyak_srlab.config)
	- This file is specific for this analysis using alignment score of -0.8 and attempting to perform the analysis first on coenv node, then on srlab node, then on ckpt-all because we already know it takes ~6 hours to align each fastq. 
- samplesheet: [samplesheet.csv](https://gannet.fish.washington.edu/metacarcinus/Mchilensis/20250731_methylseq/samplesheet.csv)
	- the format of this samplesheet is standard for nf-core pipelines and is a 3 column csv file that lists sample name, path of fastq1, and path of fastq2. An example of the first two lines is shown below. The format including the header must be exact. 

```
sample,fastq_1,fastq_2
LCo_BSr1,/gscratch/scrubbed/strigg/analyses/20250731_methylseq/raw-reads/LCo_BSr1_R1.fastq.gz,/gscratch/scrubbed/strigg/analyses/20250731_methylseq/raw-reads/LCo_BSr1_R2.fastq.gz

```

## Methods

Pipeline run on 8/5/2025

```
# start screen session

screen -S methylseq0801

# start interactive node
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=30GB --time=96:00:00
mamba activate nextflow

# zip fastq files in /gscratch/scrubbed/strigg/analyses/20250731_methylseq

gzip *.fastq

# run methylseq pipeline from /gscratch/scrubbed/strigg/analyses/20250731_methylseq

nextflow run nf-core/methylseq \
-c /gscratch/scrubbed/strigg/analyses/20250731_methylseq/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250731_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250731_methylseq \
--fasta /gscratch/scrubbed/strigg/analyses/20250731_methylseq/MchilensisGenomeV1.fa \
--accel \
-resume \
-with-report nf_report.html \
-with-trace \
--nomeseq 

rsync --archive --verbose --progress --exclude=work/ --exclude=*.bam --exclude=BismarkInde
x/ --exclude=*.fa --exclude=raw-reads/ 20250731_methylseq shellytrigg@gannet.fish.washington.edu:/volume2/web/metacarcinus/Mchilensis

```

## Results
-  multiqc report here: [multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/Mchilensis/20250731_methylseq/multiqc/bismark/multiqc_report.html)
- all pipeline output is here: [https://gannet.fish.washington.edu/metacarcinus/Mchilensis/20250731_methylseq/](https://gannet.fish.washington.edu/metacarcinus/Mchilensis/20250731_methylseq/)
	- bismark results: [bismark](https://gannet.fish.washington.edu/metacarcinus/Mchilensis/20250731_methylseq/bismark)

the trimming looks like it could be improved for all read 2 based on the multiqc.	

## Step-by-step instructions for re-running with improved trimming


If you don't have mamba in your path, run: 

```
/gscratch/srlab/nextflow/bin/miniforge/bin/mamba init
```

Then run:

``` 
source ~/.bashrc
```

create a working directory and change into it

```
# this is just an example, and you'd replace what's in '< >' with your username and the date

mkdir /gscratch/scrubbed/<username>/analyses/<date>_methylseq

cd /gscratch/scrubbed/<username>/analyses/<date>_methylseq
```

create a screen session

```
screen -S nextflow
```

request an interactive node. This run should complete in about 24 hours, but requesting 48 hours to be safe. We don't need a lot of memory here hence 30GB.

```
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=30GB --time=48:00:00
```

activate the nextflow environment 

```
mamba activate nextflow
```

run the methylseq nf-core pipeline with the improved trimming parameter. NOTE: you need to change the outdir to your working directory that you created above. You will also have to change the paths of the `.config` file, the samplesheet.csv, and the genome file if you run this after 8/21/2025 (because of the file shelflife on scrubbed). 

```
nextflow run nf-core/methylseq \
-c /gscratch/scrubbed/strigg/analyses/20250731_methylseq/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250731_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/<username>/analyses/<date>_methylseq \
--fasta /gscratch/scrubbed/strigg/analyses/20250731_methylseq/MchilensisGenomeV1.fa \
--accel \
--clip_r2 10 \
-resume \
-with-trace \
-with-timeline nf_timeline.html \
--nomeseq 
```

