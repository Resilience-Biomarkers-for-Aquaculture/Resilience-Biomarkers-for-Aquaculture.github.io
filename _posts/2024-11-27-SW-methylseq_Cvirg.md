---
layout: post
title: Run methylseq on Klone
tags: nextflow Klone denovotranscript
---
### Purpose
Related to this GitHub issue: https://github.com/RobertsLab/resources/issues/2045#issuecomment-2521028376

### Procedure

1. Create samplesheet
I created the samplesheet linked below to use in the methylseq pipeline as [described here](https://nf-co.re/methylseq/2.7.1/docs/usage). I got the paths on Klone (https://github.com/RobertsLab/resources/issues/2045#issuecomment-2521028376).

https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/samplesheet.csv

2. Run
I ran the following commands to start the pipeline but ran into a disk quote exceeded error that was resolved as described here:

```
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 15 --mem=250GB --time=2-12:00:00


mamba activate nextflow


nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20241127_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241127_methylseq \
--fasta /gscratch/scrubbed/sr320/github/ceasmallr/data/genome/Cvirginica_v300.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html
```
second attempt still errored
```
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=100GB --time=4-12:00:00

nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20241127_methylseq_test/samplesheet_test.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241127_methylseq_test \
--fasta /gscratch/scrubbed/sr320/github/ceasmallr/data/genome/Cvirginica_v300.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html
```
3rd attempt using test samplesheet.csv below still failed
```
#samplesheet.csv
sample,fastq_1,fastq_2
SRR389222_sub1,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz
SRR389222_sub2,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub2.fastq.gz
SRR389222_sub2,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub3.fastq.gz
Ecoli_10K_methylated,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R2.fastq.gz


salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=150GB --time=2-12:00:00

nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input samplesheet.csv \
--outdir . \
--genome GRCh37

```
# reinstalling mambaforge in srlab.

This was initially installed in my home dir and I want to limit my home dir storage because that is what caused my disk quota exceeded error. The main culprit was an old work directory from a previous nextflow run that was accidentally written there, but I still prefer to use /gscratch/srlab to store this stuff because it has more capacity.

Ran the following lines of code from `/gscratch/srlab/strigg/bin/`
```
wget https://github.com/conda-forge/miniforge/releases/download/24.9.0-0/Mambaforge-24.9.0-0-Linux-x86_64.sh
bash Mambaforge-24.9.0-0-Linux-x86_64.sh -b -p /gscratch/srlab/strigg/bin/mambaforge
micromamba run -p  /gscratch/srlab/strigg/bin/mambaforge mamba create -n nextflow -c bioconda -c conda-forge python=3.8 nf-core nextflow -y -q
```
ran the following to test functionality
```
1735  /gscratch/srlab/strigg/bin/mambaforge/bin/mamba activate nextflow
1736  /gscratch/srlab/strigg/bin/mambaforge/bin/mamba
```
these worked and the disk quota error was resolved


#### started at 08:50 EST Dec 3, 2024:

```
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=100GB --time=4-12:00:00

mamba activate nextflow

nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20241127_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241127_methylseq \
--fasta /gscratch/scrubbed/sr320/github/ceasmallr/data/genome/Cvirginica_v300.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html
```

I had to quit this job because there were 12 scripts that would not run on our node because it was going to conflict with scheduled maintenance of Klone. see issue: https://github.com/RobertsLab/resources/issues/2045#issuecomment-2521028376. After I modified the .config file I re-ran the above code and the pipeline completed in ~21 hours including the time it took for the previously completed steps. This is pretty amazing considering this dataset is ~469 GB.

### Results

multiqc report of results and qc for entire pipeline: https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/multiqc/bismark/multiqc_report.html

nextflow pipeline report:
https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/nf_report.html

nextflow pipeline timeline:
https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/nf_timeline.html

the bismark output is here: https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/bismark/
I only copied over the .cov files, but the alignment files, deduplicated files, methylation calls (e.g. CHH_OT_CF01-CM02-Larvae_1_val_1_bismark_bt2_pe.deduplicated.txt.gz), and m-bias files and other outputs are still on /gscratch/scrubbed/strigg/analyses/20241127
