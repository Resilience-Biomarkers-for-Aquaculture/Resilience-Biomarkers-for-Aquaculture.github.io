---
layout: post
title: Reproduce Oyster WGBS analysis with nf-core methylseq
author: Shelly Wanamaker
tags: methylseq WGBS params klone 
---
## Background

Steven ran Bismark to align 24 _C. gigas_ samples but ran into issues using `ckpt` on klone because of the time limit. He eventually was able to complete the analysis here: [https://sr320.github.io/tumbling-oysters/posts/40-HawsOA/](https://sr320.github.io/tumbling-oysters/posts/40-HawsOA/) 

He asked if I could get it to run with the nf-core pipeline [`methylseq`](https://nf-co.re/methylseq/2.6.0/)

## Results

After correctly modifying the process section of the `.config` file to match his alignment parameters, I was able to successfully reproduce his results in 16.5 hours:

- output: [20250415_methylseq](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250415_methylseq/nf_report.html)
- mutliqc:[https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250415_methylseq/multiqc/bismark/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250415_methylseq/multiqc/bismark/multiqc_report.html)
- pipeline report: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250415_methylseq/nf_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250415_methylseq/nf_report.html)

## Methods
### Modifying `.config` file to include Bismark parameters

I added the following chunck to the `.config`
file. Full `.config` file is on Klone at `/gscratch/srlab/strigg/bin/uw_hyak_srlab.config` and [here](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250415_methylseq/uw_hyak_srlab.config).

```
  withName: BISMARK_ALIGN {
        ext.args = {
          [
            '--score_min L,0,-0.8',
	    '--non_directional'
          ].join(' ')
        }
        ext.prefix = { "${meta.id}" }
    }
   
   withName: BISMARK_METHYLATIONEXTRACTOR {
        ext.args = {
	  [
	   '--merge_non_CpG',
	   '--comprehensive',
	   '--multicore 48',
	   '--buffer_size 75%'
	  ].join(' ')
	}
	ext.prefix = { "${meta.id}" }
    }

   withName: BISMARK_COVERAGE2CYTOSINE {
        ext.args = {
          [
            '--merge_CpG',
            '--zero_based'
          ].join(' ')
        }
        ext.prefix = { "${meta.id}" }
    }
}
```
### Running the methylseq pipeline on Klone

I ran the following code to reprocess the data:

```
#start screen session 
screen -S methylseq

#start interactive node
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=72:00:00

#activate nextflow conda environment
mamba activate nextflow

#run methylseq pipeline
nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250401_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250415_methylseq \
--fasta /gscratch/srlab/sr320/github/project-oyster-oa/data/Haws-11/GCF_963853765.1_xbMagGiga1.1_genomic.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html \
--skip_trimming \
--nomeseq 

```

### Other attempts that FAILED

1. In this attempt, I was using the default alignment parameters (e.g. score_min L,0,-0.2) so the %aligned was much lower than Steven saw.

	- pipeline report: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250312_methylseq/nf_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250312_methylseq/nf_report.html)
	- output: [20250312_methylseq](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250312_methylseq/)
	- multiqc: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250312_methylseq/multiqc/bismark/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250312_methylseq/multiqc/bismark/multiqc_report.html)

2. In this attempt, I had the incorrect syntax in the `.config` file

	- pipeline report: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250331_methylseq/nf_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250331_methylseq/nf_report.html)
	- output: [20250331_methylseq](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250331_methylseq/)


3. In this attempt, I still had the incorrect syntax in the `.config` file. See slack thread here: [https://nfcore.slack.com/archives/CP3RJSMF0/p1743520726837279](https://nfcore.slack.com/archives/CP3RJSMF0/p1743520726837279)

- pipeline report: [https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250401_methylseq/nf_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250401_methylseq/nf_report.html)
- output: [20250401_methylseq](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_methylseq/20250401_methylseq/)












