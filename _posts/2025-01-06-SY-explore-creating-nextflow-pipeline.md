---
layout: post
title: Creating a Nextflow pipeline for gene count analysis
author: Steve Yost
tags: NextFlow 
---

# Motivation
Given the interesting results found during the previously posted ChatGPT work,
there was some sentiment that the code I had consolidated from that session
might be reusable if it were generalized. I used this as an opportunity to
explore creating a NextFlow pipeline, taking that python script as a starting
point. 

# Results
A working pipeline is
[here](https://github.com/Resilience-Biomarkers-for-Aquaculture/gmgi-nf-gene-count-analysis).
I've run it successfully on the Seqera platform, producing result files -- heatmaps, a PCA plot, and a TSV file -- in the configured S3 bucket path.

# Method
## Reflecting on ChatGPT
Having had success with ChatGPT interaction, I started this by attempting to have it generate an entire NextFlow pipeline as a zip file. This was not as successful
as before; ChatGPT generated a config file with bad sytax (perhaps simply outdated). Working through these errors was time-consuming, especially starting with my faulty assumption that ChatGPT's code should be correct.
In retrospect
it would have been more efficient if I'd started by thoroughly reading the basic [NextFlow documentation](https://www.nextflow.io/docs/latest/), espcially to get
a better sense of its philosophy and approach. However, muddling through did
cause me to touch upon many aspects of NextFlow that I wouldn't have encountered
by starting from the basics.
## Containers
A key realization was the need to create a Docker container or other virtualization
to run this -- Sequera doesn't handle this. My final script used several packages, and finding versions that didn't have version-conflicting dependencies on lower-level package was a challenge. Seqera's container-building service was
not effective. I eventually found a solution by using `conda` locally to load and automatically resolve
versions, and building my own container and posting it to DockerHub.
## Compute environment
... to be continued.




