---
layout: post
title: Creating a Nextflow pipeline for gene count analysis
author: Steve Yost
tags: NextFlow 
---

# Motivation
Given the interesting results found during the previously posted ChatGPT work,
there was some sentiment that the python script that I had consolidated from that session
might be reusable if it were generalized. I used this as an opportunity to
explore creating a NextFlow pipeline, taking that python script as a starting
point. 

# Results
A working pipeline is
[here](https://github.com/Resilience-Biomarkers-for-Aquaculture/gmgi-nf-gene-count-analysis).
I've run it successfully on the Seqera platform, producing result files -- heatmaps, a PCA plot, and a TSV file -- in the configured S3 bucket path.

# Method
## Reflecting on ChatGPT
Having had success with ChatGPT interaction, I started this effort by attempting to have ChatGPT generate an entire NextFlow pipeline as a zip file, based on my python script. This was not as successful
as before; ChatGPT generated a `nextflow.config` file with bad sytax (perhaps simply outdated). Working through these errors was time-consuming, especially starting with my faulty assumption that ChatGPT's code would be correct.
In retrospect,
it would have been more efficient if I'd started by thoroughly reading the basic [NextFlow documentation](https://www.nextflow.io/docs/latest/), espcially to get
a better sense of its philosophy and approach. However, muddling through did
cause me to touch upon many aspects of NextFlow that I wouldn't have encountered
by starting from the basics.
## Containers
A key realization was the need to create a Docker container or other virtualization
to run this -- Sequera doesn't handle this. My python script used several major packages, and finding versions that didn't have version-conflicting dependencies on lower-level package was a challenge. Seqera's container-building service was
not capable of resolving conflicts, and I produced several failed container builds in sequence. I eventually found a solution to this classic python "versioning hell" by using `conda` locally to load and automatically resolve package
versions, and building my own Docker container and [posting it to DockerHub](https://hub.docker.com/r/journeymansix/gmgi_nf_gene_count_pca).
## Parameter schema
Creating a `nextflow_schema.json` file allowed Seqera to present the
familiar friendly interface for setting parameters, including browsing
S3 for files and directories.

<figure>
    <img src="/assets/seqera-gmgi-nextflow-params.png" alt="Seqera params page"/>
    <figcaption class="caption">Seqera Run params page</figcaption>
</figure>

## Compute environment
The current implementation makes use of my pre-configured AWS Batch compute environment. I haven't yet run it locally, though there's the promise of doing so with alternate configuration parameters.
## Writing to AWS S3
In previous work, packages that were used were able to transparently handle writing to either a local file system or an S3 bucket. This wasn't the case for `matplotlib`'s `saveFigure` function, and so I wrote a function to emulate that behavior.
## Future work
To pursue this further as both a learning experiment and something practically useful, I would do the following:
* Break the single script down into its discrete `process` steps, each of which uses a distinct set
of python packages. Each `process`'s outputs and inputs (between the steps) would be handled by NextFlow `Channels`.
* Specializing each process might allow the use of *existing* Seqera containers for each step, as each process could use its own container with few package dependencies.
* Parallelization of the mutual-information iterations via NextFlow channel multiplicities. This would be an artificial usage -- too cumbersome for this scale -- but a good learning exercise.



