---
layout: post
title: Set up AWS Batch pipeline on Seqera and run rnaseq test
tags: nextflow seqera aws rnaseq test
---

## 11-23-2024

The goal here was initial familiarization with Seqera, particularly configuring my own pipeline using the AWS Batch (one of several cloud choices available) [compute environment](https://docs.seqera.io/platform/23.2/compute-envs/overview) and running `rnaseq` using the `test` config. Success would amount to replicating the same results produced by running `rnaseq` in Seqera's community showcase, which provides a free pre-configured compute environment that also uses AWS Batch.

## Sequence
### Community Showcase
I first explored the Community Showcase. I noted there's an `rnaseq` pipeline defined there: [https://cloud.seqera.io/orgs/community/workspaces/showcase/launchpad/189046775836482](https://cloud.seqera.io/orgs/community/workspaces/showcase/launchpad/189046775836482). It uses a compute environment named `AWS_Batch_Ireland_FusionV2_NVMe`. The training video demonstrates using the `rnaseq` pipeline. I ran that pipeline successfully, following the instructions.
Here's the commandline it produced:
```
nextflow run https://github.com/nf-core/rnaseq
		 -name zen_yost
		 -params-file https://api.cloud.seqera.io/ephemeral/SGqqOCgLhH4wH_cIWvZy3w.json
		 -with-tower https://api.tower.nf
		 -r 3.14.0
		 -profile test
```

Seqera provides 100 free CPU hours. Running this pipeline used only 0.7 CPU hours, at an estimated cost of $0.024.

> One note: I found that when running a pipeline, the status page didn't update dynamically as expected. I had to refresh the browser page to get updated status.

### Creating the AWS Batch compute environment
I then wished to replicate that result with my own compute environment, which is required for any real Seqera work beyond familiarization.
[Sequera documentation](https://docs.seqera.io/platform/24.2/compute-envs/aws-batch) for setting up an AWS Batch compute environment was thorough and, with a few minor exceptions, up to date with the AWS console web interface.

My first attempt to create a Batch Forge compute environment, following the above instructions for setting up AWS resources, failed with the message "Unable to find default AWS subnets for vpc none and region us-east-1".  A Seqera AI doc query suggested (among many other things) that I needed a default VPC on AWS. It turned out I indeed didn't have a default VPC, so I created one here: [https://us-east-1.console.aws.amazon.com/vpcconsole/home?region=us-east-1#vpcs](https://us-east-1.console.aws.amazon.com/vpcconsole/home?region=us-east-1#vpcs])
`Actions -> Create default VPC`.

I then repeated the same steps to create a new Batch Forge compute environment (because it didn't allow editing the failed one), with success this time. I deleted the failed environment.

I then created an `rnaseq` pipeline using the AWS Batch-based compute environment, and the same test profile. I lannched and ran that successfully, producing output in my AWS S3 bucket. It used 0.9 CPU hours, with an estimated cost of $0.024.

The command line it produced was as follows:
```
nextflow run 'https://github.com/nf-core/rnaseq'
-name nostalgic_crick
-params-file 'https://api.cloud.seqera.io/ephe...'
-with-tower
-r 3.14.0
-profile test
```
... which is nearly identical to the one created by the commununity showcase pipeline above.


> Some naive dataset inferences: I see that in the community showcase Datasets, there's a dataset named `nf-core-rnaseq-test`, which has columns `sample, fastq_1, fastq_2, strandedness`. Not all datasheets there with "rnaseq" in their names use this schema, though some do. `rnaseq_1` has `fastq_1` and `fastq_2`, but no strandedness, and has others. Note that most of these are user provided datasets.

### Next steps
Next I'll attempt to replicate the work in [https://resilience-biomarkers-for-aquaculture.github.io/ES-RNAseq_with_reference_dataset1/](https://resilience-biomarkers-for-aquaculture.github.io/ES-RNAseq_with_reference_dataset1/) on Seqera. I'll make a separate post for that effort.


