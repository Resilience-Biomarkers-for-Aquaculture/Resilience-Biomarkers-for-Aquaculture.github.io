---
layout: post
title: Run ceasmallr (C. virginica) WGBS data through nf-core methylseq pipeline
author: Shelly Wanamaker
tags: methylseq nf-core WGBS ceasmallr virginica
---




### Data needed:

**NOTE:** the fastq files linked below were incorrect and contained errors. SJW reprocessed the reads to correct these errors on in [November 2024](https://robertslab.github.io/sams-notebook/posts/2024/2024-11-15-Trimming-and-QC---C.virginica-CEASMALLR-FastQs-Using-fastp-BBTools-and-FastQC-MultiQC/index.html).

- Fastqs: 
	- Correct fastqs: [https://gannet.fish.washington.edu/gitrepos/ceasmallr/output/00.00-trimming-fastp/](https://gannet.fish.washington.edu/gitrepos/ceasmallr/output/00.00-trimming-fastp/)
	- INCORRECT fastqs: [https://gannet.fish.washington.edu/seashell/bu-github/ceasmallr/data/](https://gannet.fish.washington.edu/seashell/bu-github/ceasmallr/data/)
- Genome: [https://gannet.fish.washington.edu/seashell/bu-github/ceasmallr/data/genome/Cvirginica_v300.fa](https://gannet.fish.washington.edu/seashell/bu-github/ceasmallr/data/genome/Cvirginica_v300.fa)
- samplesheet: 
	- Correct samplesheet: [Cvirg_methylseq/20250506_methylseq/samplesheet.csv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/20250506_methylseq/samplesheet.csv)
	- INCORRECT samplesheet: [Cvirg_methylseq/samplesheet.csv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/samplesheet.csv)

## Methods

### Step 1. Change paths in samplesheet.csv

```
# make new analysis directory and change into it

### analysis with correct reads ###
mkdir /gscratch/scrubbed/strigg/analyses/20250506_methylseq

# create samplesheet with paths of correct reads

cat /gscratch/scrubbed/strigg/analyses/20250502_methylseq/data/00.00-trimming-fastp/fastq_pairs.txt |\
awk -F"_" '{print $1",""/gscratch/scrubbed/strigg/analyses/20250502_methylseq/data/00.00-trimming-fastp/"$0}' |\ # print sample name without suffix, print a comma then the file path, then the file name
 awk '{print $1}' | \ # only print the sample name and file path
 awk -F"," '{print $1","$2","$2}' | \ # print file path in column 2 and 3
 awk 'BEGIN{FS=OFS=","}{gsub(/R1/, "R2", $3)}1'| \ #sub R2 for R1 in column 3
 awk 'BEGIN{print "sample,fastq_1,fastq_2"}1 \   #add header
 > samplesheet.csv





### analysis with incorrect reads ###
mkdir 20250502_methylseq
cd /gscratch/scrubbed/strigg/analyses/20250502_methylseq

# copy samplesheet to klone

wget https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/samplesheet.csv

# change paths of fastqs in samplesheet
sed 's/\/gscratch\/scrubbed\/sr320\/github\/ceasmallr\/data/\/gscratch\/scrubbed\/strigg\/analyses\/20250502_methylseq\/data/g' samplesheet.csv > ss.csv

#rename samplesheet
mv ss.csv samplesheet.csv
```

### Step 2. Copy input data (fastqs and genome) to analysis directory

```
#get genome file
wget https://gannet.fish.washington.edu/seashell/bu-github/ceasmallr/data/genome/Cvirginica_v300.fa


### correctly trimmed reads ###
#copy reads to klone from gannet
rsync --verbose --progress --archive shellytrigg@gannet.fish.washington.edu:/volume2/web/gitrepos/ceasmallr/output/00.00-trimming-fastp .


### incorrect reads ###
rsync --verbose --progress --archive shellytrigg@gannet.fish.washington.edu:/volume2/web/seashell/bu-github/ceasmallr/data .
```

### Step 3. Run methylseq pipeline

```
# start screen session; I already had one started so rejoined it

screen -r methylseq

# start interactive node

salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=50GB --time=2-00:00:00

# activate nextflow environment

mamba activate nextflow

# run methylseq pipeline from /gscratch/scrubbed/strigg/analyses/20250506_methylseq

nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250506_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250506_methylseq \
--fasta /gscratch/scrubbed/strigg/analyses/20250502_methylseq/data/genome/Cvirginica_v300.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html \
--skip_trimming \
--nomeseq 



### pipeline with incorrect reads ###
# ran in /gscratch/scrubbed/strigg/analyses/20250502_methylseq

nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250502_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250502_methylseq \
--fasta /gscratch/scrubbed/strigg/analyses/20250502_methylseq/data/genome/Cvirginica_v300.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html \
--skip_trimming \
--nomeseq 

```

this same error occured twice so will try pipeline without preseq
```
ERROR ~ Error executing process > 'NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP (CF01-CM01-Zygote)'                                                  
                                                                                                                                                   
Caused by:                                                                                                                                         
  Process `NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP (CF01-CM01-Zygote)` terminated with an error exit status (1)                                 
                                                                                                                                                   
                                                                                                                                                   
Command executed:                                                                                                                                  
                                                                                                                                                   
  preseq \                                                                                                                                         
      lc_extrap \                                                                                                                                  
       \                                                                                                                                           
      -pe \                                                                                                                                        
      -output CF01-CM01-Zygote.lc_extrap.txt \                                                                                                     
      CF01-CM01-Zygote.sorted.bam                                                                                                                  
  cp .command.err CF01-CM01-Zygote.command.log                                                                                                     
                                                                                                                                                   
  cat <<-END_VERSIONS > versions.yml                                                                                                               
  "NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP":                                                                                                    
      preseq: $(echo $(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*$//')                                                                       
  END_VERSIONS                                                                                                                                     
                                                                                                                                                   
Command exit status:                                                                                                                               
  1                                                                                                                                                
                                                                                                                                                   
Command output:                                                                                                                                    
  (empty)                                                                                                                                          
                                                                                                                                                   
Command error:                                                                                                                                     
  INFO:    Environment variable SINGULARITYENV_TMPDIR is set, but APPTAINERENV_TMPDIR is preferred                                                 
  INFO:    Environment variable SINGULARITYENV_NXF_TASK_WORKDIR is set, but APPTAINERENV_NXF_TASK_WORKDIR is preferred                             
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred                                           
  ERROR:        max count before zero is less than min required count (4) duplicates removed                                                       
                                                                                                                                                   
Work dir:                                                                                                                                          
  /mmfs1/gscratch/scrubbed/strigg/analyses/20250506_methylseq/work/e5/64279ad453959ee490b0430a9eb77e                                               
                                                                                                                                                   
Container:                                                                                                                                         
  /gscratch/scrubbed/srlab/.apptainer/depot.galaxyproject.org-singularity-preseq-3.2.0--hdcf5f25_6.img                                             
                                                                                                                                                   
Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`                                  
                                                                                                                                                   
 -- Check '.nextflow.log' file for details                                                                                                         
ERROR ~ Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting           
```



```
nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250506_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250506_methylseq \
--fasta /gscratch/scrubbed/strigg/analyses/20250502_methylseq/data/genome/Cvirginica_v300.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html \
--skip_trimming \
--nomeseq \
--skip_preseq 
```

that didn't work so i added the preseq ignore option to my config file as described here: [https://github.com/nf-core/methylseq/issues/161](https://github.com/nf-core/methylseq/issues/161)

```
process {
    errorStrategy = 'retry'
    maxSubmitAwait = '60 min'
    maxRetries = 2
    executor = 'slurm'
    queue  = { task.attempt < 4 ? (task.attempt < 3 ? 'ckpt' : 'cpu-g2' ) : 'cpu-g2-mem2x' }
    clusterOptions = { task.attempt < 4 ? (task.attempt < 3 ? "-A srlab" : "-A coenv" ) : "-A srlab" }
    scratch = '/gscratch/scrubbed/srlab/'
    resourceLimits = [
        cpus: 16,
        memory: '150.GB',
        time: '72.h'
   ]

   withName:preseq {
        errorStrategy = 'ignore'
    }   
  
   withName: BISMARK_ALIGN {
        ext.args = {
          [


```

Adding that parameter did not stop preseq from running unfortunately. I tried updating nextflow with `nextflow self-update` and running the pipeline again

now I got the error `terminated with an error exit status (141)`

will try running whole pipeline in a new dir tomorrow.



**2025-05-08**


```
# start interactive node

salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=50GB --time=2-00:00:00

# activate nextflow environment

mamba activate nextflow

# run methylseq pipeline from /gscratch/scrubbed/strigg/analyses/20250508_methylseq

nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250508_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250508_methylseq \
--fasta /gscratch/scrubbed/strigg/analyses/20250502_methylseq/data/genome/Cvirginica_v300.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html \
--skip_trimming \
--nomeseq 

```
got another error about job submission. there is still something weird going on because no ckpt jobs are running and hyakalloc shows thousands of idle cpus. 

I attempted to update the methylseq container because I was getting a message that the was a more recent version available.

I ran `nextflow pull https://github.com/nf-core/methylseq`.

I no longer get the note about the container but I am still getting errors every time the pipeline tries to submit to ckpt. 


```
#pipeline trace output for one sample
15	90/f379f5	NFCORE_METHYLSEQ:METHYLSEQ:FASTQC (CF08-CM04-Larvae)	FAILED	1	-	-	-	-	-	ckpt-all	-	-
48	d8/ee9d82	NFCORE_METHYLSEQ:METHYLSEQ:FASTQC (CF08-CM04-Larvae)	FAILED	1	-	-	-	-	-	ckpt-all	-	-


90/f379f5 .command.run sbatch script header
#!/bin/bash
#SBATCH -J nf-NFCORE_METHYLSEQ_METHYLSEQ_FASTQC_(CF08-CM04-Larvae)
#SBATCH -o /mmfs1/gscratch/scrubbed/strigg/analyses/20250508_methylseq/work/90/f379f5ddf821f9fba8b08d4a8baf91/.command.log
#SBATCH --no-requeue
#SBATCH --signal B:USR2@30
#SBATCH -c 6
#SBATCH -t 08:00:00
#SBATCH --mem 36864M
#SBATCH -p ckpt-all
#SBATCH -A srlab

[d8/ee9d82] .command.run sbatch script header
#!/bin/bash
#SBATCH -J nf-NFCORE_METHYLSEQ_METHYLSEQ_FASTQC_(CF08-CM04-Larvae)
#SBATCH -o /mmfs1/gscratch/scrubbed/strigg/analyses/20250508_methylseq/work/d8/ee9d8285ef39131054efec88d2535e/.command.log
#SBATCH --no-requeue
#SBATCH --signal B:USR2@30
#SBATCH -c 12
#SBATCH -t 16:00:00
#SBATCH --mem 73728M
#SBATCH -p ckpt-all
#SBATCH -A srlab



[d8/ee9d82] NOTE: Error submitting process 'NFCORE_METHYLSEQ:METHYLSEQ:FASTQC (CF08-CM04-Larvae)' for execution -- Execution is retried (
2)

```
I figured out the problem by going into the work directory with where the error occured and tried running the sbatch script by running `sbatch .command.run`. This gave the following error:

```
sbatch: error: When submitting a Checkpoint job, you cannot specify --no-requeue
sbatch: error: Batch job submission failed: Unspecified error
```

You can see in the headers of the sbatch scripts they have a `--norequeue` parameter. Interestingly this was not a problem on 2025-04-22 when I ran the methylseq pipeline and it completely all jobs submitted to ckpt. There must have been an update to the ckpt partition that is now causing this to happen. 


**2025-05-09**

Modified the .config file to get around the `--norequeue` default as suggested by Carson Miller:

A quick "hack" to workaround this for the time being would be to do include the following in your `nextflow.config` .
`process {
    clusterOptions = { "--requeue" }
    errorStrategy = { task.exitStatus in ((130..145)) ? 'ignore' : 'retry' }
}`
Basically this will include `--requeue` in the SLURM submission so that you can still use checkpoint queue. The caveat is that if a job is pre-empted, Nextflow will think that the job failed and will try to create a new job, but the SLURM scheduler will also resubmit the job (even though this won't be tracked by Nextflow), so what you'll end up with is duplication of jobs and a potentially very large overhead. So I added the errorStrategy line so that Nextflow will ignore pre-emptions, with the idea that I will just resume the pipeline to clean up tasks that were pre-empted. 


Edited .config file:

```
params {
    config_profile_description = 'UW Hyak Roberts labs cluster profile provided by nf-core/configs.'
    config_profile_contact = 'Shelly A. Wanamaker @shellywanamaker'
    config_profile_url = 'https://faculty.washington.edu/sr320/'
}
 
process {
    errorStrategy = {task.exitStatus in ((130..145)) ? 'ignore' : 'retry'}
    maxSubmitAwait = '60 min'
    maxRetries = 2
    executor = 'slurm'
    queue  = { task.attempt < 4 ? (task.attempt < 3 ? 'ckpt-all' : 'cpu-g2' ) : 'cpu-g2-mem2x' }    
    clusterOptions = { task.attempt < 4 ? (task.attempt < 3 ? "-A srlab --requeue" : "-A coenv --requeue" ) : "-A srlab --requeue" } 
    scratch = '/gscratch/scrubbed/srlab/'
    resourceLimits = [
        cpus: 16,
        memory: '150.GB',
        time: '72.h'
   ]
 
   withName:preseq { 
        errorStrategy = 'ignore' 
    }  
  
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


 
executor {
    queuesize = 50
    submitRateLimit = '1 sec'
}
 
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/gscratch/scrubbed/srlab/.apptainer'
}

trace {
    enabled = true
    file = 'pipeline_trace.txt'
    fields = 'task_id,hash,name,status,exit,submit,duration,realtime,%cpu,rss,vmem,rchar,wchar,disk,queue,hostname,cpu_model'
}
 
debug {
    cleanup = false
}
```

Re-ran the pipeline in a new directory

```
# start interactive node

salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=50GB --time=2-00:00:00

# activate nextflow environment

mamba activate nextflow

# run methylseq pipeline from /gscratch/scrubbed/strigg/analyses/20250509_methylseq

nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250509_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250509_methylseq \
--fasta /gscratch/scrubbed/strigg/analyses/20250502_methylseq/data/genome/Cvirginica_v300.fa \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html \
--skip_trimming \
--nomeseq 
```

## Results

Pipeline initally failed because the fastq files had errors as noted in Sam's [November 2024 notebook entry](https://robertslab.github.io/sams-notebook/posts/2024/2024-11-15-Trimming-and-QC---C.virginica-CEASMALLR-FastQs-Using-fastp-BBTools-and-FastQC-MultiQC/index.html).

```
Command error:
  INFO:    Environment variable SINGULARITYENV_TMPDIR is set, but APPTAINERENV_TMPDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_TASK_WORKDIR is set, but APPTAINERENV_NXF_TASK_WORKDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred
  Failed to process file EF07-EM01-Zygote_2.gz
  uk.ac.babraham.FastQC.Sequence.SequenceFormatException: Ran out of data in the middle of a fastq entry.  Your file is probably truncated
  	at uk.ac.babraham.FastQC.Sequence.FastQFile.readNext(FastQFile.java:187)
  	at uk.ac.babraham.FastQC.Sequence.FastQFile.next(FastQFile.java:129)
  	at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:77)
  	at java.base/java.lang.Thread.run(Thread.java:833)
```

- See [nf_report.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_methylseq/20250502_methylseq/pipeline_info/nf_report.html)