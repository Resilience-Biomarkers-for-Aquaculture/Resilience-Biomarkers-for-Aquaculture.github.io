---
layout: post
title: Running RNAseq pipeline with reference transcript for C.Gigas Dataset #1
tags: RNAseq nextflow Klone referencetranscript Cgigas
---

## RNASeq workflow for C.Gigas Dataset #1

This workflow is fully documented in a separate post [here](https://resilience-biomarkers-for-aquaculture.github.io/ES-RNAseq-pipeline-construction-post/).

Datasets:    
- [1: Arredondo-Espinoza et al 2023 C. gigas thermotolerance](#dataset-1-arredondo-espinoza-et-al-2023)    

Data is downloaded via nextflow pipeline [fetchNGS](https://nf-co.re/fetchngs/1.12.0/) and see Shelly's example download post [here](https://resilience-biomarkers-for-aquaculture.github.io/a-fetchNGSKlone/). 

Creating a mamba environment on my user for following workflows:

```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

bash Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge

source ~/mambaforge/bin/activate
mamba init 
## close and reopen shell for changes to take place 

mamba create -n nextflow -c bioconda -c conda-forge python=3.8 nf-core nextflow -y -q
```

Downloading [C.gigas genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1/). Using xbMagGiga1.1 published in December 2023. Ref ID: GCF_963853765.1

```
cd /mmfs1/gscratch/scrubbed/elstrand/genomes
mkdir C_gigas

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1/GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1/GCF_963853765.1_xbMagGiga1.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1/GCF_963853765.1_xbMagGiga1.1_genomic.gtf.gz
```

## Dataset 1: Arredondo-Espinoza et al 2023

Paper found [here](https://www.sciencedirect.com/science/article/pii/S1744117X23000345?via%3Dihub). Shelly downloaded this data on [20240925 with fetchNGS](https://resilience-biomarkers-for-aquaculture.github.io/a-fetchNGSKlone/). 

Supplemental data that might be helpful:  
- [Supplementary Table II](https://ars.els-cdn.com/content/image/1-s2.0-S1744117X23000345-mmc6.xls). Transcripts annotated for the heat-resistant phenotype of C. gigas whose high expression might be useful as a marker of resistance to high temperatures.     
- [Supplementary Table III](https://ars.els-cdn.com/content/image/1-s2.0-S1744117X23000345-mmc7.xls). Transcripts annotated for the heat-susceptible phenotype of C. gigas whose high expression might be useful as a marker of resistance to high temperatures. No sure what the difference is between II and III?  

This study focused on the Pacific oyster, *Crassostrea gigas*, and completed a laboratory-simulated daily oscillatory thermal challenge (26 to 34 °C) for 30 days. 

Path to fastq files: `/mmfs1/gscratch/scrubbed/<ShellynetID>/analyses/20240925/fastq`   
Path to samplesheet created by fetchNGS: `/mmfs1/gscratch/scrubbed/<ShellynetID>/analyses/20240925/samplesheet/samplesheet.csv` 

Path to Emma's analysis: `/mmfs1/gscratch/scrubbed/<EmmanetID>/Cgigas_ArredondoEspinoza2023/`
Path to C.gigas genome: `/mmfs1/gscratch/scrubbed/elstrand/genomes/GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz`

Testing this pipeline: 

```
# create screen session
screen -S rnaseq_cgigas_AE2023

# request a compute node (mem and time requests can be modified)
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=40GB --time=12:00:00

# activate mamba environment (created at beginning of this post)
mamba activate nextflow

# update rnaseq workflow in nextflow within mamba environment
nextflow pull nf-core/rnaseq

## attach back to screen 
screen -r 1266631.rnaseq_cgigas_AE2023

## detach from screen 
Press Ctrl+A, then D (the default key binding)

## export contents of screen session

```

Nextflow RNAsew workflow (`01-rnaseq.sh`):

Start conda environment before running sbatch. 

```
#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --account=srlab
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --time=1-20:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=2
#SBATCH --chdir=/mmfs1/gscratch/scrubbed/elstrand/Cgigas_ArredondoEspinoza2023/scripts

#conda activate nextflow

samplesheet="/mmfs1/gscratch/scrubbed/elstrand/Cgigas_ArredondoEspinoza2023/samplesheet_rnaseq_dataset1.csv"
output="/mmfs1/gscratch/scrubbed/elstrand/Cgigas_ArredondoEspinoza2023/"
genome="/mmfs1/gscratch/scrubbed/elstrand/genomes/C_gigas/GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz"
gff="/mmfs1/gscratch/scrubbed/elstrand/genomes/C_gigas/GCF_963853765.1_xbMagGiga1.1_genomic.gff.gz"
gtf="/mmfs1/gscratch/scrubbed/elstrand/genomes/C_gigas/GCF_963853765.1_xbMagGiga1.1_genomic.gtf.gz"

nextflow run nf-core/rnaseq \
    -resume \
    -profile singularity \
    --input ${samplesheet} \
    --outdir ${output} \
    --gtf ${gtf} \
    --gff ${gff} \
    --fasta ${genome} \
    --trimmer fastp \
    --extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10' \
    --aligner star_salmon \
    --skip_pseudo_alignment \
    --multiqc_title rnaseq_Cgigas_ArredondoEspinoza2023 \
    --deseq2_vst
```

`sbatch 01-rnaseq.sh` to run this script and `squeue -u elstrand` to check the status 

To add pseudo alignment: `--pseudo_aligner salmon`

#### 2024-10-09

Trying the above without any extra flags to see if this runs properly. Will need to adjust trimming based on multiqc output.  

Error below. The samplesheet doesn't have a strandedness column.. This samplesheet.csv was generated by fetchNGS. Is there way to default to auto?

```
ERROR ~ ERROR: Validation of 'input' file failed!

 -- Check '.nextflow.log' file for details
The following errors have been detected:

* -- Entry 1: Missing required value: strandedness
* -- Entry 2: Missing required value: strandedness
* -- Entry 3: Missing required value: strandedness
* -- Entry 4: Missing required value: strandedness
```

#### 2024-10-22 

**Samplesheet issue** 

Coming back to this samplesheet issue.. I'm exporting the .csv file from the fetchNGS output: 

```
## run this outside of klone 
scp xxx@klone.hyak.uw.edu:/mmfs1/gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet.csv C:\Users\EmmaStrand\Downloads

## moved this file to my C:\Users\EmmaStrand\MyProjects\Resilience_Biomarkers_Aquaculture folder that houses the Resilience-Biomarkers-for-Aquaculture.github.io.
## do we want this on our github? folder for files and scripts?

## this below didn't work because I can only read shelly's folder 
scp C:\Users\EmmaStrand\MyProjects\Resilience_Biomarkers_Aquaculture\samplesheet_rnaseq_dataset1.csv elstrand@klone.hyak.uw.edu:/mmfs1/gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet_rnaseq_dataset1.csv

## this worked!
scp C:\Users\EmmaStrand\MyProjects\Resilience_Biomarkers_Aquaculture\samplesheet_rnaseq_dataset1.csv elstrand@klone.hyak.uw.edu:/mmfs1/gscratch/scrubbed/elstrand/Cgigas_ArredondoEspinoza2023/
```

- Created an R script to convert this samplesheet to an rnaseq input samplesheet. We need metadata information for the SRA ID to match to sample ID from the treatments in the expertiment. How do we set this up where it will be a template script for all datasets?   
- We just need an additional column called "strandedness" that is set to "auto" to auto-detect the forward or reverse strand. Exporting file from R script to klone (code in chunk above).   

**Attaching to screen again** 

I'm not able to get back into this screen.. `screen attach-session -t 1266631.rnaseq_cgigas_AE2023`. If I use `screen -ls` I can see it activate but attach session isn't working. But this did `screen -r 1266631.rnaseq_cgigas_AE2023`. Replaced instructions at beginning of post. 

**Running rnaseq** 

When running rnaseq nextflow I got this warning: `NOTE: Your local project version looks outdated - a different revision is available in the remote repository [1f3f64dac7]`. Does this refer to the rnaseq pipeline?

Got error: `ERROR ~ Empty header columns are not allowed in CSV file`. Added `row.names = FALSE` to the export csv function in the R script. Re-uploaded to klone. Try again! 

Got error: Input validation failed. This can't find the fastq files because now that folder only has one sample not all of them.. But Shelly's last post has a screenshot of the whole set. Probably because of 30 day limit on files. Shelly redownloaded to this same folder. 

```
ERROR ~ ERROR: Validation of 'input' file failed!

 -- Check '.nextflow.log' file for details
The following errors have been detected:

* -- Entry 2 - fastq_2: the file or directory '/gscratch/scrubbed/strigg/analyses/20240925/fastq/SRX5644305_SRR8856684_2.fastq.gz' does not exist.
```

Got error: It's running! But needs 6 CPUs not the available 1. Make sbatch script instead of srun and screen. also get output this way.

```
WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path: /mmfs1/home/elstrand/work/singularity -- Use the environment variable NXF_SINGULARITY_CACHEDIR to specify a different location
ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_FASTP:FASTP (SRX5644304)'

Caused by:
  Process requirement exceeds available CPUs -- req: 6; avail: 1
```

For future datasets, can there be a base script that you just give a different config file to every time? So it's not the same sbatch script edited for every different file?

Got error: Need to run mamba init to activate/deactivate mamba or just use conda. Changed sbatch script to just use conda. 

Got error: Fastp error once running, I'm still getting this error.. Trying to add ntasks and cpus per tasks to the sbatch script.

```
ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_FASTP:FASTP (SRX5644305)'

Caused by:
  Process requirement exceeds available CPUs -- req: 6; avail: 1
```

Got error: This is better! but need 12 so upping the ntasks. 

```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ERROR ~ Error executing process > 'NFCORE_RNASEQ:PREPARE_GENOME:STAR_GENOMEGENERATE (GCF_963853765.1_xbMagGiga1.1_genomic.fna)'

Caused by:
  Process requirement exceeds available CPUs -- req: 12; avail: 10
```

Got error: Ugh, getting there. Up the memory allocated.

```
ERROR ~ Error executing process > 'NFCORE_RNASEQ:PREPARE_GENOME:STAR_GENOMEGENERATE (GCF_963853765.1_xbMagGiga1.1_genomic.fna)'

Caused by:
  Process requirement exceeds available memory -- req: 72 GB; avail: 70 GB
```

#### 2024-10-23 

This worked! ~13 hours for this dataset. Yay! Back up output on gannet.

```
rsync --archive --verbose --progress *.html emma.strand@gannet.fish.washington.edu:/volume2/web/emma.strand/rnaseq/Cgigas_ArredondoEspinoza2023
```

#### 2024-10-29 

Deciding on parameters to include for fastp, star, and salmon. 

```
    --extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10' \
```

Options for STAR that I wasn't sure if we needed? (`--extra_star_align_args ''`):    
- `--outFilterMultimapNmax`: Maximum number of multiple alignments allowed   
- `--outFilterMismatchNmax`: Maximum number of mismatches allowed   
- `--alignIntronMin`: Minimum intron size   
- `--alignIntronMax`: Maximum intron size   
- `--alignMatesGapMax`: Maximum gap between paired reads   

Options for SALMON that I wasn't sure we needed (`--extra_salmon_quant_args ''`):    
- `--numBootstraps`: Number of bootstrap samples to generate  
- `--numGibbsSamples`: Number of Gibbs sampling rounds   
- `--fldMean`: Mean fragment length for single-end reads  
- `--fldSD`: Standard deviation of fragment length for single-end reads     
- `--validateMappings`: Perform selective alignment    
- `--gcBias`: Enable GC bias correction  
- `--seqBias`: Enable sequence-specific bias correction  

We also decided we don't need the pseudo-aligning step so I'm adding a flag to skip that (`--skip_pseudo_alignment`). 

Re-running script ~11 am. This took ~12 hours last time with the pseudo alignment. The traditional alignment is the time intensive step so I'm guessing 9-10 hours for this run? Running while nextflow conda environment is activated b/c got an error to run conda init before activate within the sbatch script.

#### 2024-10-30 

Checking on this output from yesterday. Got the error below - I didn't get this issue before. Google says "The error you're encountering suggests that the singularity image pull operation is timing out. This is likely due to network issues or the size of the image". 

```
ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:BAM_MARKDUPLICATES_PICARD:PICARD_MARKDUPLICATES (13)'

Caused by:
  Failed to pull singularity image
    command: singularity pull  --name depot.galaxyproject.org-singularity-picard-3.1.1--hdfd78af_0.img.pulling.1730233288426 https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0 > /dev/null
    status : 143
    hint   : Try and increase singularity.pullTimeout in the config (current is "20m")
    message:
      INFO:    Downloading network image

 -- Check '.nextflow.log' file for details
ERROR ~ Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting

 -- Check '.nextflow.log' file for details
 ```

This also resulted in the final execution log being: 

```
294     9b/fd4a4a       54236   NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:SALMON_QUANT (SRX5644336)     ABORTED -       2024-10-29 13:50:40.680 -       -       -       -       -       -       -
296     db/8eb422       55321   NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:SALMON_QUANT (SRX5644337)     ABORTED -       2024-10-29 13:50:44.275 - 
```

I removed all the previous iterations of running this dataset to start fresh. Activated mamba nextflow outside of sbatch script and ran again ~11 am. This worked, yay!

### 2024-10-31 

Rsync the output to gannet. 

```
rsync --archive --verbose --progress *.html emma.strand@gannet.fish.washington.edu:/volume2/web/emma.strand/rnaseq/Cgigas_ArredondoEspinoza2023
```

### Comparing output to Roberto's paper 

Rysnc output to gannet 

```
rsync --archive --verbose --progress deseq2_qc/*.RData emma.strand@gannet.fish.washington.edu:/volume2/web/emma.strand/rnaseq/Cgigas_ArredondoEspinoza2023
rsync --archive --verbose --progress deseq2_qc/*.txt emma.strand@gannet.fish.washington.edu:/volume2/web/emma.strand/rnaseq/Cgigas_ArredondoEspinoza2023
rsync --archive --verbose --progress salmon.* emma.strand@gannet.fish.washington.edu:/volume2/web/emma.strand/rnaseq/Cgigas_ArredondoEspinoza2023
rsync --archive --verbose --progress tx2gene.tsv emma.strand@gannet.fish.washington.edu:/volume2/web/emma.strand/rnaseq/Cgigas_ArredondoEspinoza2023
```
