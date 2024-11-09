---
layout: post
title: UW bioinformatic resources and project set-up
tags: project-setup HPC-server UW-resources
---

## Project Set-up

Project will be completed on UW's Klone server  
- Login: `ssh xx@klone.hyak.uw.edu` with username replacing xx, provide password, and confirm DUO push notification  

Roberts lab information on Klone node: https://robertslab.github.io/resources/klone_Data-Storage-and-System-Organization/

Server structure:

```
[xx@klone-login03 ~]$ cd ../../
[xx@klone-login03 mmfs1]$ ls
admin  apsearch  data  encrypted  gscratch  home  slurmdata  ssg  sw

## path to our working space on srlab: /mmfs1/gscratch/srlab/<shelly-UW_NetID>
```

Storage:  
- Login node (`/mmfs1/home/<UW_NetID>`): 10 GB    
- Roberts lab has 1 TB in `srlab` folder (`/gscratch/srlab/`)  
- Temporary storage ("Scrubbed"): 200 TB in user folder (`/gscratch/scrubbed/<UW_NetID>`)  

To see space and file utilization: `hyakstorage`

## Data management resources

Steven Roberts lab handbook: https://robertslab.github.io/resources/Data-Management/    
UW's system: https://hyak.uw.edu/systems     
UW's how to run job: https://hyak.uw.edu/docs/hyak101/basics/jobs   

## Tips and tricks 

`squeue -A srlab`: check what jobs are being run on srlab
`squeue -A srlab -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %c %m"`: check jobs and show CPUs and memory used by each job     
`hyakalloc -g srlab`: check how many resources are in use and free for use     
`hyakalloc -p ckgt`: check what resources are currently available on all of hyak    

Shelly added config file for nextflow `/gscratch/srlab/strigg/bin/uw_hyak_srlab.config` that will use other nodes if nobody is using them.

## Hyak info

Total GB:    
Total cpu: 32    
Partition to use: `-p cpu-g2-mem2x`

We don't have access to ckpt at the moment so config file is not using other resources
