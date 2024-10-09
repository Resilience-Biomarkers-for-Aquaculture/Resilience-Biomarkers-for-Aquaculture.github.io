---
layout: post
title: UW bioinformatic resources and project set-up
tags: project-setup HPC-server UW-resources
---

## Project Set-up 

Project will be completed on UW's Klone server  
- Login: `ssh xx@klone.hyak.uw.edu` with username replacing xx, provide password, and confirm DUO push notification  

Roberts lab information on Klone node: [https://robertslab.github.io/resources/klone_Data-Storage-and-System-Organization/]

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

Steven Roberts lab handbook: [https://robertslab.github.io/resources/Data-Management/]  
UW's system: [https://hyak.uw.edu/systems]    
UW's how to run job: [https://hyak.uw.edu/docs/hyak101/basics/jobs]  

