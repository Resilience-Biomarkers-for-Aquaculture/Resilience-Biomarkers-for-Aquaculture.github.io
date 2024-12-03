---
layout: post
title: Repoint nf-core projectDir to srlab directory
tags: nextflow nf-core klone
---

## Purpose
The purpose of this post is to describe how to redirect the projectDir when running an nf-core pipeline on Klone. This came up when I attempted to run the methylseq pipeline for the first time and got a disk quota exceeded error (https://github.com/RobertsLab/resources/issues/2045). This happened because nextflow was trying to pull images into my home directory on klone `/mmfs1/home/strigg/.nextflow/assets/nf-core/` but it couldn't because there is a 10GB limit and I had a `work` directory in my home directory that was taking up a lot of space. Once I removed that directory, I was able to run the pipeline without getting the disk quota exceeded error, but it prompted the question "Is this going to continue happening each time I run a new pipeline?" and "Can I repoint nextflow to my /gscratch/srlab/strigg directory which has a lot more space?"

This turned into a discussion on the [nf-core slack under the #configs channel](https://nfcore.slack.com/archives/CGRTMASKY/p1733234847813229?thread_ts=1649954706.174089&cid=CGRTMASKY). With that guidance and following [this guide for setting environmental variables in a conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#setting-environment-variables) (since I run nextflow through mamba), I was able to change the variables `NXF_HOME` and `NXF_TEMP` in order to get nextflow to stop downloading files on my home directory and instead download them on gscratch/srlab/strigg/bin.

## Procedure

1. determine the directory that nextflow is in

  ```
  which nextflow
  /mmfs1/gscratch/srlab/strigg/bin/mambaforge/envs/nextflow/bin/nextflow
  ```

2. move into that directory and confirm the directories `activate.d` and `deactivate.d` exist in the `etc/conda` directory. If they don't, create them.

  ```
  cd /mmfs1/gscratch/srlab/strigg/bin/mambaforge/envs/nextflow/etc/conda/
  ls
  activate.d  deactivate.d
  cd activate.d
  ```

3. create a file called `env_vars.sh` in the `activate.d` directory to set the `NXF_HOME` and the `NXF_TEMP` path (I used vim)

  ```
    #!/bin/sh

    export NXF_HOME='/mmfs1/gscratch/srlab/strigg/bin'
    export NXF_TEMP='/mmfs1/gscratch/scrubbed/strigg/bin'
  ```

4. create a file called `env_vars.sh` in the `deactivate.d` directory to clear the path that contains the following (I used vim)

  ```
  cd ../deactivate.d
  ```

  ```
    #!/bin/sh

    unset NXF_HOME
    unset NXF_TEMP
  ```

5. check variable is set as specified in `env_vars.sh` script in `activate.d` directory

  ```
  mamba activate nextflow
  echo $NXF_HOME
  /mmfs1/gscratch/srlab/strigg/bin
  echo $NXF_TEMP
  /mmfs1/gscratch/scrubbed/strigg/bin
  ```

6. check variable is unset after environment has been deactivated

  ```
  mamba deactivate
  echo $NXF_HOME
  ```

7. test NXF_HOME is set

  ```
  # request a compute node (mem and time requests can be modified)
  salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=12:00:00

  # load the nextflow environment
  mamba activate nextflow

  # run nextflow pipeline
  nextflow run \
  nf-core/fetchngs \
  -resume \
  -c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
  --input /gscratch/scrubbed/strigg/analyses/20241203_FetchNGStest/ids.csv \
  --outdir /gscratch/scrubbed/strigg/analyses/20241203_FetchNGStest \
  --download_method sratools
  ```

![Screenshot 2024-12-03 122730](https://github.com/user-attachments/assets/1010bf7f-d731-40bf-9e8b-85754c1cbf0a)


## Output
https://gannet.fish.washington.edu/metacarcinus/Nf-core_tests/20241203_FetchNGStest/
[.nextflow.log file](https://gannet.fish.washington.edu/metacarcinus/Nf-core_tests/20241203_FetchNGStest/.nextflow.log)

**Other iterations that Failed**

I initially only had NXF_HOME defined in the env_vars.sh file as `export NXF_HOME='$NXF_HOME:/mmfs1/gscratch/srlab/strigg/bin'`. This gave the following error when I ran the fetchngs pipeline test (excerpt from .nextflow.log):

```
Dec-03 08:08:06.869 [main] DEBUG nextflow.plugin.PluginsFacade - Plugins declared=[nf-validation@1.1.3]
Dec-03 08:08:06.870 [main] DEBUG nextflow.plugin.PluginsFacade - Plugins default=[]
Dec-03 08:08:06.870 [main] DEBUG nextflow.plugin.PluginsFacade - Plugins resolved requirement=[nf-validation@1.1
.3]
Dec-03 08:08:06.870 [main] DEBUG nextflow.plugin.PluginUpdater - Installing plugin nf-validation version: 1.1.3
Dec-03 08:08:06.872 [main] INFO  nextflow.plugin.PluginUpdater - Downloading plugin nf-validation@1.1.3
Dec-03 08:08:07.945 [main] INFO  org.pf4j.util.FileUtils - Expanded plugin zip 'nf-validation-1.1.3.zip' in 'nf-
validation-1.1.3'
Dec-03 08:08:07.947 [main] DEBUG nextflow.plugin.PluginUpdater - Failed atomic move for plugin /tmp/pf4j-update-
downloader9126177795036055700/nf-validation-1.1.3 -> $NXF_HOME:/mmfs1/gscratch/srlab/strigg/bin/plugins/nf-valid
ation-1.1.3 - Reason: /tmp/pf4j-update-downloader9126177795036055700/nf-validation-1.1.3 -> $NXF_HOME:/mmfs1/gsc
ratch/srlab/strigg/bin/plugins/nf-validation-1.1.3: Invalid cross-device link - Fallback on safe move
Dec-03 08:08:08.022 [main] ERROR nextflow.cli.Launcher - @unknown
org.pf4j.PluginRuntimeException: Cannot find the manifest path
        at org.pf4j.ManifestPluginDescriptorFinder.readManifestFromDirectory(ManifestPluginDescriptorFinder.java
:142)
        at org.pf4j.ManifestPluginDescriptorFinder.readManifest(ManifestPluginDescriptorFinder.java:72)
        at org.pf4j.ManifestPluginDescriptorFinder.find(ManifestPluginDescriptorFinder.java:58)
```

Googling the error `Invalid cross-device link` gave the explanation that the destination directory is not on the same partition as the source directory. In this case, the /tmp/ directory and the /mmfs1/ directories are on different partitions. Why nextflow downloads to a tmp directory and then moves the files is a mystery to me. But this prompted me to identify another environment variable that needs to be redefined in addition to `NXF_HOME` and that's `NXF_TEMP`.

In my second attempt, my env_vars.sh file contained `export NXF_HOME='$NXF_HOME:/mmfs1/gscratch/srlab/strigg/bin'` and `export NXF_TEMP='$NXF_TEMP:/mmfs1/gscratch/scrubbed/strigg/bin'`. But this led to the error:

```
Dec-03 09:18:35.365 [main] DEBUG nextflow.plugin.PluginUpdater - Installing plugin nf-validation version: 1.1.3
Dec-03 09:18:35.367 [main] ERROR nextflow.cli.Launcher - $NXF_TEMP:/gscratch/scrubbed/strigg/bin:/gscratch/scrubbed/strigg/bin/nextflow-plugin-nf-validation-1.1.3.lock (No such file or directory)
java.io.FileNotFoundException: $NXF_TEMP:/gscratch/scrubbed/strigg/bin:/gscratch/scrubbed/strigg/bin/nextflow-plugin-nf-validation-1.1.3.lock (No such file or directory)

```

I realized I needed to remove the `$NXF_TEMP:` and `$NXF_HOME:` from the variable definition in the env_vars.sh file. Once I did that, the pipeline completed without errors.
