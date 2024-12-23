---
layout: post
title: PCA for gene count comparison between 2012 and 2023 reference genomes
tags: rnaseq pca 
---

## 12-23-2024

The goal here was to compare gene counts between `rnaseq` runs using
1. [The NCBI reference genome and gene annotations from the period of Roberto's study](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000297895.1/) (herein referred to as `2012`)
2. The most recent [reference genome and gene annotations](https://ncbi.nlm.nih.gov/datasets/genome/GCA_963853765.1/) (herein referred to as `2023`) for C. gigas.

## Motivation
As previously [noted](https://resilience-biomarkers-for-aquaculture.github.io/SW-CGI_ID_matching_attemp_1/), in attemping to compare gene counts between
Roberto's results and `rnaseq` results, we encountered the following hurdle:
The gene count matrix that Roberto provided overwhelmingly used gene IDs prefixed with `MSTRG_` (novel or unannotated genes/transcripts identified by
StringTie) and `CGI_` (CpG islands), whereas the avaialbe NCBI annotations, and thus gene count data, for `2023` uses gene IDs prefixed with `LOC_`,
commonly used for annotated genes or predicted genes that do not yet have an official gene symbol or name.
However, we found that the `2012` NCBI genome and annotations did use `LOC_` gene IDs, which facilitates matching by ID
rather than requiring BLAST or other indirect gene ID matching methods.

So, for this current task, we used the `2012`
reference genome (linked above) as as a proxy for Roberto's, given that it represents the state of the art at the time of his
publication.

## Method
We ran `rnaseq` on Seqera usnig `2012` using the following parameters, which are identical to the `2023` run with the obvious exceptions of reference genome and annotation paths:
```
{
    "remove_ribo_rna": false,
    "help_full": false,
    "custom_config_base": "https://raw.githubusercontent.com/nf-core/configs/master",
    "skip_deseq2_qc": false,
    "gencode": false,
    "umitools_dedup_stats": false,
    "plaintext_email": false,
    "save_reference": false,
    "skip_markduplicates": false,
    "ribo_database_manifest": "/.nextflow/assets/nf-core/rnaseq/workflows/rnaseq/assets/rrna-db-defaults.txt",
    "monochrome_logs": false,
    "aligner": "star_salmon",
    "featurecounts_group_type": "gene_biotype",
    "save_bbsplit_reads": false,
    "skip_multiqc": false,
    "skip_preseq": true,
    "skip_dupradar": false,
    "save_align_intermeds": false,
    "gtf": "https://steveyost-seqera.s3.us-east-1.amazonaws.com/Cgigas_ArredondoEspinoza2023/run_ncbi_dataset_GCA_000297895.1_GCF/genomic_removed_empty_gene_ids.gtf",
    "max_multiqc_email_size": "25.MB",
    "save_trimmed": false,
    "bracken_precision": "S",
    "min_trimmed_reads": 10000,
    "skip_fastqc": false,
    "pseudo_aligner_kmer_size": 31,
    "deseq2_vst": true,
    "umitools_extract_method": "string",
    "validate_params": true,
    "bam_csi_index": false,
    "extra_fastp_args": "--cut_mean_quality 30 --trim_front1 10 --trim_front2 10",
    "with_umi": false,
    "skip_qc": false,
    "version": false,
    "trimmer": "fastp",
    "publish_dir_mode": "copy",
    "input": "s3://steveyost-seqera/Cgigas_ArredondoEspinoza2023/fetchngs/samplesheet/samplesheet_for_rnaseq.csv",
    "skip_bigwig": false,
    "igenomes_base": "s3://ngi-igenomes/igenomes/",
    "skip_gtf_transcript_filter": false,
    "stringtie_ignore_gtf": false,
    "skip_umi_extract": false,
    "kallisto_quant_fraglen_sd": 200,
    "save_unaligned": false,
    "featurecounts_feature_type": "exon",
    "multiqc_title": "rnaseq_Cgigas_ArredondoEspinoza2023",
    "skip_gtf_filter": false,
    "custom_config_version": "master",
    "fasta": "https://steveyost-seqera.s3.us-east-1.amazonaws.com/Cgigas_ArredondoEspinoza2023/run_ncbi_dataset_GCA_000297895.1_GCF/GCF_000297895.1_oyster_v9_genomic.fna",
    "hisat2_build_memory": "200.GB",
    "skip_rseqc": false,
    "save_non_ribo_reads": false,
    "save_kraken_unassigned": false,
    "skip_alignment": false,
    "skip_qualimap": false,
    "star_ignore_sjdbgtf": false,
    "skip_stringtie": false,
    "gtf_extra_attributes": "gene_name",
    "save_merged_fastq": false,
    "save_umi_intermeds": false,
    "save_kraken_assignments": false,
    "show_hidden": false,
    "min_mapped_reads": 5,
    "rseqc_modules": "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication",
    "gtf_group_features": "gene_id",
    "skip_trimming": false,
    "kallisto_quant_fraglen": 200,
    "igenomes_ignore": false,
    "skip_pseudo_alignment": true,
    "outdir": "s3://steveyost-seqera/Cgigas_ArredondoEspinoza2023/run_ncbi_dataset_GCA_000297895.1_GCF/",
    "skip_bbsplit": true,
    "pipelines_testdata_base_path": "https://raw.githubusercontent.com/nf-core/test-datasets/7f1614baeb0ddf66e60be78c3d9fa55440465ac8/",
    "help": false,
    "stranded_threshold": 0.8,
    "unstranded_threshold": 0.1,
    "umitools_grouping_method": "directional",
    "skip_biotype_qc": false
}
```

With the resulting `salmon.merged.gene_counts.tsv` files from each run in hand, I then queried `ChatGPT  4o` for recommendations on visually comparing gene counts. Its result included a Venn diagram to indicate matching/non-matching gene IDs and scatter plots showing gene counts on X and Y axes for matching gene IDs for the respective runs.

<figure>
  <img src="https://github.com/user-attachments/assets/cc6ac25b-f86f-4424-a91d-333613750c50" alt="Venn diagram"/>
  <figcaption class="caption">Venn diagram showing common and non-common gene IDs</figcaption>
</figure>

![gene_count_scatter_plot_1_sample](https://github.com/user-attachments/assets/13213250-afa1-4d6e-8bc5-dcdbf4bcb332)
![gene_count_scatter_plot_all_samples](https://github.com/user-attachments/assets/5a5a01db-59db-438f-91e3-865afd260f5f)





