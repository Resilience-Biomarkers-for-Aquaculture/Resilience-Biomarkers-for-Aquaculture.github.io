---
layout: post
title: PCA for gene count comparison between 2012 and 2023 reference genomes (post in progress)
tags: rnaseq pca 
---

## 12-23-2024

The goal here was to compare gene counts between `rnaseq` runs using
1. [The NCBI reference genome and gene annotations from the period of Roberto's study](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000297895.1/) (herein referred to as `2012`)
2. The most recent [reference genome and gene annotations](https://ncbi.nlm.nih.gov/datasets/genome/GCA_963853765.1/) (herein referred to as `2023`) for C. gigas.

## Motivation
As previously [noted](https://resilience-biomarkers-for-aquaculture.github.io/SW-CGI_ID_matching_attemp_1/), in attemping to compare gene counts between
Roberto's results and `rnaseq` results, we encountered the following hurdle:
The gene count matrix that Roberto provided predominantly used gene IDs prefixed with `MSTRG_` (novel or unannotated genes/transcripts identified by
StringTie) and `CGI_` (CpG islands), whereas the available NCBI annotations, and thus gene count data, for `2023` uses gene IDs prefixed with `LOC_`,
commonly used for annotated or predicted genes that don't yet have an official gene symbol or name.
However, we found that the `2012` NCBI genome and annotations did use `LOC_` gene IDs, which facilitates matching by ID
rather than requiring BLAST or other indirect gene ID matching methods.

So, for this current task, we used the `2012`
reference genome (linked above) as as a proxy for Roberto's, supported by the fact that it represents the state of the art at the time of his
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

With the resulting `salmon.merged.gene_counts.tsv` files from each run in hand, I then queried `ChatGPT  4o` for recommendations on visually comparing gene counts. Its recommendations included a Venn diagram to indicate matching/non-matching gene IDs and scatter plots showing gene counts on X and Y axes for matching gene IDs for the respective runs.

Further, to visualize the gene count differences between the `rnaseq` runs using the `2012` and `2023` reference genomes and annotations, I chose PCA among the recommendations made by `ChatGPT 4o`, and interactively worked with it to produce a 2D PCA plot which, with additional metadata for the samples, differentiates the two runs (by color), and each sample's noted resilience, and time interval (day 0 vs. day 30). Satisfied with this, and noting the striking similarity between the two runs, I then had ChatGPT produce an interactive 3D PCA plot to see whether greater differences appeared given a third principal component axis. 

## Results

<figure>
  <img src="https://github.com/user-attachments/assets/cc6ac25b-f86f-4424-a91d-333613750c50" alt="Venn diagram"/>
  <figcaption class="caption">Venn diagram showing common and non-common gene IDs</figcaption>
</figure>

<figure>
    <img src="https://github.com/user-attachments/assets/13213250-afa1-4d6e-8bc5-dcdbf4bcb332" alt="Scatter plot for a single sample"/>
    <figcaption class="caption">For a single sample</figcaption>
</figure>

<figure>
    <img src="https://github.com/user-attachments/assets/5a5a01db-59db-438f-91e3-865afd260f5f" alt="Scatter plot for all samples overlaid"/>
    <figcaption class="caption">For all samples, overlaid</figcaption>
</figure>

<figure>
    <img src="https://github.com/user-attachments/assets/486d760c-346b-46ef-911c-94dc6f7dd3ba" alt="PCA plot"/>
    <figcaption class="caption">PCA of gene counts for genes IDs in common between runs, for each sample. Blue indicates `2012` reference genome/annotation, red indicates `2023`, Xs for susceptible samples, dots for resistant, smaller points for Day 0 samples and larger for Day 30.</figcaption>
</figure>

<figure>
    <img src="https://github.com/user-attachments/assets/a603cac1-0b99-4c17-8f96-e6a76ae4e7cf" alt="3D PCA plot"/>
    <figcaption class="caption">Screenshot of 3D PCA</figcaption>
</figure>

<iframe src="assets/interactive_3d_pca.html" width="100%" height="600px" frameborder="0"></iframe>

The 3rd principal component dimension reveals differences between the two runs. Given time, we will (1) display the contributions of each dimension and (2) the top contributing gene IDs for each dimension.

## Code
Below is the python script that generated the above Venn diagram and scatter plots:
```python
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import matplotlib
import numpy as np

# Use non-interactive backend for headless environments
matplotlib.use('Agg')

# Function to load gene count files (ignoring the second column)
def load_gene_counts(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df = df.drop(columns=df.columns[1])  # Drop the second column (gene_name)
    return df

# Function to prepare the comparison table
def compare_gene_counts(df1, df2):
    # Merge data on gene_id
    # The suffixes modifies the column names for the respective dataframes
    merged = pd.merge(df1, df2, on='gene_id', how='outer', suffixes=('_run1', '_run2'))

    # Mark presence in each run
    # : is the row indexer; it means select all rows
    # 1 selects the second column
    # len(df1.columns selects the first column in df2 data
    # The 'present_in' columns are appended to each row
    # Only one column of each original dataframe is checked because the outer
    # join assigns NaN to all colums if the joined-on key isn't present.
    merged['present_in_run1'] = ~merged.iloc[:, 1].isna()
    merged['present_in_run2'] = ~merged.iloc[:, len(df1.columns)].isna()

    return merged

# Function to generate visualizations
def generate_visualizations(merged_df, output_prefix, df1_columns):
    # Venn Diagram
    only_in_run1 = merged_df['present_in_run1'] & ~merged_df['present_in_run2']
    only_in_run2 = merged_df['present_in_run2'] & ~merged_df['present_in_run1']
    in_both_runs = merged_df['present_in_run1'] & merged_df['present_in_run2']

    venn2(
        subsets=(only_in_run1.sum(), only_in_run2.sum(), in_both_runs.sum()),
        set_labels=('Run 1', 'Run 2')
    )
    plt.title("Gene ID Overlap Between Runs")
    plt.savefig(f"{output_prefix}_venn.png")
    plt.clf()

    # Scatter Plot for common genes
    common_genes = merged_df[merged_df['present_in_run1'] & merged_df['present_in_run2']]

    plt.xlabel("Log2 Run 1 Counts")
    plt.ylabel("Log2 Run 2 Counts")
    plt.title(f"Scatter Plot of Log2 Gene Counts for Common Genes")
    plt.grid(True)
    for col_idx in range(1, len(df1_columns)):
        column_name = df1_columns[col_idx]
        print(f"Processing sample: {column_name}")
        sample_1_run1 = column_name + '_run1'
        sample_1_run2 = column_name + '_run2'

        # Check if columns exist before plotting
        if sample_1_run1 in common_genes.columns and sample_1_run2 in common_genes.columns:
            x = np.log2(common_genes[sample_1_run1] + 1)
            y = np.log2(common_genes[sample_1_run2] + 1)

            plt.scatter(x, y, alpha=0.3, s=5)  # Reduced dot size (s=5)
        else:
            print(f"Error: Columns {sample_1_run1} or {sample_1_run2} not found in the merged dataset.")
    plt.savefig(f"{output_prefix}_scatter.png")

# Main Function
def main(file1, file2, output_prefix):
    # Load data
    df1 = load_gene_counts(file1)
    df2 = load_gene_counts(file2)

    # Perform comparison
    merged_df = compare_gene_counts(df1, df2)

    # Save tabular comparison
    merged_df.to_csv(f"{output_prefix}_comparison.csv", index=False)

    # Generate visualizations
    generate_visualizations(merged_df, output_prefix, df1.columns)

# Example usage:
# main("run1_counts.tsv", "run2_counts.tsv", "gene_comparison")


# Example usage:
main("ncbi_dataset_GCA_000297895.1/ncbi_dataset/data/GCF_000297895.1/salmon.merged.gene_counts.tsv", "seqera_roberto/salmon.merged.gene_counts_seqera.tsv", "gene_comparison")
```






