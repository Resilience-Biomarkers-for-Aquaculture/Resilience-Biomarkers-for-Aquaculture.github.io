---
layout: post
title: Run Seqera rnaseq pipeline for C. gigas study
tags: nextflow seqera aws rnaseq
---

## 11-24-2024

The goal here is to run the `rnaseq` on Seqera, replicating the process described in <https://resilience-biomarkers-for-aquaculture.github.io/ES-RNAseq_with_reference_dataset1/>. Referring to that post, I used the Seqera pipleine input form for the `rnaseq` pipeline using the parameters described there. For the input sample sheet, I converted the output of <https://resilience-biomarkers-for-aquaculture.github.io/SY-Sequera_fetchngs_Cgigas/> to strip quotes and to include only the four columns specified. To do that, I had ChatGPT write the following python script:
```
import csv

def extract_columns(input_csv, output_csv):
    # Define the desired columns
    desired_columns = ['sample', 'fastq_1', 'fastq_2', 'strandedness']

    try:
        with open(input_csv, 'r', newline='', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            # Ensure all desired columns exist in the input
            if not all(column in reader.fieldnames for column in desired_columns):
                raise ValueError("Input CSV file does not contain all required columns.")

            with open(output_csv, 'w', newline='', encoding='utf-8') as outfile:
                writer = csv.DictWriter(outfile, fieldnames=desired_columns)
                writer.writeheader()

                for row in reader:
                    # Extract desired columns and remove quotes
                    cleaned_row = {key: row[key].strip('"').strip("'") for key in desired_columns}
                    writer.writerow(cleaned_row)

        print(f"Extraction complete. Output written to {output_csv}.")

    except FileNotFoundError:
        print(f"Error: The file {input_csv} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
input_csv = 'input.csv'  # Replace with your input file path
output_csv = 'output.csv'  # Replace with your desired output file path
extract_columns(input_csv, output_csv)
```

Below are the resulting parameters:
```
{
  "help": false,
  "rseqc_modules": "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication",
  "gtf_extra_attributes": "gene_name",
  "umitools_grouping_method": "directional",
  "gtf_group_features": "gene_id",
  "validate_params": true,
  "bracken_precision": "S",
  "ribo_database_manifest": "${projectDir}/workflows/rnaseq/assets/rrna-db-defaults.txt",
  "help_full": false,
  "custom_config_version": "master",
  "aligner": "star_salmon",
  "featurecounts_group_type": "gene_biotype",
  "igenomes_base": "s3://ngi-igenomes/igenomes/",
  "publish_dir_mode": "copy",
  "umitools_extract_method": "string",
  "show_hidden": false,
  "featurecounts_feature_type": "exon",
  "min_trimmed_reads": 10000,
  "kallisto_quant_fraglen": 200,
  "kallisto_quant_fraglen_sd": 200,
  "pipelines_testdata_base_path": "https://raw.githubusercontent.com/nf-core/test-datasets/7f1614baeb0ddf66e60be78c3d9fa55440465ac8/",
  "skip_preseq": true,
  "stranded_threshold": 0.8,
  "hisat2_build_memory": "200.GB",
  "skip_bbsplit": true,
  "deseq2_vst": true,
  "trimmer": "fastp",
  "min_mapped_reads": 5,
  "pseudo_aligner_kmer_size": 31,
  "max_multiqc_email_size": "25.MB",
  "unstranded_threshold": 0.1,
  "custom_config_base": "https://raw.githubusercontent.com/nf-core/configs/master",
  "gencode": false,
  "igenomes_ignore": false,
  "remove_ribo_rna": false,
  "with_umi": false,
  "umitools_dedup_stats": false,
  "bam_csi_index": false,
  "star_ignore_sjdbgtf": false,
  "stringtie_ignore_gtf": false,
  "save_merged_fastq": false,
  "save_umi_intermeds": false,
  "save_non_ribo_reads": false,
  "save_bbsplit_reads": false,
  "save_reference": false,
  "save_trimmed": false,
  "save_align_intermeds": false,
  "save_unaligned": false,
  "save_kraken_assignments": false,
  "save_kraken_unassigned": false,
  "skip_gtf_filter": false,
  "skip_gtf_transcript_filter": false,
  "skip_umi_extract": false,
  "skip_trimming": false,
  "skip_alignment": false,
  "skip_pseudo_alignment": true,
  "skip_markduplicates": false,
  "skip_bigwig": false,
  "skip_stringtie": false,
  "skip_fastqc": false,
  "skip_dupradar": false,
  "skip_qualimap": false,
  "skip_rseqc": false,
  "skip_biotype_qc": false,
  "skip_deseq2_qc": false,
  "skip_multiqc": false,
  "skip_qc": false,
  "version": false,
  "plaintext_email": false,
  "monochrome_logs": false,
  "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1/GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz",
  "gtf": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1/GCF_963853765.1_xbMagGiga1.1_genomic.gtf.gz",
  "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1/GCF_963853765.1_xbMagGiga1.1_genomic.gff.gz",
  "extra_fastp_args": "--cut_mean_quality 30 --trim_front1 10 --trim_front2 10",
  "multiqc_title": "rnaseq_Cgigas_ArredondoEspinoza2023",
  "input": "s3://steveyost-seqera/Cgigas_ArredondoEspinoza2023/fetchngs/samplesheet/samplesheet_for_rnaseq.csv",
  "outdir": "s3://steveyost-seqera/Cgigas_ArredondoEspinoza2023/rnaseq/"
}
```

It terminated with an error:
```
Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_FASTP:FASTQC_RAW (2)'

Caused by:
  Wave container request image cannot start with URL like prefix - offending value: https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0
```

Tried again with the `mamba` config rather than the `singularity` config.

