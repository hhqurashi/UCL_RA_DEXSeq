# UCL RA DEXSeq Analysis

This repository contains scripts used for downstream exon-usage analysis with `DEXSeq`.

The workflow is designed for alternative splicing / differential exon usage analysis from RNA-seq data. It includes setup steps for preparing the environment, checking annotation files, generating DEXSeq-compatible flattened annotation, creating exon-bin count files, and running DEXSeq comparisons.

## Repository contents

```text
1_environment.sh
2_set_paths_and_check_gtf.sh
3_prepare_DEXSeq_flattened_annotation.sh
4_merge_lanes_into_one_BAM_per_sample.sh
5_build_sample_table.sh
6_generate_DEXSeq_exon_bin_count_files.sh
7_DEXSeq_comparisons.R
dexseq_complete_script_all_in_one.sh
dexseq_complete_script_all_in_one_FINAL.sh
export_dexseq_exons_from_rds.R
run_DEXSeq_4_comparisons.R
run_DEXSeq_4_comparisons_QC_passed.R
```

## Workflow summary

The analysis workflow broadly consists of:

1. Setting up the software environment
2. Defining input/output paths and checking the GTF annotation
3. Preparing a flattened annotation file for DEXSeq
4. Merging lane-level BAM files into one BAM per sample
5. Building a DEXSeq sample table
6. Generating exon-bin count files
7. Running DEXSeq differential exon usage comparisons
8. Exporting exon-level results from saved DEXSeq objects

## Main script groups

### Setup and annotation preparation

```text
1_environment.sh
2_set_paths_and_check_gtf.sh
3_prepare_DEXSeq_flattened_annotation.sh
```

These scripts set up the analysis environment, define paths, check annotation files, and generate the flattened GFF/GTF-style annotation required by DEXSeq.

### BAM and sample preparation

```text
4_merge_lanes_into_one_BAM_per_sample.sh
5_build_sample_table.sh
```

These scripts prepare sample-level BAM files and generate the sample table used for DEXSeq input.

### Exon-bin counting

```text
6_generate_DEXSeq_exon_bin_count_files.sh
```

This script generates DEXSeq-compatible exon-bin count files from aligned RNA-seq BAM files.

### DEXSeq analysis

```text
7_DEXSeq_comparisons.R
run_DEXSeq_4_comparisons.R
run_DEXSeq_4_comparisons_QC_passed.R
```

These R scripts run differential exon usage analysis across defined comparisons and export result tables.

### All-in-one workflow scripts

```text
dexseq_complete_script_all_in_one.sh
dexseq_complete_script_all_in_one_FINAL.sh
```

These scripts combine multiple workflow steps into a single runnable pipeline-style script.

### Result export

```text
export_dexseq_exons_from_rds.R
```

This script exports exon-level DEXSeq results from saved RDS objects for downstream review and reporting.

## Input

Expected inputs include:

- aligned RNA-seq BAM files
- GTF genome annotation
- sample metadata
- DEXSeq flattened annotation files
- exon-bin count files

Raw FASTQ files, BAM files, and large pipeline outputs are not tracked in this repository.

## Output

Outputs include:

- merged BAM file references
- DEXSeq sample tables
- exon-bin count files
- DEXSeq result tables
- saved RDS analysis objects
- exported exon-level summaries
- HTML reports or downstream plots, depending on script settings

Large results and intermediate files should remain outside GitHub.
