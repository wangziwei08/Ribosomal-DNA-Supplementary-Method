# Supplementary Computational Methods for 5S and 45S rDNA Analyses

## Overview

This document summarizes the computational workflow used to characterize 5S and 45S ribosomal DNA (rDNA) in the two haplotype-resolved assemblies (`hap1` and `hap2`). The workflow combines covariance model-based ncRNA annotation, BLAST-based sequence profiling, rule-based identification of complete 45S units, and custom plotting scripts used to generate the summary figures included in this repository.

All analyses were performed in the directory `/home/dell/disk2_16T/fojia/Rfam`, and all derived intermediate and final output files described below were written to the local `tmp/` directory.

## 1. 5S rDNA analysis

### 1.1 Identification of the representative 5S query sequence

Initial ncRNA annotation was performed with Infernal using Rfam covariance models. From the high-confidence 5S rRNA annotations, we identified the most abundant completely identical 5S sequence that was shared between `hap1` and `hap2`. This shared dominant sequence was saved as:

- `shared_most_common_5S.fa`

This sequence is 119 bp long and was used as the representative 5S query for downstream genome-wide sequence profiling.

### 1.2 Genome-wide 5S profiling

The representative 5S query was searched against the haplotype assemblies with local BLASTN. For the main analysis, hits were retained only when they satisfied both of the following criteria:

- sequence identity `>= 98%`
- query coverage `= 100%`

These filtered matches were treated as high-confidence 5S-like loci. Although the query was searched against the whole assemblies, the majority of retained loci were concentrated on:

- `hap1`: `chr2_1` and `chr12_1`
- `hap2`: `chr2_2` and `chr12_2`

The final high-confidence copy numbers were:

- `hap1 chr2_1`: 591
- `hap1 chr12_1`: 1,687
- `hap2 chr2_2`: 1,875
- `hap2 chr12_2`: 1,275

### 1.3 5S clustering and plotting

High-confidence 5S hits were further grouped into local clusters by merging neighboring hits separated by `<= 1 kb`. For figure generation, hit density was summarized in fixed `0.25 Mb` windows, and chromosome coordinates were displayed with `1 Mb` tick marks.

The main scripts used for the 5S analysis were:

- `tmp/build_5s_98_outputs.py`
- `tmp/test_build_5s_98_outputs.py`

The main output files were:

- `tmp/main_table_5S_98id100cov.tsv`
- `tmp/figure_5S_98id100cov_chr2_chr12.svg`
- `tmp/figure_5S_98id100cov_chr2_chr12.png`
- `tmp/scientific_data_methods_results_5S_98id100cov.md`

## 2. 45S rDNA analysis

### 2.1 Annotation source and interpretation of LSU labels

Large-subunit and small-subunit rRNA annotations were obtained from barrnap output:

- `hap1.rRNA.gff`
- `hap2.rRNA.gff`

In these barrnap annotations, the large-subunit rRNA feature is reported as `28S_rRNA`. For the present dicot genome analysis, this feature was treated as the 26S-equivalent large-subunit annotation when defining complete 45S units.

### 2.2 Rule-based detection of complete 45S units

Complete 45S units were identified with a custom rule-based parser. A locus was counted as one complete 45S unit when all of the following conditions were met:

- the annotations `18S_rRNA`, `5_8S_rRNA`, and `28S_rRNA` occurred on the same sequence
- all three annotations were on the same strand
- the three annotations were contained within a local span of `<= 10 kb`
- the 5.8S annotation was positioned between the 18S and large-subunit annotations
- overlap between the annotated 5.8S and the large-subunit boundary was allowed, because this pattern occurs in barrnap output

This implementation allowed both forward and reverse strand configurations, provided that the three component annotations remained co-oriented on the same strand.

### 2.3 Chromosomal distribution of complete 45S units

Complete 45S units were found on multiple sequences, but the major chromosome-level arrays were concentrated on:

- `hap1`: `chr6_1`, `chr8_1`, and `chr16_1`
- `hap2`: `chr6_2`, `chr8_2`, and `chr16_2`

The corresponding unit counts on these chromosomes were:

- `hap1 chr6_1`: 46
- `hap1 chr8_1`: 47
- `hap1 chr16_1`: 50
- `hap2 chr6_2`: 48
- `hap2 chr8_2`: 50
- `hap2 chr16_2`: 69

For visualization, complete 45S units were summarized in `0.1 Mb` windows. The plotting script used a shared chromosome-length scale across all displayed panels so that chromosome bar lengths remained proportional to their actual lengths.

The main scripts used for the 45S analysis were:

- `tmp/find_45s_units.py`
- `tmp/test_find_45s_units.py`
- `tmp/plot_45s_chr6_8_16.py`
- `tmp/test_plot_45s_chr6_8_16.py`

The main output files were:

- `tmp/hap1_45S_units.tsv`
- `tmp/hap2_45S_units.tsv`
- `tmp/hap1_45S_summary.tsv`
- `tmp/hap2_45S_summary.tsv`
- `tmp/45S_summary.txt`
- `tmp/figure_45S_chr16_chr8_chr6.svg`
- `tmp/figure_45S_chr16_chr8_chr6.png`

## 3. Representative 45S sequences

Representative complete 45S sequences were extracted from the major chromosome-level arrays (`chr6`, `chr8`, and `chr16`) in each haplotype. Exact full-length sequence types were counted after strand normalization, and the most abundant exact 45S sequence in each haplotype was exported as a representative sequence.

The script used for this step was:

- `tmp/extract_representative_45s.py`

Outputs:

- `tmp/representative_45S.fa`
- `tmp/representative_45S_summary.txt`

We also identified the most abundant exact 45S sequence shared by both haplotypes across these major chromosome-level arrays.

The script used for the shared representative sequence was:

- `tmp/extract_shared_representative_45s.py`

Outputs:

- `tmp/shared_representative_45S.fa`
- `tmp/shared_representative_45S_summary.txt`

## 4. Testing and verification

Small unit-style tests were written for the custom parsing and plotting scripts to verify critical behaviors, including:

- midpoint-based bin assignment
- strand-aware sequence extraction
- exact shared-sequence selection
- fixed y-axis labels used in the figures
- compatibility of complete 45S detection with barrnap-style overlapping 5.8S/LSU boundaries

The corresponding test files are:

- `tmp/test_build_5s_98_outputs.py`
- `tmp/test_find_45s_units.py`
- `tmp/test_plot_45s_chr6_8_16.py`
- `tmp/test_extract_representative_45s.py`
- `tmp/test_extract_shared_representative_45s.py`

At the time of generation, all local tests for these scripts passed successfully.
