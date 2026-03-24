# Ribosomal DNA Supplementary Package

This package contains the scripts, figures, tables, and sequence files used for the 5S and 45S rDNA analyses in the haplotype-resolved assemblies.

## Contents

- `README.md`
  Main entry point for the package.
- `methods/`
  Supplementary computational methods in English.
- `scripts/`
  Analysis scripts used to derive the 5S and 45S results.
- `figures/`
  Final PNG and SVG figures for 5S and 45S distributions.
- `tables/`
  Main summary tables and chromosome-level summaries.
- `sequences/`
  Representative and shared 5S/45S FASTA files.

## Main files

- `methods/github_supplementary_5S_45S_methods.md`
- `figures/figure_5S_98id100cov_chr2_chr12.png`
- `figures/figure_45S_chr16_chr8_chr6.png`
- `tables/main_table_5S_98id100cov.tsv`
- `tables/hap1_45S_summary.tsv`
- `tables/hap2_45S_summary.tsv`
- `sequences/shared_most_common_5S.fa`
- `sequences/representative_45S.fa`
- `sequences/shared_representative_45S.fa`

## Notes

- The 5S analysis used the shared predominant 5S sequence identified from Infernal/Rfam annotations.
- The 45S analysis used barrnap rRNA annotations and treated `28S_rRNA` as the 26S-equivalent LSU feature for dicot 45S rDNA.
- This directory is organized for convenient upload to GitHub as supplementary material.
