# Supplementary Computational Methods for 5S and 45S rDNA Analyses

This repository contains the scripts, result tables, representative sequences, and figure files used to summarize 5S and 45S rDNA features in two haplotype-resolved assemblies (`hap1` and `hap2`).

The repository is organized so that the content can be rerun locally if the same assembly FASTA files and annotation inputs are available.

## Repository layout

- `methods/`
  - Supplementary methods text for the 5S and 45S analyses.
- `scripts/`
  - Python scripts used to extract loci, summarize counts, and draw figures.
- `tests/`
  - Small `pytest` tests for key parsing and plotting logic.
- `figures/`
  - Final figure files in `svg` and `png` format.
- `tables/`
  - Summary tables and detailed text summaries.
- `sequences/`
  - Representative 5S and 45S FASTA sequences used or extracted in the analysis.

## External input files required for rerunning

The following large input files were used locally but are not included in this repository:

- `fojia.hap1.filtered.fa`
- `fojia.hap2.filtered.fa`
- `hap1.rfam.gff`
- `hap2.rfam.gff`
- `hap1.rRNA.gff`
- `hap2.rRNA.gff`
- `tmp/hap1_vs_shared_5S.blast.tsv`
- `tmp/hap2_vs_shared_5S.blast.tsv`

These files contain the haplotype assemblies, the Infernal/Rfam-based ncRNA annotation results, the barrnap rRNA annotation results, and the BLAST tabular outputs used by the summary scripts.

## Software used

The analyses in this repository depend on the following software:

- Infernal / `cmscan` for Rfam-based ncRNA annotation of 5S candidates
- Rfam covariance model database
- BLAST+ for sequence similarity searches
  - verified local version used for downstream filtering: `blastn 2.5.0+`
- barrnap for rRNA annotation of 18S, 5.8S, and 28S features
- Python for custom parsing and plotting
  - verified local version: `Python 3.9.21`
- `pytest` for lightweight script verification
  - verified local version: `pytest 8.4.2`
- `rsvg-convert` for converting `svg` figures to `png`

The retained annotation files clearly record that Infernal/Rfam annotation was performed with `cmscan` and that rRNA features were derived from barrnap output. However, the exact executable versions of `cmscan` and barrnap were not embedded in the retained local result files and are therefore not claimed here.

## Analysis summary

### 1. 5S rDNA

The 5S analysis used the most abundant completely identical 5S sequence shared by `hap1` and `hap2` as the representative query:

- `sequences/shared_most_common_5S.fa`

This 119-bp sequence was searched against the haplotype assemblies with BLASTN. Hits were retained only when:

- sequence identity was `>= 98%`
- query coverage was `100%`

The retained high-confidence 5S-like loci were summarized and plotted for the chromosomes with the dominant signals (`chr2` and `chr12` in both haplotypes).

Main outputs:

- `tables/main_table_5S_98id100cov.tsv`
- `figures/figure_5S_98id100cov_chr2_chr12.svg`
- `figures/figure_5S_98id100cov_chr2_chr12.png`

### 2. 45S rDNA

Complete 45S units were defined from barrnap-derived annotations by requiring co-oriented `18S_rRNA`, `5_8S_rRNA`, and `28S_rRNA` features on the same sequence within a 10-kb span. In this dicot analysis, the barrnap `28S_rRNA` feature was treated as the 26S-equivalent LSU annotation.

Major chromosome-level 45S arrays were summarized on `chr6`, `chr8`, and `chr16` in both haplotypes.

Main outputs:

- `tables/hap1_45S_units.tsv`
- `tables/hap2_45S_units.tsv`
- `tables/hap1_45S_summary.tsv`
- `tables/hap2_45S_summary.tsv`
- `tables/45S_summary.txt`
- `figures/figure_45S_chr16_chr8_chr6.svg`
- `figures/figure_45S_chr16_chr8_chr6.png`

### 3. Representative 45S sequences

Representative exact 45S sequence types were extracted from the major chromosome-level arrays, together with the dominant exact 45S type shared by both haplotypes.

Outputs:

- `sequences/representative_45S.fa`
- `tables/representative_45S_summary.txt`
- `sequences/shared_representative_45S.fa`
- `tables/shared_representative_45S_summary.txt`

## Scripts included

- `scripts/build_5s_98_outputs.py`
  - Filters BLAST hits at `identity >= 98%` and `coverage = 100%`, summarizes 5S loci, and renders the main 5S figure.
- `scripts/find_45s_units.py`
  - Parses barrnap GFF files and detects complete 45S units using fixed structural rules.
- `scripts/plot_45s_chr6_8_16.py`
  - Draws the chromosome-scale 45S distribution figure.
- `scripts/extract_representative_45s.py`
  - Extracts the dominant exact 45S type in each haplotype.
- `scripts/extract_shared_representative_45s.py`
  - Extracts the dominant exact 45S type shared by both haplotypes.

## Example commands for reproducing the repository outputs

The commands below assume the repository root is the working directory and that the required external input files are available at the paths expected by the scripts.

### Rebuild 5S summary outputs

```bash
python scripts/build_5s_98_outputs.py
rsvg-convert figures/figure_5S_98id100cov_chr2_chr12.svg \
  -o figures/figure_5S_98id100cov_chr2_chr12.png
```

### Recompute 45S units and summaries

```bash
python scripts/find_45s_units.py
python scripts/plot_45s_chr6_8_16.py
rsvg-convert figures/figure_45S_chr16_chr8_chr6.svg \
  -o figures/figure_45S_chr16_chr8_chr6.png
```

### Re-extract representative 45S sequences

```bash
python scripts/extract_representative_45s.py
python scripts/extract_shared_representative_45s.py
```

### Run included tests

```bash
pytest -q tests/test_build_5s_98_outputs.py \
  tests/test_find_45s_units.py \
  tests/test_plot_45s_chr6_8_16.py \
  tests/test_extract_representative_45s.py \
  tests/test_extract_shared_representative_45s.py
```

## Notes on paths

The current scripts were written from the original local analysis workspace and therefore use absolute local base paths internally. If the repository is rerun in a different environment, the `BASE` and `TMP` path definitions at the top of each script should be adjusted to match the local directory layout.

## Methods text

The repository also includes a prose supplementary methods file:

- `methods/github_supplementary_5S_45S_methods.md`

This file is intended for manuscript or repository documentation, whereas the scripts and tables are the executable and derived components used to regenerate the summarized results.
