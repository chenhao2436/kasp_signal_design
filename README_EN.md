# KASP Signal Design Workflow

**简体中文：** [README.md](README.md) | **English**

The **KASP Signal Design Workflow** prioritizes candidate SNPs from BSA results and prepares their flanking sequences for downstream KASP assay design. It connects SNP-index ranking, VCF genotype filtering, candidate reporting, and reference-sequence extraction. Snakemake is used to coordinate region-level jobs.

> **Status: research workflow / alpha.** The repository contains Snakemake rules and analysis scripts, but configuration files still include server-specific absolute paths. Replace them and validate the workflow with your own data before deployment. The legacy Bash driver expects an INI file that is not included in the repository and is not recommended as the default entry point.

## Features

- Rank candidate SNPs on the left and right sides of a BSA peak region using `DELTA_SNP_INDEX`.
- Support flank-based and bin-based candidate selection.
- Inspect parental and pool genotypes in a joint VCF.
- Filter non-SNPs, multiallelic sites, non-homozygous or non-differing parental genotypes, and nearby variants.
- Prepare recommended-marker tables and input for sequence extraction.
- Extract flanking reference sequence and mark the variant as `[REF/ALT]`.
- Run multiple target regions through Snakemake.

## Workflow

```text
SNP-index / DELTA_SNP_INDEX table
              |
              v
    01_rank_snpindex_candidates
              |
              v
       ranked candidate SNPs
              |
              +---- Joint VCF + sample genotypes
              |                 |
              v                 v
       02_filter_isolated_parental_snps
              |
              v
         PASS / FAIL sites
              |
              v
     03_prepare_extract_fa_input
              |
              v
       extract_fa.py + reference FASTA
              |
              v
     KASP flanking sequence table
```

## Repository layout

```text
kasp_signal_design/
├── Snakefile
├── environment.yml             # Snakemake driver environment
├── data/
│   ├── *_config.yaml           # Analysis configurations
│   ├── *_regions.tsv           # Target regions
│   └── kasp_rule_env.yaml      # Per-rule Conda environment
├── scripts/
│   ├── 20260418_kasp_01_rank_snpindex_candidates.py
│   ├── 20260418_kasp_02_filter_isolated_parental_snps.py
│   ├── 20260418_kasp_03_prepare_extract_fa_input.py
│   ├── extract_fa.py
│   └── 20260418_kasp_00_run_pipeline.sh  # Legacy Bash driver
├── doc/
├── results/
└── temp/
```

## Requirements and installation

The workflow uses two Conda environments.

### Snakemake driver

- Linux (recommended for server use).
- Conda or Mamba.
- Snakemake.

Create the driver environment:

```bash
conda env create -f environment.yml
conda activate kasp-snakemake
snakemake --version
```

### Rule environment

The per-rule environment is defined in `data/kasp_rule_env.yaml` and includes:

- Python 3.12
- bcftools
- htslib
- pyfaidx

Snakemake creates this environment when run with `--use-conda`. Most ranking, filtering, and table-preparation code uses the Python standard library. VCF genotype queries call `bcftools`, and reference sequence extraction uses `pyfaidx`.

A separate pip `requirements.txt` is not necessary for the current implementation because dependencies are already separated into the Snakemake driver and rule environments.

## Input data

Each analysis configuration provides:

| Key | Description |
|---|---|
| `snpindex_path` | BSA SNP-index / plot result table |
| `vcf_path` | bgzip-compressed and indexed joint VCF |
| `reference_fasta` | Reference FASTA corresponding to the VCF |
| `parent1_sample` / `parent2_sample` | Parent sample names in the VCF |
| `high_pool_sample` / `low_pool_sample` | High- and low-pool sample names |
| `flank_no_variant_bp` | Nearby-variant exclusion distance, commonly 50 bp |
| `target_markers_per_side` | Target number of markers per side |
| `initial_rank_limit` | Initial candidate inspection limit |
| `expansion_step_bp` | Expansion step toward the peak center |
| `regions_file` | Region ID to genomic interval table (optional; see Snakefile default) |

The SNP-index table should contain columns such as `CHROM`, `POS`, and `DELTA_SNP_INDEX`. Headerless input is interpreted using fixed default column positions. The VCF must contain the configured sample names and have a usable `.tbi` or `.csi` index. Chromosome naming must be consistent across the SNP-index table, VCF, and reference FASTA.

**Important:** The committed configuration files retain absolute paths and sample names from the original server environment. Review and replace these values before running the workflow.

## Quick start

### 1. Review the active configuration

The current `Snakefile` loads this configuration by default:

```text
data/70_Chr11_255900000_260500000_config.yaml
```

Edit the file before running. Update the SNP-index path, VCF path, reference FASTA, sample names, and target regions.

To use another region configuration, update the `configfile` directive at the top of the Snakefile and confirm that its regions TSV exists.

### 2. Perform a dry run

```bash
snakemake -n
```

Review the planned rules and output paths. A dry run does not execute analysis commands.

### 3. Run the workflow

```bash
snakemake --use-conda -j 1
```

Start with one job. Increase parallelism only after confirming that paths, indexes, and VCF sample names are correct.

### 4. Outputs

For each region, files are written to `results/<region_id>/`:

| File | Description |
|---|---|
| `01_ranked_candidates.tsv` | Ranked candidate SNPs for the left and right sides |
| `01_rank_summary.tsv` | Candidate counts, windows, and expansion levels |
| `02_checked_sites.tsv` | All inspected sites and filtering reasons |
| `02_pass_sites.tsv` | Candidates passing hard filters |
| `03_extract_fa_input.tsv` | Headerless four-column input for sequence extraction |
| `03_recommended_markers.tsv` | Recommended markers and shortage reasons |
| `04_kasp_sequences.tsv` | Flanking sequences with REF/ALT annotation |

Run logs are written to `temp/<region_id>/`.

## Filtering logic

The current candidate filter checks:

1. Whether the candidate position exists in the VCF.
2. Whether the site is a biallelic SNP rather than an indel or multiallelic site.
3. Whether both parents are homozygous.
4. Whether the two parents have different genotypes.
5. Whether genotypes can be read for the high and low pools.
6. Whether other variants occur within the configured flank around the candidate.

The current example configuration commonly uses a 50 bp nearby-variant window and targets two markers per side. The YAML file is authoritative.

**A passing candidate is not an experimentally validated KASP assay.** The generated sequences are input material for primer design and subsequent laboratory validation.

## Known limitations

- Configuration files contain server-specific absolute paths.
- The Snakefile currently defaults to a specific Chr11 configuration rather than a universal region set.
- The legacy `scripts/20260418_kasp_00_run_pipeline.sh` expects `20260418_kasp_pipeline.ini`, which is not present in the repository. Use the Snakemake entry point instead.
- No small, distributable test VCF, SNP-index table, or reference genome is included. Cross-server end-to-end validation is therefore not claimed.
- The sequence extractor tries alternative chromosome naming formats, but the reference assembly must still match the VCF.
- Primer specificity, amplification efficiency, and genotype-cluster quality require separate assay-design and experimental validation.

## Validation checklist

For each run, check:

- `snakemake -n` shows the expected regions and outputs.
- The VCF is readable, indexed, and contains the four configured samples.
- Failure reasons in `02_checked_sites.tsv` are biologically and technically plausible.
- `03_recommended_markers.tsv` meets the target count on each side.
- The reference allele in each extracted sequence matches the reference FASTA.

## Roadmap

- [ ] Replace server-specific paths with reusable configuration templates.
- [ ] Unify parameter definitions between Snakemake and the legacy Bash entry point.
- [ ] Add a small synthetic dataset and automated tests.
- [ ] Add configuration validation and clearer run reports.
- [ ] Add optional downstream checks for KASP primer design results.

## Citation and support

When using this research workflow, cite the repository URL, the Git commit, and the actual software versions used.

For issues, provide the configuration with private paths removed, the command, relevant log excerpts, and input format details.
