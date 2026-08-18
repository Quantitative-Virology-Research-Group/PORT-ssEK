# PORT-ssEK: a derivative of PORT-EK and PORT-EK-v2
Pathogen Origin Recognition Tool using strand-specific Enriched *k*-mers
version 1.0

A tool for the identification of *strand-asymmetric* genomic regions of viruses that arise in different multi-genomes.
Based on *k*-mer counting, does not require MSA.
Can highlight changes independent of viral phylogeny.

Citation\
[]

## Installation

**Requirements:** Python 3.12 or higher

1. Clone this repository
2. Navigate to the repository directory `cd PORT-EK2`
2. Ensure you have Python 3.12 installed (consider using [pyenv](https://github.com/pyenv/pyenv) for version management)
3. Create a Python virtual environment: `python3.12 -m venv env`
4. Activate the environment: `source env/bin/activate` (Linux/Mac) or `env\Scripts\activate` (Windows)
5. Install PORTEK: `pip install -e .`

For detailed installation instructions, see [INSTALL.md](INSTALL.md).

## Usage Guide

### Quick Start

1. **Create a new project:**
   ```bash
   portek new $project_directory
   ```

2. **Configure your project:**
   - Edit `$project_directory/config.yaml` (see Configuration section below)
   - Copy your FASTA files to `$project_directory/input/`

3. **Find optimal k-mer length and get k-mers:**
   ```bash
   portek find_k $project_directory --max_k 31
   ```

4. **Identify enriched k-mers:**
   ```bash
   portek find_enriched $project_directory -k 15
   ```

5. **Map k-mers to reference sequence:**
   ```bash
   portek map $project_directory -k 15
   ```

6. **Construct phylogenetic tree:**
   ```bash
   portek tree $project_directory -k 15
   ```

### Available Commands

#### `portek new`
Create a new PORTEK project directory structure.

```bash
portek new PROJECT_DIR
```

**Arguments:**
- `PROJECT_DIR`: Path to the new project directory (must not exist)

Creates the directory structure and a template `config.yaml` file.

---

#### `portek find_k`
Test different k-mer lengths to find the optimal k value for your dataset.

```bash
portek find_k PROJECT_DIR [OPTIONS]
```

**Arguments:**
- `PROJECT_DIR`: Path to the project directory

**Options:**
- `--min_k K`: Minimum k-mer length to test (default: 5)
- `--max_k K`: Maximum k-mer length to test (default: 31)
- `--n_jobs N`: Number of parallel processes (default: 4)
- `-v, --verbose`: Print detailed progress information

Tests all odd k values from `min_k` to `max_k` and identifies the optimal length based on k-mer entropy and uniqueness.

---

#### `portek find_enriched`
Identify k-mers enriched in specific sample groups.

```bash
portek find_enriched PROJECT_DIR -k K [OPTIONS]
```

**Arguments:**
- `PROJECT_DIR`: Path to the project directory
- `-k K`: K-mer length to analyze (required)

**Options:**
- `--fdr`: Use false discovery rate correction with Benjamini-Yekutieli procedure
- `-m, --max_mem GB`: Maximum memory (in GB) for count matrix (default: 2.0)
- `--n_jobs N`: Number of parallel processes (default: 4)
- `-v, --verbose`: Print detailed progress information

Performs statistical analysis to identify k-mers that are significantly enriched in specific groups. Outputs include volcano plots, PCA visualization, and k-mer statistics.

---

#### `portek map`
Map enriched k-mers to a reference genomic sequence.

```bash
portek map PROJECT_DIR -k K [OPTIONS]
```

**Arguments:**
- `PROJECT_DIR`: Path to the project directory
- `-k K`: K-mer length to map (required)

**Options:**
- `-d D`: Maximum edit distance for mapping (default: 2)
- `-v, --verbose`: Print detailed progress information

Maps k-mers to reference sequence positions, accounting for mismatches, insertions, and deletions up to the specified edit distance.

---

#### `portek tree`
Construct a phylogenetic tree based on k-mer count profiles.

```bash
portek tree PROJECT_DIR -k K [OPTIONS]
```

**Arguments:**
- `PROJECT_DIR`: Path to the project directory
- `-k K`: K-mer length to use (required)

**Options:**
- `-n, --tree_subsample_size N`: Number of samples to subsample (default: None)
- `-b, --balance_groups`: Balance group sizes when subsampling
- `-v, --verbose`: Print detailed progress information

Constructs a neighbor-joining phylogenetic tree using Euclidean distances between k-mer count profiles. Outputs a Newick format tree file.

---

### Configuration

Before running analysis tools, you must edit the `config.yaml` file in your project directory. Here's what each field means:

```yaml
sample_groups: [group1, group2, group3]
```
**Required.** List of sample group names. These identify the different categories in your dataset (e.g., different host species, viral subtypes, etc.). Names are case-sensitive.

```yaml
input_files: [group1.fasta, group2.fasta, group3.fasta]
```
**Required.** List of FASTA files containing sequences for each group. The order must match `sample_groups`. Files should be placed in `$project_dir/input/`.

```yaml
header_format: [plain, gisaid, ncbi]
```
**Optional.** Format of FASTA headers in each input file. Options:
- `plain`: Use entire header as sequence ID (default)
- `gisaid`: Extract GISAID accession from header (format: `...||accession|...`)
- `ncbi`: Extract NCBI accession from header (format: `accession|...`)

If omitted, assumes `plain` for all files.

```yaml
mode: ava
```
**Required.** Comparison mode:
- `ava`: All-vs-all comparison (compares every group against every other group)
- `ovr`: One-vs-rest comparison (compares one group of interest against all others combined)

```yaml
goi: group1
```
**Required only for `ovr` mode.** Specifies the group of interest. Must match one of the names in `sample_groups`.

```yaml
ref_seq: reference.fasta
```
**Optional but necessary for mapping.** Name of FASTA file containing the reference genomic sequence. Required for the `map` command. File should be placed in `$project_dir/input/`.

#### Example Configuration

```yaml
sample_groups: [human, pig, bat]
input_files: [human_sequences.fasta, pig_sequences.fasta, bat_sequences.fasta]
header_format: [ncbi, ncbi, ncbi]
mode: ava
goi: 
ref_seq: reference_genome.fasta
```

---

### Output Files

PORTEK generates various output files in `$project_dir/output/`:

- **find_k**: `k_selection_results.csv` - Statistics for each tested k value
- **find_enriched**:
  - `enriched_Kmers_stats.csv` - Statistics for enriched k-mers
  - `Kmer_counts_for_classifier.csv` - K-mer count matrix
  - `Kmer_PCA.svg` - PCA visualization
  - `volcano_plot_*.svg` - Volcano plots for each group comparison
- **map**:
  - `mapping_Kmers_max_D_mismatches.tsv` - K-mer positions on reference
  - `coverage_Kmers_max_D_mismatches.tsv` - Coverage statistics
- **tree**: `Kmer_phylo_tree.nwk` - Phylogenetic tree in Newick format
