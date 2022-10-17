# __MGG prophages__

<div align="justify">
MGG snakemake workflow for confident prophage detection in annotated bacterial genomes.
</div> <br>

## __Steps__ (Linux)

### Clone repo & install snakemake

```sh
git clone https://github.com/bioinf-mcb/mgg_prophages
conda install -c conda-forge -c bioconda snakemake mamba biopython=1.79 pathlib=1.0.1 pandas datetime
```

### Run prophage detection

<div align="justify">
To run workflow customize config file (config.yml), download dependencies (SetupTools) and execute workflow (PROPHAGES). Make sure that bacterial and protein ids are unique! Avoid white characters in bacterial contig ids!
</div> <br>


**1. Customize config file (config.yml)**

1. Directory to download and install dependencies (SETUP_TOOLS_DIR).
2. File with annotated bacterial genome(s) (recommended: via PATRIC server) (GENBANK_FILE). <br><br>

**2. Download dependencies (SetupTools)**

<div align="justify">
To setup dependencies run command below. If you already did that it is enough to point in config file (SETUP_TOOLS_DIR) to directory with dependencies.
</div> <br>

```sh
snakemake --use-conda --cores all --snakefile SetupTools
```
<br>

**3. Run test (optional)**

<div align="justify">
After setting up dependencies (point 2) run prophage detection with default config file paths.
</div> <br>

```sh
snakemake --use-conda --cores all --snakefile PROPHAGES
```
<br>

**4. Execute prophage detection**

```sh
snakemake --use-conda --cores all --snakefile PROPHAGES
```
<br>


## Results

<div align="justify">
Sequences and metadata of prophages and prophage-like elements can be found in prophages.fasta and prophages.tsv, respectively.
To obtain high-confidence prophages one should consider phages full-filling criteria of completeness above 90% and at least medium confidence.
</div> <br>


## __Details__ (preliminary)

<div align="justify">
Tool detects prophages in bacterial genomes by using two complementary tools: VirSorter and PhiSpy. Primary prophage detection is any region detected the tools as a prophage. Subsequently, detections are decontaminated from bacterial DNA and their completeness is estimated (CheckV). Based, on completeness prophages are filtered to obtain only high-confidence prophages.

Prophage circurality is detected only for prophages found as whole contigs and it's evalueated by VirSorter.
</div> <br>





## __DEV info__


1 To solve problems with installing env by snakemake (Error: cannot read json file (?) ).

```sh
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda config --env --add channels defaults
conda config --env --set channel_priority true
```
<br>

Or set ~/.condarc to contain: <br>

```sh
channel_priority: true
channels:
  - bioconda
  - conda-forge
  - defaults
```
