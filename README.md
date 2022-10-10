# __MGG prophages__

<div align="justify">
MGG snakemake workflow for confident prophage detection in annotated bacterial genomes.
</div> <br>

## __Steps__ (Linux)

### Clone repo & install snakemake

```sh
git clone https://github.com/bioinf-mcb/mgg_prophages
conda install -c conda-forge -c bioconda snakemake mamba biopython=1.79 pathlib=1.0.1
```

### Run prophage detection

<div align="justify">
To run workflow customize config file (config.yml), download dependencies (SetupTools) and execute workflow (PROPHAGES).
</div> <br>


**1. Customize config file (config.yml)**

1. Directory to download and install dependencies (SETUP_TOOLS_DIR).
2. File with annotated bacterial genome(s) (recommended: via PATRIC server) (GENBANK_FILE). <br>

**Make sure that phage and protein ids are unique!** <br> <br>

**2. Download dependencies (SetupTools)**

<div align="justify">
Executing command below will guide you through setting up dependencies. If you already did that it is enough to point in config file (SETUP_TOOLS_DIR) to directory with dependencies.
</div> <br> <br>

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


## __Details__ (preliminary)

PROPHAGES - detects prophages in bacterial genomes (VirSorter and PhiSpy curated with checkV)<br>
