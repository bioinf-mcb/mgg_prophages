# __MGG prophages__

<div align="justify">
MGG snakemake workflow for confident prophage detection. To run workflow customize config file (config.yml), download dependencies (SetupTools) and execute workflow (PROPHAGES).
</div> <br>

## __Steps__ (Linux)

### Clone repo & install snakemake

```sh
git clone https://github.com/bioinf-mcb/mgg_prophages
conda install -c conda-forge -c bioconda snakemake mamba biopython=1.79 pathlib=1.0.1
```

### Setup MGG prophage detection

<div align="justify">
Download dependencies of prophage workflow by executing SetupTools to specified directory in config file.
Execute prophage detection by running PROPHAGES snakefile.
</div> <br>

**1. Customize config file (config.yml)**

Config file description:
* Directory to download and install dependencies (SETUP_TOOLS_DIR)
* Params for NCBI workflow (INPHARED).
* Params for prophage detection workflow (PROPHAGES).

**2. SetupTools**

<div align="justify">
Executing command below will guide you through seting up dependencies. If you already did that it is enough to point in config file (SETUP_TOOLS_DIR) to directory with dependencies.
</div> <br>

```sh
snakemake --use-conda --cores all --snakefile SetupTools
```

**3. Execute prophage workflow**

<div align="justify">
Execute workflow.
</div> <br>

```sh
snakemake --use-conda --cores all --snakefile PROPHAGES
```


## __Details__ (preliminary)

PROPHAGES - detects prophages in bacterial genomes (VirSorter and PhiSpy curated with checkV)<br>
