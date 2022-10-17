# __MGG prophages__

<div align="justify">
Automatic & scalable prophage detection meta-tool. Exploits two tools for prophage detection (VirSorter & Phispy) and performs decontamination of bacterial DNA (CheckV). Primary prophages are collapsed (union) and analyzed by CheckV. Based on assigned completeness and confidence prophages can be filtered, accordingly to subsequent analyses.
</div>

## __Steps__ (Linux)

### Clone repo & install snakemake

```sh
git clone https://github.com/bioinf-mcb/mgg_prophages
conda install -c conda-forge -c bioconda snakemake mamba biopython=1.79 pathlib=1.0.1 pandas datetime
```

### Run prophage detection

<div align="justify">
To run workflow customize config file, download dependencies and execute workflow. Setup config file by providing paths to download dependencies and annotated bacterial genomes in genbank file (BV-BRC recommended) and run analysis. If dependencies were already downloaded just provide path to the folder. To run test provide path to the test.gb in config file and execute workflow (optional). <br><br>

```sh
snakemake --use-conda --cores all --snakefile SetupTools # setup dependencies
snakemake --use-conda --cores all --snakefile PROPHAGES # execute workflow
```


</div> <br>
** comment: Make sure that bacterial and protein ids are unique and avoid white characters in bacterial contig ids. **


## Results

<div align="justify">
Sequences and metadata of prophages and prophage-like elements are final output of the workflow ("prophages" files).
To obtain high-confidence prophages one should consider prophages full-filling criteria of completeness at least:  >= 90% and >= medium confidence.
</div> <br>


## __Details__ (preliminary)

<div align="justify">
Tool detects prophages in bacterial genomes by using two complementary tools: VirSorter and PhiSpy. Primary prophage detection is any region detected the tools as a prophage. Subsequently, detections are collapsed (union) decontaminated from bacterial DNA and their completeness (with confidence) is estimated (CheckV).
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
