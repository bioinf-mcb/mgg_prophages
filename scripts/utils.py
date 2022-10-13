from Bio import SeqIO
from pathlib import Path
import pandas as pd
import itertools


# Define colors for printing.
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def choose_databases(SETUP_TOOLS_DIR):
    """ SetupTools:
    Get targets (paths) to snakemake workflow depending on user input. """

    databases2download = []
    prophages = input(f'Setup prophage detection in {SETUP_TOOLS_DIR} (PROPHAGES) [21G]? [Y/N]: ')
    print('')

    # setup prophage detection
    if prophages == 'Y':
        print(f'{bcolors.OKGREEN}Setting up prophage detection!{bcolors.ENDC}')
        pvogs = Path(SETUP_TOOLS_DIR, 'pVOGs.hmm')
        checkvdb = Path(SETUP_TOOLS_DIR, 'checkvDB')
        virsorter = Path(SETUP_TOOLS_DIR, 'virsorter_tool')

        databases2download.append(pvogs)
        databases2download.append(checkvdb)
        databases2download.append(virsorter)

    elif prophages == 'N':
        pass
    else:
        print(f'Choose Y or N. Wrong option {prophages}! Rerun!')
        exit()


    ### Check what to download and what has been already done
    print(f'\n{bcolors.OKCYAN}TO DO:{bcolors.ENDC}')
    databases2download_filtered = []
    for db in databases2download:
        if Path(db).exists():
            print(f'{bcolors.WARNING}- {str(db)} (done){bcolors.ENDC}')
        else:
            print(f'{bcolors.OKGREEN}- {str(db)} (download){bcolors.ENDC}')
            databases2download_filtered.append(db)
    print('\n')
    return databases2download_filtered


def process_input(GENBANK_FILE, OUTPUT_DIR):
    """ PROPHAGES:
    Preprocess bacterial genomes from input file: get metadata & fasta file.
    """

    # paths
    records = list(SeqIO.parse(GENBANK_FILE, 'genbank'))

    genbank = Path(OUTPUT_DIR, 'bacteria.gb')
    fasta = Path(OUTPUT_DIR, 'bacteria.fasta')
    metadata = Path(OUTPUT_DIR, 'bacteria.tsv')

    # create folder
    Path(OUTPUT_DIR).mkdir(exist_ok=True, parents=True)

    if Path(genbank).exists() or Path(fasta).exists() or Path(metadata).exists():
        print(f'{bcolors.WARNING}Seems that preprocessing was already runned! Continue with snakemake only! {bcolors.ENDC}', end='')
        return fasta, genbank, metadata

    # save fasta & gebank
    n = SeqIO.write(records, fasta, 'fasta')
    n = SeqIO.write(records, genbank, 'genbank')

    # get metadata & records
    recordIDs, recordDESCs, recordLEN = [], [], []
    for record in records:
        recordIDs.append(record.id)
        recordDESCs.append(record.description)
        recordLEN.append(len(record.seq))

    # save metadata & records
    metadata_df = pd.DataFrame({'contigID': recordIDs,
                                'contigDESC': recordDESCs,
                                'contigLEN': recordLEN})

    metadata_df.to_csv(metadata, sep='\t', index=False)

    return fasta, genbank, metadata
