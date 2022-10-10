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
    """ Get targets (paths) to snakemake workflow depending on user input. """

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


def process_input(INPUT_DIR, OUTPUT_DIR, EXTENSION):
    """
    Preprocess bacterial genomes from input directory.
    Get one fast file with curated record identifiers.
    Save metadata table with processed data.
    """

    username, userpasswors = '', ''

    fasta, metadata_df = load_bacteria_records(INPUT_DIR, OUTPUT_DIR, EXTENSION)
    genbank = submit2PATRIC(fasta, username, userpasswors)

    return fasta, genbank, metadata_df


def load_bacteria_records(INPUT_DIR, OUTPUT_DIR, EXTENSION):
    """
    Load fasta files from input directory to a list of records and curate IDs.
    """

    # warning
    print('Add your own curated ID & N50 calculation!')

    # paths
    paths = list(Path(INPUT_DIR).glob(f'*.{EXTENSION}'))
    fasta = Path(OUTPUT_DIR, 'bacteria.fasta')
    metadata_table = Path(OUTPUT_DIR, 'bacteria.tsv')

    # create folder
    Path(OUTPUT_DIR).mkdir(exist_ok=True, parents=True)

    # get metadata & records
    records_all = []
    recordIDs, recordDESCs, ncontigs, paths2frame = [], [], [], []
    for path in paths:
        records = list(SeqIO.parse(path, 'fasta'))
        for record in records:
            recordIDs.append(record.id)
            recordDESCs.append(record.description)
            ncontigs.append(len(records))
            paths2frame.append(path)
            records_all.append(record)

    # save metadata & records
    metadata_df = pd.DataFrame({'bacteriumID': recordIDs,
                                'bacteriumDESC': recordDESCs,
                                'ncontigs': ncontigs,
                                'path': paths2frame})

    n = SeqIO.write(records_all, fasta, 'fasta')
    metadata_df.to_csv(metadata_table, sep='\t', index=False)

    return fasta, metadata_df


def submit2PATRIC(fasta, username, userpasswors):
    """ ... """
    return 'test/bacterium.gb'
