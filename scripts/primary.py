"""
Get union of detections from VirSorter and PhiSpy.
Extend detections by extend variable from config file (default 2kb).

Extract prophages to fasta files and genbank files (input bacterial annotation).

Putative VirSorter prophages:
categories: 1,2 (whole contig)
categories: 4,5 (within contig)

Putative PhiSpy prophages: all
"""

# modules
import re
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# functions
def collapse_overlapping(input_detections):
    """ Collapse overlapping detections.

    INPUT: list of tuples/lists with starts and ends
    OUTPUT: list of tuples/lists with starts and ends, but merged if were overlapping.

    eg,
    INPUT: [[10, 20], [18, 30], [100, 1000], [120, 800], [120, 1001], [0, 8]]
    OUPUT: [(0, 8), (10, 30), (100, 1001)]

    """

    # covert to sets
    detections = []
    for d in input_detections:
        start, end = int(d[0]), int(d[1])
        if start <= end:
            detections.append(set(list(range(start, end+1))))
        else:
            detections.append(set(list(range(end, start+1))))

    # identify overlapping & collapse
    merged = True
    while merged:
        merged = False
        results = []
        while detections:
            common, rest = detections[0], detections[1:]
            detections = []
            for x in rest:
                if x.isdisjoint(common):
                    detections.append(x)
                else:
                    merged = True
                    common = common.union(x)

            results.append(common)
        detections = results

    # get locations from sets
    detections = [(min(d), max(d)) for d in detections]
    detections.sort(key=lambda d: d[0]) # sort

    return detections


def extract_phages(row, records):
    """ Extract fragment of bacterial genome
    records: bacterial records
    row: phage metadata from dataframe
    """

    contigID, start, end = row['contigID'], row['start'], row['end']

    for record in records:
        print(record.id, contigID, start, end)
        if record.id == contigID: break

    seq = Seq(record.seq[start:end+1])

    return pd.Series([seq])


def get_records(row):
    """ Convert dataframe to records. """

    contigID, start, end, seq = row['contigID'], row['start'], row['end'], row['seq']
    record = SeqRecord(seq=Seq(seq), id=f'{contigID}_{start}_{end}', description=f'primary_prophage_length={len(seq)}')
    return record


# paths & params
phispy_table = Path(snakemake.input.phispy)
virsorter_dir = Path(snakemake.input.virsorter, 'Predicted_viral_sequences')
metadata_table = Path(snakemake.input.metadata)
genbank = Path(snakemake.input.genbank)

fasta_output = snakemake.output.fasta  # prophages fasta
union_output = Path(snakemake.output.union)
primary_output = Path(snakemake.output.primary)
phispy_output = Path(snakemake.output.phispy)
virsorter_output = Path(snakemake.output.virsorter)

PRIMARY_EXTEND = snakemake.params.PRIMARY_EXTEND
log = Path(str(snakemake.log))


# #### testing
# phispy_table = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/raw/phispy.tsv')
# virsorter_dir = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/raw/virsorter', 'Predicted_viral_sequences')
# metadata_table = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/0_input/bacteria.tsv')
# genbank = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/0_input/bacteria.gb')
#
# fasta_output = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/tmp/prophages.fasta')  # prophages fasta
# union_output = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/tmp/union.tsv')
# primary_output = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/tmp/primary.tsv')
# phispy_output = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/tmp/phispy.tsv')
# virsorter_output = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/tmp/virsorter.tsv')
# virsorter_raw = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/tmp/virsorter_raw.tsv')
#
# PRIMARY_EXTEND = 2000
# log = Path('/home/MCB/jkoszucki/phagedb/PROPHAGES_2022-10-14/1_primary/tmp/log')




                        ##############################
                        ####### PROCESS PHISPY #######
                        ##############################

### process phispy results
phispy_df = pd.read_csv(phispy_table, sep='\t', header=None)
phispy_df = phispy_df.iloc[:, 1:4] # only important columns
phispy_df.columns = ['contigID', 'start', 'end']

phispy_df['contigID'] = phispy_df.apply(lambda row: row['contigID'].replace('.', '_'), axis=1) # curate dot in IDs :D irony
phispy_df['tool'] = 'phispy'
phispy_df.sort_values('start', ascending=True, inplace=True)
phispy_df.to_csv(phispy_output, sep='\t', index=False)

                        #################################
                        ####### PROCESS VIRSORTER #######
                        #################################

### process virsorter results
# load fasta files
virsorter_fasta = virsorter_dir.glob('*[1254].fasta') # only specified categories

headers = []
for fasta in virsorter_fasta:
    records = SeqIO.parse(fasta, 'fasta')
    for record in records:
        headers.append(record.id) # extract headers


### extract localisation & contigID
contigIDs, starts, ends, circurality = [], [], [], []
for header in headers:
    # phage within contig
    try:
        localisation = re.search('\d{1,12}-\d{1,12}-cat_[45]', header).group() # localisation
        start, end = localisation.split('-')[:2]
        start = int(start) + 1

        contigID = re.search('VIRSorter_.*_gene', header).group() # contigID
        contigID = contigID.strip('VIRSorter_')
        contigID_split = [e for e in contigID.split('_')[:-3] if e]
        contigID = '_'.join(contigID_split)

    # phage as a whole contig
    except AttributeError:
        try:
            # check if contig is circular
            if '-circular' in header:
                circular = True
                header = header.replace('-circular', '')
            else: circular = False

            contigID = re.search('VIRSorter_.*-cat', header).group() # contigID
            contigID = contigID.strip('VIRSorter_').strip('-cat')
            start, end = 1, 'contig_len'

        # something went wrong
        except:
            contigID, start, end = 'NotFoundID', 0, 0

    contigIDs.append(contigID)
    starts.append(start)
    ends.append(end)
    circurality.append(circular)

# virsorter table
virsorter_df = pd.DataFrame({'contigID': contigIDs,
                             'start': starts,
                             'end': ends,
                             'circular': circurality,
                             'header': headers})

virsorter_df['tool'] = 'virsorter'
virsorter_df.sort_values('start', ascending=True, inplace=True)
virsorter_df.to_csv(virsorter_output, sep='\t', index=False)


                        ########################################
                        ####### GET FINAL DETECTIONS ###########
                        ########################################

### primary detections
# load tables
primary_df = pd.concat([phispy_df, virsorter_df]) # primary detections
metadata_df = pd.read_csv(metadata_table, sep='\t') # bacterial metadata

# omg, dots...
metadata_df['contigID'] = metadata_df.apply(lambda row: row['contigID'].replace('.', '_'), axis=1) # make pretty
primary_df.reset_index(drop=True, inplace=True) # make pretty

# add contig info
primary_df = primary_df.merge(metadata_df, on='contigID', how='left')

# get end for whole contig prophages
filt = (primary_df['end'] == 'contig_len')
primary_df.loc[filt, 'end'] = primary_df.loc[filt, 'contigLEN']

# sort & save primary detections
primary_df.sort_values(['contigLEN', 'start'], ascending=[False, True], inplace=True)
primary_df.fillna('None', inplace=True)
primary_df.to_csv(primary_output, sep='\t', index=False)

### union
# union of detections (collapse overlapping ones) per bacterial contig
groups = primary_df.groupby('contigID')
unions, contigIDs = [], []
new_start, new_end = [], []
for contigID, group in groups:
    starts, ends = group['start'].to_list(), group['end'].to_list()
    locations = zip(starts, ends)

    unions_per_contig = collapse_overlapping(locations)

    for start, end in unions_per_contig:
        contigIDs.append(contigID)
        new_start.append(start)
        new_end.append(end)

union_df = pd.DataFrame({'contigID': contigIDs,
                         'start_union': new_start,
                         'end_union': new_end})

union_df = union_df.merge(metadata_df, on='contigID', how='left')
# union_df = union_df.merge(primary_df, on='contigID', how='left', suffixes=('_union', '_primary'))


### extend
# extend union of detections
union_df['start'] = (union_df['start_union'] - PRIMARY_EXTEND)
union_df['end'] = (union_df['end_union'] + PRIMARY_EXTEND)

filt_start = (union_df['start'] <= 0)
filt_end = (union_df['end'] >= union_df['contigLEN'])

union_df.loc[filt_start, 'start'] = 1
union_df.loc[filt_end, 'end'] = union_df['contigLEN']

# union_df.to_csv(primary_output, sep='\t', index=False) # overwrite primary (nicer primary table)


                        ###################################
                        ####### EXTRACT FASTA  ############
                        ###################################

bacterial_records = list(SeqIO.parse(genbank, 'genbank')) # load bacterial genomes
union_df.drop_duplicates(subset=['contigID', 'start', 'end'], inplace=True)
union_df['seq'] = union_df.apply(extract_phages, args=[bacterial_records], axis=1) # extract seq

cols = ['contigID', 'start', 'end', 'start_union', 'end_union', \
        'contigDESC', 'contigLEN', 'seq']

union_df.sort_values(['contigLEN', 'start'], ascending=[False, True], inplace=True)
union_df[cols].to_csv(union_output, sep='\t', index=False) # save table

# convert to fasta
prophage_records = union_df.apply(get_records, axis=1) # convert to records
prophage_records = prophage_records.to_list() # convert series2list
n = SeqIO.write(prophage_records, fasta_output, 'fasta') # save fasta

# get & save log
log_message = []
log_message.append(f'Prophages in primary predictions: {primary_df.shape[0]}')
log_message.append(f'Prophages collapsed (union): {union_df.shape[0]}')
log_message.append(f'Prophages saved to fasta: {n}.')

with open(log, 'w+') as l:
    l.write('\n'.join(log_message) + '\n')
