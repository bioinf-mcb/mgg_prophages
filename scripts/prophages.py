"""
Load CheckV results on union of primary detections (collapsed overlaps).
Decontaminate by using CheckV results and assing prophage completeness & confidence (of completeness).
Saves table of detections & fasta file.
"""


# modules
from pathlib import Path
import pandas as pd
import numpy as np
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# functions
def get_confidence(row):
    """ Based on completeness_method column from checkv quality_summary.tsv
    file get confidence of completeness estimation """

    # get confidence info
    confidence_raw = row['confidence']
    if confidence_raw != 'None': confidence_raw = confidence_raw.split(' ')[-1]

    # re-define confidence
    if confidence_raw == 'None': confidence = 'undetermined'
    elif confidence_raw == '(lower-bound)': confidence = 'low'
    elif confidence_raw == '(medium-confidence)': confidence = 'medium'
    elif confidence_raw == '(high-confidence)': confidence = 'high'
    else: confidence = 'Warning! Confidence not found!'

    return pd.Series([confidence])


def split_primary_prophageID(row):
    """ split primary prophageID and to contigID, start_primary, end_primary """

    try:
        contigID = '_'.join(str(row['primary_prophageID']).split('_')[:2])
        start_primary, end_primary = str(row['primary_prophageID']).split('_')[2:]
    except ValueError: 
        print('ERROR! Remove "_" characters from bacterial genome identifiers!')
        exit()

    return pd.Series([contigID, start_primary, end_primary])


def relocate_phages(row):
    """ Based on location of prophages on bacterial contigs (primary) and
    relative location of CheckV curated prophages (relative to location of extracted primary prophages)
    define new curated location of prophages on bacterial contigs. """

    # get curated location if any
    if row['end_relative'] != 0:
        start = row['start_primary'] + row['start_relative'] - 1
        end = row['start_primary'] + row['end_relative'] - 1
    else:
        start = row['start_primary']
        end = row['end_primary']

    return pd.Series([start, end])


def extend(row, EXTEND):
    """ Extend final detection of prophages by user defined value """

    start = row['start'] - EXTEND
    end = row['end'] + EXTEND
    contig_length = row['contigLEN']

    if start <= 1: start = 1
    if end >= contig_length: end = contig_length

    return pd.Series([start, end])


def get_prophageID(n_prophageIDs, prefix=''):
    """Generate just simple prophage IDs: suffix_PHAGEnumber"""

    format = len(str(n_prophageIDs))
    prophageIDs = [f'{prefix}_PHAGE' + f'{num}'.zfill(format) for num in list(range(1,n_prophageIDs+1))]
    return prophageIDs


def add_size_category(row):
    """ get prophage size category"""

    prophageID = row['prophageID']
    start, end = int(row['start']), int(row['end'])
    length = end - start + 1

    if length <= 10000: size_category = 'S'
    elif length > 10000 and length <= 90000: size_category= 'M'
    else: size_category = 'L'

    return f'{prophageID}_{size_category}'


def extract_phages(row, bacteria_records):
    """ Extract fragment of bacterial genome
    records: bacterial records
    row: phage metadata from dataframe
    """

    contigID, start, end = row['contigID'], int(row['start']), int(row['end'])

    seq = f'Error in extracting seq ({contigID}, {start}, {end})'
    for record in bacteria_records:
        if record.id == contigID:
            seq = record.seq[start-1:end+1]

    return pd.Series([str(seq)])


def get_records(row):
    """ Convert dataframe to records. Assing phage size category"""

    prophageID, seq = row['prophageID'], row['seq']
    length = len(seq)

    record = SeqRecord(seq=Seq(seq), id=f'{prophageID}', description='')
    return record



# paths & params
# quality = Path(snakemake.input.checkv_dir, 'quality_summary.tsv')
# contamination = Path(snakemake.input.checkv_dir, 'contamination.tsv')
# bacteria = Path(snakemake.input.fasta)
# metadata = Path(snakemake.input.metadata)

# prophages_fasta = snakemake.output.fasta
# prophages_tsv = snakemake.output.tsv
# PREFIX = snakemake.params.PREFIX
# EXTEND = snakemake.params.FINAL_EXTEND

### testing
working_dir = '/net/pr2/projects/plgrid/plggmgmcb/team/jkoszucki/code/mgg_kleb/_data/2_INTERMEDIATE/08_PROPHAGES/KAG/BATCH1'
input_dir = Path(working_dir, '0_input')
checkv_dir = Path(working_dir, '2_checkv')

quality = Path(checkv_dir, 'quality_summary.tsv')
contamination = Path(checkv_dir, 'contamination.tsv')
bacteria = Path(input_dir, 'bacteria.fasta')
metadata = Path(input_dir, 'bacteria.tsv')

prophages_fasta = Path(working_dir, 'prophages.fasta')
prophages_tsv = Path(working_dir, 'prophages.tsv')
PREFIX = 'T1'
EXTEND = 2000


# quality = Path('/home/MCB/jkoszucki/storage/dbmgg/databases/bacteria/KASPAH_2022-11-08/PROPHAGES_2022-11-14/2_checkv','quality_summary.tsv')
# contamination = Path('/home/MCB/jkoszucki/storage/dbmgg/databases/bacteria/KASPAH_2022-11-08/PROPHAGES_2022-11-14/2_checkv', 'contamination.tsv')
#
# prophages_fasta = Path('/home/MCB/jkoszucki/tmp/prophages.fasta')
# prophages_tsv = Path('/home/MCB/jkoszucki/tmp/prophages.tsv')
# PREFIX = 'TEST'
# EXTEND = snakemake.params.FINAL_EXTEND

# load tables
quality_df = pd.read_csv(quality, sep='\t')
contamination_df = pd.read_csv(contamination, sep='\t')
metadata_df = pd.read_csv(metadata, sep='\t')


# filter & rename columns
quality_cols = ['contig_id', 'completeness', 'completeness_method']
contamination_cols = ['contig_id', 'region_types', 'region_coords_bp']

quality_df = quality_df[quality_cols]
contamination_df = contamination_df[contamination_cols]

rename_map = {'contig_id': 'primary_prophageID',
              'completeness_method': 'confidence'}

quality_df = quality_df.rename(rename_map, axis=1)
contamination_df = contamination_df.rename(rename_map, axis=1)

# completeness & confidence
quality_df['completeness'] = quality_df['completeness'].fillna(0)
quality_df['completeness'] = pd.to_numeric(quality_df['completeness'], downcast='integer')

quality_df['confidence'] = quality_df['confidence'].fillna('None')
quality_df['confidence'] = quality_df.apply(get_confidence, axis=1)

### main checkv table
checkv_df = quality_df.merge(contamination_df, on='primary_prophageID', how='outer')

# checkpoint
if quality_df.shape[0] == contamination_df.shape[0] & contamination_df.shape[0] == checkv_df.shape[0]: pass
else: print('WARNING! CheckV tables are inconsistent! (prophages.py)')

# relative locations of checkv viral regions to primary detections
checkv_df = checkv_df.set_index(['primary_prophageID', 'completeness', 'confidence']).apply(lambda row: row.str.split(',').explode()).reset_index()
checkv_df['region_types'] = checkv_df['region_types'].fillna('viral')
checkv_df['region_coords_bp'] = checkv_df['region_coords_bp'].fillna('0-0')

filt_viral_regions = (checkv_df['region_types'] == 'viral')
checkv_df = checkv_df.loc[filt_viral_regions]
checkv_df[['start_relative', 'end_relative']] = checkv_df['region_coords_bp'].str.split('-', expand=True)
checkv_df = checkv_df.drop(['region_types', 'region_coords_bp'], axis=1)

# primary detections start & end
checkv_df[['contigID', 'start_primary', 'end_primary']] = checkv_df.apply(split_primary_prophageID, axis=1)

# convert values to numeric
for col in ['start_primary', 'end_primary', 'start_relative', 'end_relative']:
    checkv_df[col] = pd.to_numeric(checkv_df[col], downcast='integer')

# correct prophage location based on relative checkv postions
checkv_df[['start', 'end']] = checkv_df.apply(relocate_phages, axis=1)

# extend prophages by user defined value
checkv_df = checkv_df.merge(metadata_df, on='contigID', how='left')
checkv_df[['start', 'end']] = checkv_df.apply(extend, args=[EXTEND], axis=1)
checkv_df = checkv_df.drop(['contigDESC', 'contigLEN'], axis=1)

# sort
confidence_order = ['high', 'medium', 'low', 'undetermined']
checkv_df['confidence'] = pd.Categorical(checkv_df['confidence'], categories = confidence_order)
checkv_df = checkv_df.sort_values(['confidence', 'completeness'], ascending=[True, False])

### give prophage IDs
# generate IDs
n_prophageIDs = len(checkv_df)
prophageIDs = get_prophageID(n_prophageIDs, prefix=PREFIX)
checkv_df['prophageID'] = prophageIDs # assign IDs

# add size category to prophageIDs
checkv_df['prophageID'] = checkv_df.apply(add_size_category, axis=1)

### extract decontaminated sequences
bacteria_records = list(SeqIO.parse(bacteria, 'fasta'))
checkv_df['seq'] = checkv_df.apply(extract_phages, args=[bacteria_records], axis=1) # extract seq

cols = ['prophageID', 'contigID', 'start', 'end', 'completeness', 'confidence', \
        'primary_prophageID', 'start_primary', 'end_primary', 'start_relative', 'end_relative', 'seq']

checkv_df = checkv_df[cols]
checkv_df.to_csv(prophages_tsv, sep='\t', index=False)

# convert to fasta
prophages_records = checkv_df.apply(get_records, axis=1) # convert to records
prophages_records = prophages_records.to_list() # convert series2list
n = SeqIO.write(prophages_records, prophages_fasta, 'fasta') # save fasta
