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
def get_decontaminated_prophage(row):
    """
    Get new viral location (decontaminated) from processed checkv contamination.tsv table.
    This function assumes that in region_types columns string viral occurs only once!
    But maybe on large scalen analysis it can get more complicated like: host,viral,host,viral
    """

    region_types = row['region_types'].split(',')
    region_coords_bp = row['region_coords_bp'].split(',')

    for region, coords in zip(region_types, region_coords_bp):
        if region == 'viral':
            start, end = coords.split('-')
            return pd.Series([start, end])


def get_prophageID(n_prophageIDs, prefix=''):
    """Generate just simple prophage IDs: suffix_PHAGEnumber"""

    format = len(str(n_prophageIDs))
    prophageIDs = [f'{prefix}_PHAGE' + f'{num}'.zfill(format) for num in list(range(1,n_prophageIDs+1))]
    return prophageIDs


def extract_phages(row, records):
    """ Extract fragment of bacterial genome
    records: bacterial records
    row: phage metadata from dataframe
    """

    primary_prophageID, start, end = row['primary_prophageID'], int(row['start']), int(row['end'])

    for record in records:
        if record.id.replace('.', '_') == primary_prophageID: break # remove dots ...

    seq = Seq(record.seq[start-1:end+1])
    return pd.Series([seq])


def get_records(row):
    """ Convert dataframe to records. """

    contigID, start, end, seq = row['contigID'], row['start'], row['end'], row['seq']
    record = SeqRecord(seq=Seq(seq), id=f'{contigID}_{start}_{end}', description=f'primary_prophage_length={len(seq)}')
    return record



# paths & params
quality = Path(snakemake.input.checkv_dir, 'quality_summary.tsv')
contamination = Path(snakemake.input.checkv_dir, 'contamination.tsv')
union_prophages = Path(snakemake.input.union_prophages)

prophages_fasta = snakemake.output.fasta
prophages_tsv = snakemake.output.tsv
PREFIX = snakemake.params.PREFIX

# load tables
quality_df = pd.read_csv(quality, sep='\t')
contamination_df = pd.read_csv(contamination, sep='\t')


# filter & rename columns
quality_cols = ['contig_id', 'provirus', 'completeness', 'completeness_method']
contamination_cols = ['contig_id', 'provirus', 'region_types', 'region_coords_bp']
rename_map = {'contig_id': 'primary_prophageID',
              'completeness_method': 'completeness_method_and_confidence'}

quality_df = quality_df[quality_cols]
contamination_df = contamination_df[contamination_cols]

quality_df.rename(rename_map, inplace=True, axis=1)
contamination_df.rename(rename_map, inplace=True, axis=1)

### format columns
# confidence
quality_df['confidence'] = quality_df['completeness_method_and_confidence'].str.split(' ', expand=True)[1]
quality_df['confidence'].fillna(False, inplace=True)

filt_high = (quality_df['confidence'] == '(high-confidence)')
filt_medium = (quality_df['confidence'] == '(medium-confidence)')
filt_low = (quality_df['confidence'] == '(low-confidence)')
filt_undetermined = (quality_df['confidence'] == False)

quality_df.loc[filt_high, 'confidence'] = 'high'
quality_df.loc[filt_medium, 'confidence'] = 'medium'
quality_df.loc[filt_low, 'confidence'] = 'low'
quality_df.loc[filt_undetermined, 'confidence'] = 'undetermined'

# completeness
quality_df['completeness'].fillna(0, inplace=True)
quality_df['completeness'] = pd.to_numeric(quality_df['completeness'], downcast='integer')

### get final table
filt_decontaminate = (contamination_df['provirus'] == 'Yes')
filt_clean = (contamination_df['provirus'] == 'No')

clean_df = contamination_df.loc[filt_clean].copy() # clean prophages (no contamination)
decontaminate_df = contamination_df.loc[filt_decontaminate].copy() # contaminated prophages

# remove contamination
# print('get_decontaminated_prophage function assumes that viral string occurs only once in region_types (test on bigger dataset)!!!!')
decontaminate_df[['start', 'end']] = decontaminate_df.apply(get_decontaminated_prophage, axis=1)

# get clean prophages location
if len(clean_df):
    clean_df[['start', 'end']] = clean_df.apply(lambda row: pd.Series(row['primary_prophageID'].split('_')[-2:]), axis=1)
    # add clean phages
    decontaminate_df = pd.concat([decontaminate_df, clean_df])

decontaminate_df['contigID'] = decontaminate_df.apply(lambda row: '_'.join(row['primary_prophageID'].split('_')[:-2]), axis=1)
decontaminate_df.sort_values(['contigID', 'start'], inplace=True, ascending=[False, True])

### give prophage IDs
# generate IDs
n_prophageIDs = len(decontaminate_df)
prophageIDs = get_prophageID(n_prophageIDs, prefix=PREFIX)
decontaminate_df['prophageID'] = prophageIDs # assign IDs

### merge completeness & decontamination
checv_df = quality_df.merge(decontaminate_df, on=['primary_prophageID', 'provirus'], how='outer')
checkv_cols = ['prophageID', 'primary_prophageID', 'contigID',
               'completeness', 'confidence', 'start', 'end',
               'completeness_method_and_confidence', 'region_types', 'region_coords_bp']

checv_df = checv_df[checkv_cols]
checv_df.sort_values(['contigID', 'start'], inplace=True, ascending=[False, True])
checv_df.reset_index(drop=True, inplace=True)
checv_df.fillna('None', inplace=True)


### extract decontaminated sequences
union_prophages_records = list(SeqIO.parse(union_prophages, 'fasta'))
checv_df['seq'] = checv_df.apply(extract_phages, args=[union_prophages_records], axis=1) # extract seq
checv_df.to_csv(prophages_tsv, sep='\t', index=False)


# convert to fasta
prophages_records = checv_df.apply(get_records, axis=1) # convert to records
prophages_records = prophages_records.to_list() # convert series2list
n = SeqIO.write(prophages_records, prophages_fasta, 'fasta') # save fasta
