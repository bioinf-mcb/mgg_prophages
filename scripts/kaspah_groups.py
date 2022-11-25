from Bio import SeqIO
from pathlib import Path
import pandas as pd
import itertools

work_dir = '/home/MCB/jkoszucki/storage/dbmgg/databases/bacteria/KASPAH'
kaspah_groups_table = '/home/MCB/jkoszucki/Code/mgg_bacteria/dependencies/tables/kaspah_coding.tsv'

genbank_dir = Path(work_dir, 'ONE_FILE_VRS_BACTERIA_2022-11-05/1_PATRIC_batches')
metadata_tsv = Path(work_dir, 'ONE_FILE_VRS_BACTERIA_2022-11-05/bacteria.tsv')

genbank_out = Path(work_dir, 'BACTERIA_2022-11-25/bacteria.gb')
fasta_out = Path(work_dir, 'BACTERIA_2022-11-25/bacteria.fasta')
metadata_out = Path(work_dir, 'BACTERIA_2022-11-25/bacteria.tsv')

# load tables
df = pd.read_csv(kaspah_groups_table, sep='\t')
metadata_df = pd.read_csv(metadata_tsv, sep='\t')

genbanks = list(Path(genbank_dir).glob('*.gb'))

records = [list(SeqIO.parse(genbank, 'genbank')) for genbank in genbanks]
records = list(itertools.chain(*records))

# rename contigs to new group
groups_mapper = {}
for record in records:
    contigID = record.id
    old_group = contigID.split('_')[-2]
    new_group = df.loc[df['groupID'] == old_group, 'new_groupID'].iloc[0]

    record.id = record.id.replace(old_group, new_group)
    record.name = ''
    record.description = ''

    groups_mapper[contigID] = record.id.replace(old_group, new_group)

# rename contigs in metadata
metadata_df['contigID'] = metadata_df['contigID'].map(groups_mapper)


SeqIO.write(records, genbank_out, 'genbank')
SeqIO.write(records, fasta_out, 'fasta')
metadata_df.to_csv(metadata_out, sep='\t', index=False)
