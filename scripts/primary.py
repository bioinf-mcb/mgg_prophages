"""
Get union of detections from VirSorter and PhiSpy.
Extend detections by extend variable (default 2kb).

Putative VirSorter prophages:
categories: 1,2 (whole contig)
categories: 4,5 (within contig)

Putative PhiSpy prophages: all
"""

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


# modules
import re
import pandas as pd
from Bio import SeqIO
from pathlib import Path

# paths & params
phispy = Path(snakemake.input.phispy)
virsorter = Path(snakemake.input.virsorter)
metadata = Path(snakemake.input.metadata)

union_output = Path(snakemake.output.union)
primary_output = Path(snakemake.output.primary)
phispy_output = Path(snakemake.output.phispy)
virsorter_output = Path(snakemake.output.virsorter, 'Predicted_viral_sequences')
virsorter_raw = Path(snakemake.output.virsorter_raw)

# # paths
# phispy_table = Path('/home/MCB/jkoszucki/phagedb/others/PROPHAGES_2022-10-13/1_primary/raw/phispy.tsv')
# virsorter_dir = Path('/home/MCB/jkoszucki/phagedb/others/PROPHAGES_2022-10-13/1_primary/raw/virsorter/Predicted_viral_sequences')
# metadata_table = Path('/home/MCB/jkoszucki/phagedb/others/PROPHAGES_2022-10-13/0_input/bacteria.tsv')
#
# phispy_output = Path('/home/MCB/jkoszucki/phagedb/others/PROPHAGES_2022-10-13/1_primary/results/1_phispy.tsv')
# virsorter_raw = Path('/home/MCB/jkoszucki/phagedb/others/PROPHAGES_2022-10-13/1_primary/results/0_raw_virsorter.tsv')
# virsorter_output = Path('/home/MCB/jkoszucki/phagedb/others/PROPHAGES_2022-10-13/1_primary/results/2_virsorter.tsv')
# primary_output = Path('/home/MCB/jkoszucki/phagedb/others/PROPHAGES_2022-10-13/1_primary/results/3_primary.tsv')
# union_output = Path('/home/MCB/jkoszucki/phagedb/others/PROPHAGES_2022-10-13/1_primary/results/4_union.tsv')


                        ##############################
                        ####### PROCESS PHISPY #######
                        ##############################

### process phispy results
phispy_df = pd.read_csv(phispy_table, sep='\t', header=None)
phispy_df = phispy_df.iloc[:, 1:4] # only important columns
phispy_df.columns = ['contigID', 'start', 'end']

phispy_df['contigID'] = phispy_df.apply(lambda row: row['contigID'].replace('.', '_'), axis=1) # curate dot in IDs :D irony
phispy_df['tool'] = 'phispy'
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

# save raw headers
with open(virsorter_raw, 'w+') as f:
    f.write('\n'.join(headers) + '\n')

### extract localisation & contigID
contigIDs, starts, ends, circurality = [], [], [], []
for header in headers:
    # phage within contig
    try:
        localisation = re.search('\d{1,12}-\d{1,12}-cat_[45]', header).group() # localisation
        print(localisation)
        start, end = localisation.split('-')[:2]

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
                             'circular': circurality})

virsorter_df['tool'] = 'virsorter'
virsorter_df.to_csv(virsorter_output, sep='\t', index=False)


                        ########################################
                        ####### GET UNION OF RESULTS ###########
                        ########################################

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

# sort
primary_df.sort_values(['contigID', 'contigLEN'], ascending=[False, False], inplace=True)
primary_df.fillna('None', inplace=True)
primary_df.to_csv(primary_output, sep='\t', index=False)


# union of detections (collapse overlapping ones) per bacterial contig
groups = primary_df.groupby('contigID')
unions, contigIDs = [], []
new_start, new_end = [], []
for contigID, group in groups:
    starts, ends = group['start'].to_list(), group['end'].to_list()
    locations = zip(starts, ends)

    unions_per_contig = collapse_overlapping(locations)

    contigIDs.append(contigID)
    for start, end in unions_per_contig:
        new_start.append(start)
        new_end.append(end)

union_df = pd.DataFrame({'contigID': contigIDs,
                             'start': new_start,
                             'end': new_end})

union_df = union_df.merge(primary_df, on='contigID', how='left', suffixes=('_union', '_primary'))
union_df.to_csv(union_output, sep='\t', index=False)
