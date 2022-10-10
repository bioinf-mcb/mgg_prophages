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
# phispy_table = Path(snakemake.input.phispy)
# virsorter_dir = Path(snakemake.input.virsorter)
#
# phispy = Path(snakemake.output.phispy)
# virsorter = Path(snakemake.output.virsorter)
# primary = Path(snakemake.output.primary)
# extend = Path(snakemake.params.PRIMARY_EXTEND # bp to extend union of detections)

# paths
phispy_table = Path('/home/MCB/jkoszucki/PROPHAGES_2022-10-10/phispy.tsv')
virsorter_dir = Path('/home/MCB/jkoszucki/PROPHAGES_2022-10-10/virsorter/Predicted_viral_sequences')
metadata_table = Path('/home/MCB/jkoszucki/PROPHAGES_2022-10-10/bacteria.tsv')

phispy_output = Path('/home/MCB/jkoszucki/PROPHAGES_2022-10-10/tmp1.tsv')
headers_output = Path('/home/MCB/jkoszucki/PROPHAGES_2022-10-10/headers.tsv')
virsorter_output = Path('/home/MCB/jkoszucki/PROPHAGES_2022-10-10/tmp2.tsv')
primary_output = Path('/home/MCB/jkoszucki/PROPHAGES_2022-10-10/tmp3.tsv')


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
with open(headers_output, 'w+') as f:
    f.write('\n'.join(headers))

### extract localisation & contigID
contigIDs, starts, ends = [], [], []
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
            contigID = re.search('VIRSorter_.*-cat', header).group() # contigID
            contigID = contigID.strip('VIRSorter_').strip('-cat')
            start, end = 1, 'contig_len'

        # something went wrong
        except:
            contigID, start, end = 'NotFoundID', 0, 0

    contigIDs.append(contigID)
    starts.append(start)
    ends.append(end)

# virsorter table
virsorter_df = pd.DataFrame({'contigID': contigIDs,
                             'start': starts,
                             'end': ends})

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

# get whole virsorter contigs (end of prophage is length of the contig)
filt_contigs = (primary_df['end'] == 'contig_len')
contigs_df = primary_df.loc[filt_contigs]
contigs_df.merge(metadata_df, on='contigID', how='left')
contigs_df['end'] = contigs_df['contigLEN']

# combine whole contigs & rest
rest_df = primary_df.loc[~filt_contigs]
rest_df.merge(metadata_df, on='contigID', how='left')
primary_df = pd.concat(contigs_df, rest_df)

primary_df.to_csv(primary_output, sep='\t', index=False)

# primary_df.groupby('contigID')
