"""
Get union of detections from VirSorter and PhiSpy.
Extend detections by extend variable (default 2kb).

Putative VirSorter prophages:
categories: 1,2 (whole contig)
categories: 4,5 (within contig)

Putative PhiSpy prophages: all
"""

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


phispy_table = Path('/home/MCB/jkoszucki/Code/mgg_prophages/PROPHAGES_2022-10-10/phispy.tsv')
virsorter_dir = Path('/home/MCB/jkoszucki/Code/mgg_prophages/PROPHAGES_2022-10-10/virsorter/Predicted_viral_sequences')
phispy = Path('/home/MCB/jkoszucki/Code/mgg_prophages/PROPHAGES_2022-10-10/tmp1.tsv')
virsorter = Path('/home/MCB/jkoszucki/Code/mgg_prophages/PROPHAGES_2022-10-10/tmp2.tsv')
primary = Path('/home/MCB/jkoszucki/Code/mgg_prophages/PROPHAGES_2022-10-10/tmp3.tsv')


                        ##############################
                        ####### PROCESS PHISPY #######
                        ##############################

### process phispy results
phispy_df = pd.read_csv(phispy_table, sep='\t', header=None)
phispy_df = phispy_df.iloc[:, 1:4] # only important columns
phispy_df.columns = ['contigID', 'start', 'end']

phispy_df['contigID'] = phispy_df.apply(lambda row: row['contigID'].replace('.', '_'), axis=1) # curate dot in IDs :D irony
phispy_df['tool'] = 'phispy'
phispy_df.to_csv(phispy, sep='\t', index=False)


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
contigIDs, starts, ends = [], [], []
for header in headers:
    # phage within contig
    try:
        localisation = re.search('\d{1,9}-\d{1,9}-cat_[45]', header).group() # localisation
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
virsorter_df.to_csv(virsorter, sep='\t', index=False)


                        ########################################
                        ####### GET UNION OF RESULTS ###########
                        ########################################

pd.concat(phispy_df)
phispy_df.to_csv(primary, sep='\t', index=False)
