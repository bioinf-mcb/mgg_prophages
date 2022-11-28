""" Draw plots from mgg_prophages output table (prophages.tsv) """

from pathlib import Path
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import string

work_dir = '/home/MCB/jkoszucki/storage/dbmgg/databases/bacteria/KASPAH/PROPHAGES_2022-11-25'
prophages_table = Path(work_dir, 'prophages.tsv')
output_dir = Path(work_dir, 'plots')

# create output folder
output_dir.mkdir(exist_ok=True, parents=True)

# load table
df = pd.read_csv(prophages_table, sep='\t')

# calculate length
df['length [kb]'] = df['seq'].str.len() / 1000

## prophages lenght per confidence
nconfidence_df = df.groupby('confidence').size()
nhigh, nmedium, nlow, nundetermined = nconfidence_df['high'], nconfidence_df['medium'], nconfidence_df['low'], nconfidence_df['undetermined']

sns.histplot(data=df, x='length [kb]', y='confidence', hue='confidence', legend=False)
plt.title(f"Histogram of prophages' length by confidence.\n nhigh={nhigh}  nmedium={nmedium}  nlow={nlow}  nundetermined={nundetermined} ")
plt.savefig(Path(output_dir, 'length_per_confidence.png'))

### prophages lenght per contig size
fig = plt.figure()
df['contig_size'] = df.apply(lambda row: row['contigID'].split('_')[-1].strip(string.digits), axis=1)
ncontig_size_df = df.groupby('contig_size').size()
nxl, nl, nm, ns, nxs = ncontig_size_df['XL'], ncontig_size_df['L'], ncontig_size_df['M'], ncontig_size_df['S'], ncontig_size_df['XS']

sns.histplot(data=df, x='length [kb]', y='contig_size', hue='contig_size', legend=False)
plt.title(f"Histogram of prophages' length by contig size category.\n nxl={nxl}  nl={nl}  nm={nm}  ns={ns}   nxs={nxs} ")
plt.savefig(Path(output_dir, 'length_per_contigsize.png'))
