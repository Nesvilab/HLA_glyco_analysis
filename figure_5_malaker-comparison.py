import os
import re
from glob import glob

import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib_venn import venn2

# for visual settings and compatibility with adobe illustrator
import plot_params
from general_params_glycFDR import database, out

tsv_path_1 = os.path.join(out, 'tables/psm_glyco.tsv')
tsv_path_3 = os.path.join(out, 'protein_glycosylation//Malaker2017.tsv')

tsv_files = glob(os.path.join(out, 'protein_glycosylation/*.csv'))

psm_glyco = pd.read_csv(tsv_path_1, sep='\t', header=0, low_memory=False)
psm_glyco = psm_glyco[psm_glyco['Glycan q-value'] <= 0.05]

# comparison to Malaker et al. 2017
malaker = pd.read_csv(tsv_path_3, sep='\t', header=0, low_memory=False)
malaker['Start'] = malaker.Peptide.str.extract("^([0-9]+)")[0]
malaker['End'] = malaker.Peptide.str.extract("([0-9]+)$")[0]
malaker['Peptide'] = malaker.Peptide.str.extract("([A-Z]+)",
                                                 flags=re.IGNORECASE)[0]

malaker['Position in peptide'] = malaker['Peptide'].map(
    lambda x: [i for i, c in enumerate(x) if c == 'n'])

malaker = malaker.explode('Position in peptide')
oglyco = malaker['Position in peptide'].isnull()
malaker['Glycosylation type'] = 'N-linked'
malaker.loc[oglyco, 'Glycosylation type'] = 'O-linked'

# >>> malaker['Glycosylation type'].value_counts()
# N-linked    84
# O-linked     6

malaker.dropna(subset=['Position in peptide'], inplace=True)
malaker['Position in protein'] = malaker['Start'].astype(
    int) + malaker['Position in peptide'].astype(int)

db = {}
for record in SeqIO.parse(database, format='fasta'):
    if record.id.startswith('rev_'):
        continue
    key = record.description.split(' ', 1)[1].split(' OS=')[0]
    value = record.id.split('|')[-1]
    db[key] = value

db['HLA class II DP alpha 1 chain'] = 'DPA1_HUMAN'
db['Interleukin 6 receptor subunit beta'] = 'IL6RB_HUMAN'
malaker['Entry Name'] = malaker['Source Protein'].map(db)
malaker['Glyco site'] = malaker['Entry Name'] + '_' + malaker[
    'Position in protein'].map(str)

psm_glyco['inMalaker'] = psm_glyco['Glyco site'].isin(
    malaker['Glyco site'].unique())

malaker['inMSFragger'] = malaker['Glyco site'].isin(
    psm_glyco['Glyco site'].unique())

set1 = set(psm_glyco['Glyco site'])

set2 = set(malaker['Glyco site'])

fig, ax = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
venn2([set1, set2],
      set_labels=('MSFragger-Glyco', 'Malaker et al. 2017'),
      ax=ax)
plt.savefig(os.path.join(out, 'plots/figure_2/panel_malaker-comparison.pdf'),
            transparent=True)

malaker.to_csv(os.path.join(out, 'tables/Malaker_comparison.tsv'),
               sep='\t',
               header=True,
               index=False)
