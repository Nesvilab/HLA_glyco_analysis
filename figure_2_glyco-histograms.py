# plot sample types in the dataset
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

import plot_params
from general_params_glycFDR import out, tables

ann = pd.read_csv(tables['annotation'], sep='\t', header=0)
ann = ann[ann['HLA Class'] == 2]
dtypes = {
    'Study': str,
    'Spectrum': str,
    'Peptide': str,
    'Glycan Score': np.float64,
    'Entry Name': str
}
psm = pd.read_csv(tables['psm'],
                  sep='\t',
                  header=0,
                  low_memory=False,
                  usecols=dtypes.keys(),
                  dtype=dtypes)
pep = pd.read_csv(tables['peptide'], sep='\t', header=0)
# pep = pep[~pep['Mapped Proteins'].str.contains('rev_', na=False)]

if any(psm['Study'].str.contains('_glycFDR')):
    ann['Study'] = ann['Study'] + '_glycFDR'

psm['File'] = psm.Spectrum.str.rsplit('.', 3, expand=True)[0] + '.mzML'
# add glycosylation flag
psm['Glycosylated'] = ~psm['Glycan Score'].isnull()
# remove decoys
psm = psm[~psm['Entry Name'].str.contains('rev_')]
pep = pep[~pep['Entry Name'].str.contains('rev_')]
# retreive glycosylation flag in the peptide dataframe
pep = pep.merge(psm[['Study', 'Peptide', 'Glycosylated']].drop_duplicates())
# add annotation to psm dataframe
psm = psm.merge(
    ann[['File', 'Sample Name', 'Disease', 'Tissue Location', 'Study']])
psm_glyco = psm[~psm['Glycan Score'].isnull()].copy(deep=True)

# create dataframe for plot 1
# histogram of psm glycosylation percentage
df_plot_1 = psm.groupby([
    'Disease', 'Sample Name'
]).apply(lambda x: x['Glycan Score'].isnull().value_counts(normalize=True))

df_plot_1 = df_plot_1.unstack(level=2).fillna(0)
df_plot_1.columns = df_plot_1.columns.map({False: 'Glyco PSM', True: 'PSM'})
df_plot_1.reset_index(inplace=True)

# relove peptides from psm table that didn't make it to the peptide.tsv table
psm_pep = psm.groupby('Study').apply(
    lambda x: x[x.Peptide.isin(pep[pep.Study == x.name].Peptide.unique())])

# create dataframe for plot 2
# histograme of peptide glycosylation percentage
df_plot_2 = psm_pep.groupby(
    ['Disease', 'Sample Name']).apply(lambda x: pd.DataFrame(
        {
            'Glycosylated peptides':
            x[x['Glycosylated']]['Peptide'].nunique(),
            'Glycosylated peptides percentage':
            x[x['Glycosylated']]['Peptide'].nunique() / x['Peptide'].nunique()
        },
        index=[0])).reset_index(level=2, drop=True).reset_index()

# facet per disease
# commented since it doesn't give any additional information
# g = sb.FacetGrid(df_plot_1,
#                  col="Disease",
#                  margin_titles=True,
#                  height=3,
#                  sharey=False,
#                  sharex=False,
#                  aspect=3)
# g.map_dataframe(sb.histplot,
#                 x="Glyco PSM",
#                 binrange=(0, 0.1),
#                 binwidth=0.01,
#                 edgecolor='white',
#                 linewidth=1.2)
# for ax in g.axes.flatten():
#     vals = ax.get_xticks()
#     ax.set_xticklabels(['{:,.0%}'.format(x) for x in vals])
#     ax.set_xlabel("Percentage of glyco PSMs")

# glycosylation percentage (PSM and peptide levels) panel
fig, axes = plt.subplots(1, 2, figsize=(3.5, 1.35), constrained_layout=True)
axes = axes.flatten()
ax = axes[0]
df_plot_1['Glyco PSM'].plot(kind='hist',
                            ax=ax,
                            edgecolor='white',
                            linewidth=0.5,
                            range=(0, 0.1))
vals = ax.get_xticks()
ax.set_xticklabels(['{:,.0%}'.format(x) for x in vals])
ax.set_xlabel("Glycosylated PSMs")
ax.set_ylabel("Frequency of samples")
vals = ax.get_yticks()
ax.set_ylim([-max(vals) * 0.01, max(vals)])

ax = axes[1]
df_plot_2['Glycosylated peptides percentage'].plot(kind='hist',
                                                   ax=ax,
                                                   edgecolor='white',
                                                   linewidth=0.5,
                                                   range=(0, 0.1))
vals = ax.get_xticks()
ax.set_xticklabels(['{:,.0%}'.format(x) for x in vals])
ax.set_xlabel("Glycosylated peptides")
ax.set_ylabel("Frequency of samples")
vals = ax.get_yticks()
ax.set_ylim([-max(vals) * 0.01, max(vals)])

plt.savefig(os.path.join(out, 'plots/figure_2/panel_glyco-histograms.pdf'),
            transparent=True)
