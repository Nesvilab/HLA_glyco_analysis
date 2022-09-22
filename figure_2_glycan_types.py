import os
from collections import defaultdict

import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sb
from matplotlib import pyplot as plt

# for visual settings and compatibility with adobe illustrator
import plot_params
from general_params_glycFDR import out

cm = 1 / 2.54  # centimeters in inches

df_bindings = pd.read_csv(os.path.join(
    out, 'tables/binding_predictions_best-allele.tsv'),
                          sep='\t',
                          header=0)

alleles = pd.read_csv(os.path.join(out, 'tables/bindings_alleles.tsv'),
                      sep='\t',
                      header=0)
complete_annotation_subject = alleles[
    alleles['Alleles annotation notes'].isnull()].Subject.tolist()
df_bindings = df_bindings[df_bindings['Binding Type'] == 'Non-trash']

df_bindings = df_bindings[df_bindings['Subject'].isin(
    complete_annotation_subject)]

# generated by './figure_2_glycosites.py'
glyco_path = os.path.join(out, 'tables/glycosites.tsv')

# reading tables for all studies
psm_glyco = pd.read_csv(glyco_path, sep='\t', header=0, low_memory=False)
psm_glyco['Glycosylated'] = True

# fetching glycan type annotation
glyc_ann_path = os.path.join(out, 'tables/NGlycans_groups.tsv')

df_glycan_annotation = pd.read_csv(glyc_ann_path,
                                   sep='\t',
                                   header=0,
                                   comment='#')
df_glycan_annotation.columns = [
    'Mass', 'Observed Modifications', 'Glycan Type'
]
df_glycan_annotation['Observed Modifications'] = df_glycan_annotation[
    'Observed Modifications'].str.replace(' %.+$', '', regex=True)

converter = {
    'O-Mannose': "Other",
    'O-glycan': "Other",
    "Truncated": "Truncated",
    "High Mannose": "High Mannose",
    "High Mannose w/Fucose": "High Mannose",
    "Complex/Hybrid w/Fucose & NeuAc": "Complex/Hybrid",
    "Truncated w/Fucose": "Truncated",
    "Truncated Complex/Hybrid": "Truncated Complex/Hybrid",
    "Truncated Complex/Hybrid w/Fucose": "Truncated Complex/Hybrid",
    "Complex/Hybrid": "Complex/Hybrid",
    "Complex/Hybrid w/Fucose": "Complex/Hybrid",
    "Complex/Hybrid w/NeuAc": "Complex/Hybrid"
}
df_glycan_annotation['Glycan Type'] = df_glycan_annotation['Glycan Type'].map(
    converter)
psm_glyco['Observed Modifications'] = psm_glyco[
    'Observed Modifications'].str.replace(' %.+$', '', regex=True)
psm_glyco = psm_glyco.merge(df_glycan_annotation, how='left')

# processing tabels and merging
df_bindings = df_bindings.merge(psm_glyco[[
    'Peptide', 'Subject', 'Glyco Position', 'Glycosylated', 'Glycan Type',
    'Mass', 'Glyco site'
]].drop_duplicates(),
                                how='left')
df_bindings['Glycosylated'] = df_bindings['Glycosylated'].fillna(False)

df_glyco = df_bindings[df_bindings['Glycosylated']].copy(deep=True)
df_glyco['Core start'] = df_glyco.apply(lambda x: x.Peptide.index(x.Core),
                                        axis=1) + 1
df_glyco['Glyco position within motif'] = df_glyco[
    'Glyco Position'] - df_glyco['Core start']
df_glyco['HLA Group'] = df_glyco['NetMHCIIpan Allele'].str.extract(
    '(D[P|Q|R])')[0]
# add inside vs outside motif column
df_glyco['Glycosylation Location Type'] = (
    df_glyco['Glyco position within motif'] >=
    0) & (df_glyco['Glyco position within motif'] <= 8)
df_glyco['Glycosylation Location Type'] = df_glyco[
    'Glycosylation Location Type'].map({
        True: 'Within HLA motif',
        False: 'Outside HLA motif'
    })
# Check glycoyslation position (inside / outside HLA motif) per HLA group
# Groups: DP, DQ, DR

df_plot_1 = df_glyco.groupby('HLA Group').apply(lambda x: x[
    'Glycosylation Location Type'].value_counts(normalize=True)).unstack(
        level=1)

df_plot_1 = round(df_plot_1, 4) * 100
cols_order = ["Within HLA motif", 'Outside HLA motif']

fig, axes = plt.subplots(1,
                         2,
                         figsize=(7, 4),
                         gridspec_kw={'width_ratios': [1, 10]},
                         constrained_layout=True)
axes = axes.flatten()

ax = axes[0]
df_plot_1[cols_order].plot(kind='bar',
                           ax=ax,
                           stacked=True,
                           color=['red', 'gray'])
ax.set_ylabel('Glyco peptides (%)')
ax.set_xlabel('HLA group')

ax = axes[1]
df_plot_2 = df_glyco.groupby('Glyco position within motif').agg({
    'Mass': 'mean'
}).unstack(level=1)
df_plot_2.plot(kind='bar', ax=ax)
ax.set_ylabel('Average glycan mass (Da)')
ax.set_xlabel('HLA group')

plt.savefig(os.path.join(
    out, 'plots/figure_2/panel_relativePositionPerHlaGroup.pdf'),
            transparent=True)

sp_width = df_glyco.groupby('HLA Group')['NetMHCIIpan Allele'].nunique()

fig, axes = plt.subplots(2,
                         3,
                         figsize=(18 * cm, 8 * cm),
                         constrained_layout=True,
                         gridspec_kw={'width_ratios': sp_width})
group_index = 0
for group, sdf in df_glyco.groupby('HLA Group'):
    # plot 1: percentage of glycan inside / outside motif
    ax = axes[0, group_index]
    df_plot_3 = sdf.groupby('NetMHCIIpan Allele').apply(lambda x: x[
        'Glycosylation Location Type'].value_counts(normalize=True)).unstack(
            level=1).fillna(0)
    df_plot_3 = round(df_plot_3, 4) * 100
    df_plot_3.sort_values(cols_order[0], inplace=True)
    cols_order = ["Within HLA motif", 'Outside HLA motif']
    y1_val = df_plot_3[cols_order[0]].values
    y2_val = df_plot_3[cols_order[1]].values
    x_val = df_plot_3.index
    ax.bar(x_val, y1_val, color='red', label='Within HLA motif')
    ax.bar(x_val,
           y2_val,
           bottom=y1_val,
           color='gray',
           label='Outside HLA motif')
    ax.set_ylabel('Glyco peptides (%)')
    ax.set_xlabel('HLA group')
    ax.set_xticklabels([])
    ax.legend()
    # plot 2: average glycan mass
    ax = axes[1, group_index]
    df_plot_4 = sdf.groupby('NetMHCIIpan Allele')['Mass'].mean()
    df_plot_4 = df_plot_4.loc[x_val]
    ax.bar(x_val, df_plot_4.values)
    x_val = x_val.str.replace('HLA-', '')
    x_val = x_val.str.replace('_', '')
    ax.set_xticklabels(x_val, rotation=90, fontsize=8)
    ax.set_ylabel('Average glycan mass (Da)')
    ax.set_xlabel('HLA group')
    min1, max1 = axes[0, group_index].get_xlim()
    min2, max2 = axes[1, group_index].get_xlim()
    axes[0, group_index].set_xlim(min(min1, min2), max(max1, max2))
    axes[1, group_index].set_xlim(min(min1, min2), max(max1, max2))
    group_index += 1

plt.savefig(os.path.join(
    out, 'plots/figure_2/panel_relativePositionPerHlaAllele.pdf'),
            transparent=True)

sp_width = df_glyco.groupby('HLA Group').apply(
    lambda x: x[['Subject', 'NetMHCIIpan Allele']].drop_duplicates().shape[0])
fig, axes = plt.subplots(2,
                         3,
                         figsize=(50 * cm, 15 * cm),
                         constrained_layout=True,
                         gridspec_kw={'width_ratios': sp_width})
first_plot = True
group_index = 0
min_sequences = 10
for group, sdf in df_glyco.groupby('HLA Group'):
    # plot 1: percentage of glycan inside / outside motif
    ax = axes[0, group_index]
    df_plot_3 = sdf.groupby([
        'Subject', 'NetMHCIIpan Allele'
    ]).apply(lambda x: x['Glycosylation Location Type'].value_counts(
        normalize=True)).unstack(level=2).fillna(0)
    df_plot_3 = round(df_plot_3, 4) * 100
    df_plot_3.sort_values(cols_order[0], inplace=True)
    df_plot_3.index = df_plot_3.index.map(lambda x: ' - '.join(x))
    cols_order = ["Within HLA motif", 'Outside HLA motif']
    y1_val = df_plot_3[cols_order[0]].values
    y2_val = df_plot_3[cols_order[1]].values
    x_val = df_plot_3.index
    ax.bar(x_val, y1_val, color='red', label='Within HLA motif')
    ax.bar(x_val,
           y2_val,
           bottom=y1_val,
           color='gray',
           label='Outside HLA motif')
    ax.set_ylabel('Glyco peptides (%)')
    ax.set_xlabel('HLA group')
    ax.set_xticklabels([])
    if first_plot:
        ax.legend()
        first_plot = False
    # plot 2: average glycan mass
    ax = axes[1, group_index]
    df_plot_4 = sdf.groupby(['Subject', 'NetMHCIIpan Allele'])['Mass'].median()
    df_plot_4.index = df_plot_4.index.map(lambda x: ' - '.join(x))
    df_plot_4 = df_plot_4.loc[x_val]
    ax.bar(x_val, df_plot_4.values)
    x_val = x_val.str.replace('HLA-', '')
    x_val = x_val.str.replace('_', '')
    ax.set_xticklabels(x_val, rotation=90, fontsize=8)
    ax.set_ylabel('Median glycan mass (Da)')
    ax.set_xlabel('HLA group')
    min1, max1 = axes[0, group_index].get_xlim()
    min2, max2 = axes[1, group_index].get_xlim()
    axes[0, group_index].set_xlim(min(min1, min2), max(max1, max2))
    axes[1, group_index].set_xlim(min(min1, min2), max(max1, max2))
    group_index += 1

plt.savefig(os.path.join(
    out, 'plots/figure_2/panel_relativePositionPerSubject.pdf'),
            transparent=True)

# Check glycan types per HLA group
df_plot_5 = round(
    df_glyco.groupby('HLA Group').apply(lambda x: x[
        'Glycan Type'].value_counts(normalize=True)).unstack(level=1), 4) * 100
cols_order = [
    "Truncated", "High Mannose", "Truncated Complex/Hybrid", "Complex/Hybrid",
    "Other"
]
fig, ax = plt.subplots(1,
                       1,
                       figsize=(3 * cm, 5.8 * cm),
                       constrained_layout=True)
colors = ["#e63946", "#f1faee", "#a8dadc", "#457b9d", "#1d3557"]
ax.set_ylabel('Glycan type (%)')
df_plot_5[cols_order].plot(kind='bar',
                           ax=ax,
                           stacked=True,
                           legend=True,
                           color=colors,
                           align='center')
plt.savefig(os.path.join(out, 'plots/figure_2/panel_glycanTypes.pdf'),
            transparent=True)

df_plot_6 = round(
    df_glyco.groupby([
        'HLA Group', 'Glyco position within motif'
    ]).apply(lambda x: x['Glycan Type'].value_counts(normalize=True)).unstack(
        level=[0, 2]), 4) * 100
cols_order = [
    "Truncated", "High Mannose", "Truncated Complex/Hybrid", "Complex/Hybrid",
    "Other"
]

df_plot_6 = df_plot_6.fillna(0)

fig, ax = plt.subplots(1, 1, figsize=(12, 4), constrained_layout=True)
ax.set_ylabel('Glycan type (%)')

plots = []
for group, pos in zip(['DP', 'DQ', 'DR'], [-0.2, 0, 0.2]):
    sdf = df_plot_6.loc[:, group][cols_order]
    for i in range(0, len(cols_order)):
        x_vals = sdf.index.astype(int).values + pos
        y_vals = sdf[cols_order[i]].values
        bot = sdf[cols_order[:i]].sum(axis=1).values
        if i == 0:
            plots.append(
                ax.bar(x_vals,
                       y_vals,
                       width=0.1,
                       align='center',
                       color=colors[i]))
        else:
            plots.append(
                ax.bar(x_vals,
                       y_vals,
                       bottom=bot,
                       width=0.1,
                       align='center',
                       color=colors[i]))

min_val = int(df_plot_6.index.min())
max_val = int(df_plot_6.index.max())
ax.legend(plots[:len(cols_order)], cols_order)
ax.set_xticks(np.arange(min_val, max_val + 1, 1))
ax.set_xticklabels([str(x) for x in range(min_val, max_val + 1)])
plt.savefig(os.path.join(out, 'plots/figure_2/panel_glycanTypes_2.pdf'),
            transparent=True)
