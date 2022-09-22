# plot sample types in the dataset
import os

import matplotlib.pyplot as plt
import pandas as pd

import plot_params
from general_params import out, tables

df = pd.read_csv(tables['annotation'], sep='\t', header=0)

df = df[df['HLA Class'] == 2]

df['Sample Type'] = df['Sample Type'].str.replace(
    'Tumor infiltrating lymphocytes', 'TILs')
df['Disease'] = df['Disease'].str.replace(
    'EBV transformed lymphoblastoid cell line', 'EBV TLCL')
fig, axes = plt.subplots(1, 3, figsize=(6, 2), constrained_layout=True)

# colors used in adobe illustrator
colors = ['#6c757d', '#212529', '#adb5bd']

# plot hla typing availability
ax = axes[0]
df_plot_1 = (
    ~df.drop_duplicates("Sample Name")['HLA Allele'].isnull()).value_counts(
        normalize=True)
df_plot_1.plot(kind='pie',
               ax=ax,
               explode=[0, 0.1],
               startangle=290,
               colors=[colors[1]],
               labels=["With HLA typing", "Without HLA typing"],
               autopct='%1.1f%%')
ax.set_ylabel("")
ax.set_xlabel("")
# making percentage values text white
for text in ax.texts:
    if '%' in text.get_text():
        text.set_color('w')

# plot sample types (cell line or patient)
ax = axes[1]
df_plot_2 = df.drop_duplicates("Sample Name")['Sample Type'].value_counts(
    normalize=True)
df_plot_2.plot(kind='pie',
               ax=ax,
               explode=[0, 0.05, 0.1],
               startangle=300,
               colors=[colors[1]],
               labels=df_plot_2.index,
               autopct='%1.1f%%')
ax.set_ylabel("")
ax.set_xlabel("")
# making percentage values text white
for text in ax.texts:
    if '%' in text.get_text():
        text.set_color('w')

# plot cancer type per disease
ax = axes[2]
df_plot_3 = df.drop_duplicates("Sample Name").groupby(
    'Disease')['Sample Type'].value_counts().unstack('Sample Type').fillna(0)
df_plot_3 = df_plot_3 / df_plot_3.sum().sum()
df_plot_3 = df_plot_3.loc[df_plot_3.sum(axis=1).sort_values(
    ascending=False).index, :]
df_plot_3.plot(kind='bar', ax=ax, stacked=True, color=colors)
ax.set_xlabel('Disease type')
vals = ax.get_yticks()
ax.set_yticks(vals)
# lifting the bars a bit from the x axis
ax.set_ylim([-0.01, max(vals)])
# formatting the y-axis in percentage values
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.set_xlabel("")
ax.set_ylabel("Percentage of samples")

# Make some labels.
values = df_plot_3.sum(axis=1)
labels = list(df_plot_3.index)

for x, (y, label) in enumerate(zip(values, labels)):
    ax.text(x, y + 0.01, label, ha="center", va="bottom", rotation='vertical')

ax.set_xticklabels([])
plt.savefig(os.path.join(out, 'plots/figure_1/first_panel_piecharts.pdf'),
            transparent=True)

# Alleles panel
alleles = df['HLA Allele'].str.split(';').explode()
alleles_parsed = alleles.str.extract("(^[^\*]+\*[0-9]+:[0-9]+)")
alleles = pd.concat([alleles, alleles_parsed], axis=1)
alleles.columns = ["Full HLA Allele", 'HLA Allele']
alleles['HLA Type'] = alleles['HLA Allele'].str.extract('^([A-Z]{2,3})')
alleles.dropna(inplace=True)
alleles = alleles.sort_values('HLA Type')

width = 0.45
fig, axes = plt.subplots(3, 2, figsize=(6.7, 2.9), constrained_layout=True)
axes = axes.flatten()

sdf_plots = []
for index, (hlatype, sdf) in enumerate(alleles.groupby('HLA Type')):
    ax = axes[index]
    sdf_plots.append(sdf['HLA Allele'].value_counts())
    ax.set_title(sdf_plots[-1].index[0][:3], fontsize=8)
    sdf_plots[-1].index = sdf_plots[-1].index.map(lambda x: x[3:])
    sdf_plots[-1].plot(kind="bar", ax=ax)
    # lifting the bars a bit from the x axis
    ax.set_ylim([-sdf_plots[-1].max() * 0.1, sdf_plots[-1].max()])
    ax.set_xlim([-width, (2 + sdf_plots[-1].shape[0]) * width])

# set the same width for all axes
for i in range(len(sdf_plots)):
    axes[i].set_xlim([-.5, max([x.shape[0] for x in sdf_plots]) - 1 + width])

plt.savefig(os.path.join(out, 'plots/figure_1/second_panel.pdf'),
            transparent=True)
