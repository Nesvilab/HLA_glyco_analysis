import os
from collections import defaultdict

import logomaker as lm
import pandas as pd
import seaborn as sns
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib_venn import venn2

# for visual settings and compatibility with adobe illustrator
import plot_params
from general_params_glycFDR import out, pwm_annotation_path, tables
from helper_functions import get_hla_binding

# /storage/data/HLA/_Results/reports_glycFDR/deconvolution/Subjects/3830-NJF/Responsibilities


def read_deconvolution():
    pwm_annotation = pd.read_csv(pwm_annotation_path,
                                 sep='\t',
                                 header=0,
                                 low_memory=False)

    df_list = []
    pwm_annotation = pwm_annotation[~pwm_annotation['Number Of Motifs'].isnull(
    )]
    files = list(out + '/deconvolution/Subjects/' + pwm_annotation['Subject'] +
                 '/Responsibilities/bestPepResp_K' +
                 pwm_annotation['Number Of Motifs'].astype(int).astype(str) +
                 '.txt')
    pwm_annotation['File Location'] = files
    for file_path in files:
        loop_df = pd.read_csv(file_path, sep='\t', header=0, low_memory=False)
        loop_df['File Location'] = file_path
        df_list.append(loop_df)

    df_pwm_peptides = pd.concat(df_list)
    df_pwm_peptides.drop(['resp_2nd', 'Motif_2nd', 'P1_2nd', 'Core_2nd'],
                         axis=1,
                         inplace=True)
    df_pwm_peptides['Motif ID'] = df_pwm_peptides['Motif_1st'].tolist()
    df_pwm_peptides = df_pwm_peptides.merge(
        pwm_annotation[['File Location', 'Motif ID', 'Motif Allele']],
        how='left')
    df_pwm_peptides['Motif Allele'].fillna('flat motif', inplace=True)
    df_pwm_peptides['Subject'] = df_pwm_peptides['File Location'].str.rsplit(
        '/', expand=True)[8]
    return df_pwm_peptides


# read psm table
psm = pd.read_csv(tables['Purcell_20210_PXD025877'], sep='\t', header=0)
psm['Study'] = 'Purcell_2021_PXD025877'

# End of temporary solution #

# add glyco site information
psm_glyco = psm[~psm["Glycan Score"].isnull()].copy(deep=True)
psm_glyco['Assigned Modifications'] = psm_glyco[
    'Assigned Modifications'].str.split(', ')
psm_glyco = psm_glyco.explode('Assigned Modifications')
psm_glyco = psm_glyco[psm_glyco['Assigned Modifications'].str.contains(
    '^[0-9]+N')].copy(deep=True)
psm_glyco['Glyco Position'] = psm_glyco['Assigned Modifications'].str.extract(
    '^([0-9]+)N').astype(int)
psm_glyco['Protein glyco position'] = psm_glyco['Protein Start'] + psm_glyco[
    'Glyco Position'] - 1
psm_glyco['Glyco site'] = psm_glyco['Entry Name'] + '_' + psm_glyco[
    'Protein glyco position'].map(str)

psm_glyco = psm_glyco[psm_glyco['Glycan q-value'] <= 0.05]

psm = pd.concat([psm[psm["Glycan Score"].isnull()], psm_glyco])

# add sample information
ann = pd.read_csv(tables['annotation'], sep='\t', header=0)
ann = ann[ann['HLA Class'] == 2]

psm['File'] = psm.Spectrum.str.rsplit('.', 3, expand=True)[0] + '.mzML'

if any(psm['Study'].str.contains('_glycFDR')):
    ann['Study'] = ann['Study'] + '_glycFDR'

psm = psm.merge(ann)

psm = psm[~psm['Entry Name'].str.contains('rev_')]

df_deconv = read_deconvolution()

# information about the mono-allelic cell lines
# C1R_DR  DRA1*01:02;DRB1*12:01:01;DRB3*02:02:01
# C1R_DQ                DQA1*05:05;DQB1*03:01:01
# C1R_DP     DPB1*04:01:01;DPA1*02:01;DPA1*02:02

temp = psm.merge(df_deconv)
c1r = temp[temp['Subject'].str.contains('C1R')].copy(deep=True)
c1r['Core start'] = c1r.apply(lambda x: x.Peptide.index(x['Core_1st']),
                              axis=1) + 1

# add relative glucosylation position tot the HLA motif
c1r['Glyco position within motif'] = c1r['Glyco Position'] - c1r['Core start']

# yapf: disable
comp = c1r.groupby(['Study', 'Sample Name']).apply(
    lambda x: pd.concat(
        [
            x.drop_duplicates('Spectrum')['Motif Allele'].value_counts(normalize=True),
            x[~x['Glycan Score'].isnull()].drop_duplicates('Spectrum')['Motif Allele'].value_counts(normalize=True),
            x.drop_duplicates('Peptide')['Motif Allele'].value_counts(normalize=True),
            x[~x['Glycan Score'].isnull()].drop_duplicates('Peptide')['Motif Allele'].value_counts(normalize=True)
        ],
    axis=1)
)

# yapf: enable

comp.columns = ['PSM', 'Glyco PSM', 'Peptide', 'Glyco peptide']

# Plot the percetnage of peptides vs glyco-peptides that belong the HLA motifs
# along with the HLa motifs themselves
df_plot_1 = comp.loc[comp.index.map(lambda x: x[2] != 'flat motif'
                                    ), ].reset_index(level=0, drop=True) * 100

temp1 = df_plot_1.index.get_level_values(0).str.extractall(
    "(D[P|Q|R])rep([0-9]+)")
temp2 = df_plot_1.index.get_level_values(1)
samples = temp1[0] + ' rep ' + temp1[1]
df_plot_1['Motif Allele'] = temp2.tolist()
df_plot_1.index = samples
df_plot_1['Subject'] = ['C1R_' + x for x in temp1[0].tolist()]

fig, axes = plt.subplots(4, 1, figsize=(1.5, 5.6), constrained_layout=True)
axes = axes.flatten()

ax_index = 0
for (subject, allele), sdf in df_plot_1.groupby(['Subject', 'Motif Allele']):
    if allele == 'flat motif':
        continue
    ax = axes[ax_index]
    ax.set_title(allele, fontsize=8)
    sdf[["Peptide", "Glyco peptide"]].plot(kind='bar',
                                           ax=ax,
                                           color=["#58A4B0", "#1B1B1E"],
                                           legend=False)
    ax.set_ylim(0, 100)
    ax.set_xlabel('')
    ax_index += 1
    ax.set_ylabel('Peptides wihin\nthe HLA motif (%)')

plt.savefig(os.path.join(out,
                         'plots/figure_2/panel_binding_deconvolution.pdf'),
            transparent=True)

# Fisher test as additional test to check if the glycosylation status
# and HLA motifs are dependant

c1r['Motif type'] = c1r['Motif Allele'].str.contains('flat motif',
                                                     regex=True).map({
                                                         False:
                                                         'HLA motif',
                                                         True:
                                                         'flat motif'
                                                     })

# yapf: disable
comp_motif = c1r.groupby(['Study', 'Sample Name']).apply(lambda x:
        pd.concat([
            x.drop_duplicates('Spectrum')['Motif type'].value_counts(),
            x[~x['Glycan Score'].isnull()].drop_duplicates('Spectrum')['Motif type'].value_counts(),
            x.drop_duplicates('Peptide')['Motif type'].value_counts(),
            x[(~x['Glycan Score'].isnull()) ].drop_duplicates('Peptide')['Motif type'].value_counts()
        ], axis=1)
                                                        )

# yapf: enable
from scipy import stats

comp_motif.columns = ['PSM', 'Glyco PSM', 'Peptide', 'Glyco peptide']

comp_motif_stat = comp_motif[['Peptide', 'Glyco peptide']].fillna(0).groupby(
    'Sample Name').apply(lambda x: stats.fisher_exact(x)).apply(pd.Series)

comp_motif_stat.columns = ['oddRatio', 'p-value']
comp_motif_stat['H0'] = (comp_motif_stat['p-value'] > 0.05).map({
    True:
    'Not rejected',
    False:
    'Rejected'
})

### without a-value filter ###
# >>> comp_motif_stat
#                    oddRatio   p-value            H0
# Sample Name
# C1R_DPrep1_HLA-II  1.537538  0.297880  Not rejected
# C1R_DPrep2_HLA-II  0.256948  0.250347  Not rejected
# C1R_DPrep3_HLA-II  1.703258  0.202349  Not rejected
# C1R_DPrep4_HLA-II  0.376774  0.230825  Not rejected
# C1R_DPrep5_HLA-II  2.301692  0.028162      Rejected
# C1R_DQrep1_HLA-II  3.596905  0.065038  Not rejected
# C1R_DQrep2_HLA-II  1.263228  0.557876  Not rejected
# C1R_DQrep3_HLA-II  1.916667  0.422365  Not rejected
# C1R_DQrep4_HLA-II  0.000000  1.000000  Not rejected
# C1R_DQrep5_HLA-II  0.000000  1.000000  Not rejected
# C1R_DQrep6_HLA-II  3.958533  0.006253      Rejected
# C1R_DRrep1_HLA-II  0.000000  0.650053  Not rejected
# C1R_DRrep2_HLA-II  0.000000  1.000000  Not rejected
# C1R_DRrep3_HLA-II  0.988687  1.000000  Not rejected
# C1R_DRrep4_HLA-II  0.364277  0.037611      Rejected
# C1R_DRrep5_HLA-II  0.924388  1.000000  Not rejected
# C1R_DRrep6_HLA-II  1.360831  0.197725  Not rejected
# C1R_DRrep7_HLA-II  1.176897  0.650026  Not rejected

### with q-value filter ###
#                    oddRatio   p-value            H0
# Sample Name
# C1R_DPrep1_HLA-II  0.788162  1.000000  Not rejected
# C1R_DPrep2_HLA-II  0.000000  0.251400  Not rejected
# C1R_DPrep3_HLA-II  1.728587  0.304790  Not rejected
# C1R_DPrep4_HLA-II  0.000000  0.068860  Not rejected
# C1R_DPrep5_HLA-II  0.389085  0.507035  Not rejected
# C1R_DQrep1_HLA-II  4.289125  0.096800  Not rejected
# C1R_DQrep2_HLA-II  0.000000  1.000000  Not rejected
# C1R_DQrep3_HLA-II  0.000000  1.000000  Not rejected
# C1R_DQrep4_HLA-II  0.000000  1.000000  Not rejected
# C1R_DQrep5_HLA-II  0.000000  1.000000  Not rejected
# C1R_DQrep6_HLA-II  1.071863  0.612612  Not rejected
# C1R_DRrep1_HLA-II  0.000000  1.000000  Not rejected
# C1R_DRrep2_HLA-II  0.000000  1.000000  Not rejected
# C1R_DRrep3_HLA-II  0.000000  0.633408  Not rejected
# C1R_DRrep4_HLA-II  0.231214  0.132233  Not rejected
# C1R_DRrep5_HLA-II  0.000000  0.421321  Not rejected
# C1R_DRrep6_HLA-II  1.005580  0.865762  Not rejected
# C1R_DRrep7_HLA-II  0.804314  1.000000  Not rejected

# plot relative glycosylation position
cond = c1r['Motif type'] == 'HLA motif'
rel_psm = c1r[cond].groupby('Motif Allele').apply(
    lambda x: x['Glyco position within motif'].value_counts()).unstack(level=1)

rel_pep = c1r[cond].groupby('Motif Allele').apply(lambda x: x.drop_duplicates(
    'Peptide')['Glyco position within motif'].value_counts()).unstack(level=1)

fig, axes = plt.subplots(8,
                         1,
                         figsize=(2.1, 6),
                         sharex=True,
                         constrained_layout=True)
ax_index = 0


def highlight_motif(ax):
    xticks = ax.get_xticklabels()
    highlight_colors = []
    for i in xticks:
        if i._x >= 0 and i._x < 9:
            highlight_colors.append("#E71D36")
        else:
            highlight_colors.append("black")
    for xtick, color in zip(ax.get_xticklabels(), highlight_colors):
        xtick.set_color(color)
    for xtick, color in zip(ax.get_xticklines(), highlight_colors):
        xtick.set_color(color)


for allele in rel_psm.index:
    # plot motif
    ax = axes[ax_index]
    ax.set_title(allele, fontsize=8)
    sdf1 = c1r[(c1r['Motif Allele'] == allele)
               & (c1r['Motif type'] == 'HLA motif')]
    sequences = sdf1.drop_duplicates('Peptide')['Core_1st'].tolist()
    sdf1 = sdf1.dropna(subset=['Glyco position within motif'])
    sdf1['Glyco position within motif'] = sdf1[
        'Glyco position within motif'].astype(int)
    # generate matrix for the motif plot
    d = defaultdict(lambda: [0] * 9)  # d[char] = [pos0, pos12, ...]
    for seq in sequences:
        for i, char in enumerate(seq):
            d[char][i] += 1
    count_mat = pd.DataFrame(d)
    count_mat.index.names = ["pos"]
    info_mat = lm.transform_matrix(count_mat,
                                   from_type='counts',
                                   to_type='information')
    plot = lm.Logo(info_mat, ax=ax)
    ax.set_title(allele, fontsize=8)
    ax.get_xaxis().tick_bottom()
    # ax = axes[ax_index + 1]
    # df_plot1 = rel_psm.loc[allele]  # .apply(lambda x: x / x.sum())
    # df_plot1.index = df_plot1.index.astype(int)
    # ax.bar(x=df_plot1.index, height=df_plot1.values, width=0.5)
    # ax.tick_params(labelbottom=True)
    ax_index += 1
    ax = axes[ax_index]
    df_plot2 = rel_pep.loc[allele]  # .apply(lambda x: x / x.sum())
    df_plot2.index = df_plot2.index.astype(int)
    ax.bar(x=df_plot2.index, height=df_plot2.values, width=0.5)
    ax.tick_params(labelbottom=True)
    # make motif aligned with the bar plots
    ax_index += 1

min_val = int(c1r['Glyco position within motif'].min())
max_val = int(c1r['Glyco position within motif'].max())
axes[-1].set_xlim(min_val - 0.5, max_val)
axes[-1].set_xticks(range(min_val, max_val))
ticklabels = [str(x) if x % 2 == 0 else '' for x in range(min_val, max_val)]
for ax_index in range(1, 8, 2):
    axes[ax_index].set_xticklabels(ticklabels, fontsize=8, rotation=90)
    highlight_motif(axes[ax_index])
    axes[ax_index].set_ylabel('Peptide count')
    axes[ax_index].set_xlabel('Glycosylation relative position')

for ax_index in range(0, 8, 2):
    axes[ax_index].axes.get_xaxis().set_visible(False)
    axes[ax_index].tick_params(labelbottom=True)
    axes[ax_index].set_ylabel('Bits')

out_file = os.path.join(out,
                        f'plots/figure_2/panel_deconv_relative_position.pdf')
plt.savefig(out_file, transparent=True)

# plot absolute glycosylation positions
fig, axes = plt.subplots(4, 1, figsize=(2, 5.6), constrained_layout=True)
axes = axes.flatten()
c1r_glyco = c1r.dropna(subset=['Glyco Position']).copy(deep=True)
c1r_glyco['Glyco Position'] = c1r_glyco['Glyco Position'].astype(int)

ax_index = 0

for group, sdf in c1r_glyco.groupby('Motif Allele'):
    if group == 'flat motif':
        continue
    ax = axes[ax_index]
    matrix = sdf.drop_duplicates('Peptide').groupby(
        c1r_glyco.Peptide.str.len())['Glyco Position'].value_counts().unstack(
            'Peptide')
    start_pos = min(matrix.index.min(), matrix.columns.min())
    end_pos = max(matrix.index.max(), matrix.columns.max())
    set_index = set(matrix.index)
    set_cols = set(matrix.columns)
    rx = range(start_pos, end_pos + 1)
    idx = list(set(rx).difference(set_index))
    cols = list(set(rx).difference(set_cols))
    matrix.loc[:, cols] = 0
    for i in idx:
        matrix.loc[i, :] = 0
    matrix.sort_index(inplace=True)
    matrix = matrix[matrix.columns.sort_values()]
    matrix = matrix.T
    sb.heatmap(matrix, ax=ax, cmap="Blues")
    ax.set_ylabel('Peptide length')
    ax.set_xlabel('Glycosylation\nabsolute\nposition')
    ax.set_title(group, fontsize=8)
    cterm_color = '#011627'
    nterm_color = '#A2a7A5'
    cterm = []
    nterm = []
    for i in range(end_pos):
        cterm.append([i + 1, i])
        cterm.append([i + 1, i + 1])
    ax.plot([x[0] for x in cterm], [x[1] for x in cterm],
            lw=0.5,
            color=cterm_color)
    ax.vlines(0,
              ymin=ax.get_ylim()[0],
              ymax=ax.get_ylim()[1],
              lw=0.75,
              color=nterm_color)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_yticks(range(1, 26, 4))
    ax.set_xticks(range(1, 26, 4))
    ax.set_yticklabels([str(x) for x in range(1, 26, 4)])
    ax.set_xticklabels([str(x) for x in range(1, 26, 4)])
    ax.collections[0].colorbar.set_label("Peptide count", fontsize=8)
    ax_index += 1

out_file = os.path.join(out,
                        f'plots/figure_2/panel_deconv_absolute_position.pdf')
plt.savefig(out_file, transparent=True)
