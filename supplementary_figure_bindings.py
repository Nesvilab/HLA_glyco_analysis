import os
from collections import defaultdict

import logomaker as lm
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# for visual settings and compatibility with adobe illustrator
import plot_params
from general_params_glycFDR import out

mpl.rcParams['font.size'] = 5

cm = 1 / 2.54  # centimeters in inches


def plot_motif(sequences, ax_):
    d = defaultdict(lambda: [0] * 9)  # d[char] = [pos0, pos12, ...]
    for seq in sequences:
        for i, char in enumerate(seq):
            d[char][i] += 1
    count_mat = pd.DataFrame(d)
    count_mat.index.names = ["pos"]
    info_mat = lm.transform_matrix(count_mat,
                                   from_type='counts',
                                   to_type='information')
    lm.Logo(info_mat, ax=ax_, color_scheme='weblogo_protein')


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


def adjust_axes(ax1, ax2, min_val, max_val, col_index, xlab):
    ax1.set_xlim(min_val - 0.5, max_val)
    ax2.set_xlim(min_val - 0.5, max_val)
    ticklabels = [str(x) for x in [-13, 0, 8, 16]]
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    ax1.axes.get_xaxis().set_visible(False)
    ax1.spines.bottom.set_visible(False)
    ax2.set_xticks([-13, 0, 8, 16])
    ax2.set_xticklabels(ticklabels)
    if col_index == 0:
        ax1.set_ylabel('Bits')
        ax2.set_ylabel('Glyco\nPeptide')
    ax2.set_xlabel(f'Glyco relative position\n{xlab}')


def plot_figures(plot_df, psm_glyco, out_dir):
    nrows = 10  # must be a multiple of 2
    ncols = 6
    nfig = 1
    width = 18 * cm
    height = 18.5 * cm
    fig, axes = plt.subplots(nrows,
                             ncols,
                             figsize=(width, height),
                             constrained_layout=True)
    row_index = 0
    col_index = 0
    prev_allele = ""
    min_sequences = 10
    allele_index = 0
    for allele, sdf in plot_df.groupby('Allele'):
        print(f"Allele: {allele}")
        skipped = True
        for subject, ssdf in sdf.groupby('Subject'):
            if ssdf.shape[0] < min_sequences:
                print(f"Allele: {allele}, Subject: {subject} skipped")
                continue
            # fetch glycosylation relative positions
            glyco_df = psm_glyco[psm_glyco.Subject == subject].merge(
                ssdf[['Peptide', 'Core']])
            if glyco_df.empty:
                continue
                print(f"Allele: {allele}, Subject: {subject} skipped")
            skipped = False
            print(f"Allele: {allele}, Subject: {subject} processed")
            glyco_df['Core start'] = glyco_df.apply(
                lambda x: x.Peptide.index(x.Core), axis=1) + 1
            glyco_df['Glyco position within motif'] = glyco_df[
                'Glyco Position'] - glyco_df['Core start']
            df_plot = glyco_df['Glyco position within motif'].value_counts()
            df_plot = df_plot.sort_index()
            df_plot.index = df_plot.index.astype(int)
            # plot motif
            ax1 = axes[row_index, col_index]
            peptides = list(ssdf.Core.values)
            title = allele.replace('-', '\n').replace('HLA', '')
            if prev_allele != title:
                ax1.set_title(title,
                              fontsize=5,
                              fontweight="bold",
                              bbox=dict(facecolor='none',
                                        edgecolor='black',
                                        boxstyle='round'))
            plot_motif(peptides, ax1)
            # plot glycosylation relative positions
            ax2 = axes[row_index + 1, col_index]
            if allele_index % 2 == 0:
                print(
                    f"{allele} with index: {allele_index} set grey background")
                ax1.patch.set_facecolor('grey')
                ax1.patch.set_alpha(0.2)
                ax2.patch.set_facecolor('grey')
                ax2.patch.set_alpha(0.2)
            colors = [
                'red' if x >= 0 and x <= 8 else 'gray' for x in df_plot.index
            ]
            ax2.bar(df_plot.index, df_plot.values, color=colors)
            xlabel = subject + f" n: {ssdf.shape[0]}"
            adjust_axes(ax1, ax2, -13, 16, col_index, xlabel)
            prev_allele = title
            if row_index < nrows - 2:
                row_index += 2
            elif col_index < ncols - 1:
                row_index = 0
                col_index += 1
            else:
                print(f"Saving figure {nfig}")
                plt.savefig(os.path.join(out_dir, f'figure_{nfig}.png'),
                            dpi=300)
                plt.savefig(os.path.join(out_dir, f'figure_{nfig}.pdf'),
                            transparent=True)
                nfig += 1
                fig, axes = plt.subplots(nrows,
                                         ncols,
                                         figsize=(width, height),
                                         constrained_layout=True)
                row_index = 0
                col_index = 0
        if not skipped:
            print("incrementing index")
            allele_index += 1
        print('')
    for j in range(col_index, ncols):
        if j == col_index:
            for i in range(row_index, nrows):
                axes[i, j].remove()
        else:
            for i in range(0, nrows):
                axes[i, j].remove()

    print(f"Saving figure {nfig}\n")
    plt.savefig(os.path.join(out_dir, f'figure_{nfig}.png'), dpi=300)
    plt.savefig(os.path.join(out_dir, f'figure_{nfig}.pdf'), transparent=True)


def main():
    df_bindings = pd.read_csv(os.path.join(
        out, 'tables/binding_predictions_best-allele.tsv'),
                              sep='\t',
                              header=0)

    df_bindings = df_bindings[df_bindings['Binding Type'] == 'Non-trash']

    # generated by './figure_2_glycosites.py'
    glyco_path = os.path.join(out, 'tables/glycosites.tsv')

    # reading tables for all studies
    psm_glyco = pd.read_csv(glyco_path, sep='\t', header=0, low_memory=False)
    psm_glyco['Glycosylated'] = True

    df_bindings = df_bindings.merge(
        psm_glyco[['Peptide', 'Subject', 'Glycosylated']].drop_duplicates(),
        how='left')
    df_bindings['Glycosylated'] = df_bindings['Glycosylated'].fillna(False)
    output_path = os.path.join(out, 'plots/sup_fig')
    plot_figures(df_bindings, psm_glyco, output_path)


if __name__ == '__main__':
    main()
