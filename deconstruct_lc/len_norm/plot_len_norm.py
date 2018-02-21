import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from deconstruct_lc import read_config


class PlotLenNorm(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        pdb_an_dp = os.path.join(data_dp, 'pdb_analysis')
        self.fpi = os.path.join(pdb_an_dp, 'pdb_len_norm.tsv')
        self.lc_m = 0.066213297264721263
        self.lc_b = 1.7520712972708843
        self.lc_b_up = 16.5
        self.grey_b = 36.5

    def plot_all(self):
        fig = plt.figure(figsize=(10, 3))
        ax1 = fig.add_subplot(121)
        self.plot_scatter(ax1)
        ax2 = fig.add_subplot(122)
        self.plot_nomiss(ax2)
        ax2.yaxis.set_label_position("right")
        fig.subplots_adjust(hspace=0, wspace=0)
        fig.tight_layout()
        plt.show()

    def plot_scatter(self, ax):
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        lens = list(df['Length'])
        raw_scores = list(df['score'])
        ax.scatter(lens, raw_scores, alpha=0.1, color='darkblue')
        self.plot_lines()
        ax.set_xlim([0, 1500])
        ax.set_ylim([0, 150])
        ax.set_xlabel('Protein sequence length', size=12)
        ax.set_ylabel('LC Motifs', size=12)

    def plot_nomiss(self, ax):
        """Show that the PDB norm dataset moves below the trendline when you
        don't count missing residues"""
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        lens = list(df['Length'])
        raw_scores = list(df['nomiss_score'])
        ax.scatter(lens, raw_scores, alpha=0.1, color='darkblue')
        self.plot_lines()
        ax.set_xlim([0, 1500])
        ax.set_ylim([0, 150])
        ax.set_xlabel('Protein sequence length', size=12)
        ax.set_ylabel('LC Motifs - No Missing Residues', size=12)
        plt.tick_params(axis='both', left='on', top='on', right='on',
                        bottom='on', labelleft='off', labeltop='off',
                        labelright='on', labelbottom='on')

    def plot_lines(self):
        x = np.arange(0, 1500, 0.01)
        y1 = self.plot_line(self.lc_m, self.lc_b, x)
        y2 = self.plot_line(self.lc_m, self.lc_b_up, x)
        y3 = self.plot_line(self.lc_m, self.grey_b, x)
        plt.plot(x, y1, color='black', lw=2)
        plt.plot(x, y2, color='black', lw=2, linestyle='--')
        plt.plot(x, y3, color='grey', lw=2)
        plt.xlim([0, 1500])
        plt.ylim([0, 150])

    def plot_line(self, m, b, x):
        y = m*x + b
        return y


def main():
    pln = PlotLenNorm()
    pln.plot_all()


if __name__ == '__main__':
    main()