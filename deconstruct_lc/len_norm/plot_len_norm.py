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
        plt.subplot(1, 2, 1)
        self.plot_scatter()
        plt.subplot(1, 2, 2)
        self.plot_nomiss()
        plt.subplots_adjust(hspace=0.5)
        plt.show()

    def plot_scatter(self):
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        lens = list(df['Length'])
        raw_scores = list(df['score'])
        plt.scatter(lens, raw_scores, alpha=0.1, color='darkblue')
        self.plot_lines()
        plt.xlim([0, 1500])
        plt.ylim([0, 150])
        plt.xlabel('Protein sequence length', size=12)
        plt.ylabel('LC Motifs', size=12)

    def plot_nomiss(self):
        """Show that the PDB norm dataset moves below the trendline when you
        don't count missing residues"""
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        lens = list(df['Length'])
        raw_scores = list(df['nomiss_score'])
        plt.scatter(lens, raw_scores, alpha=0.1, color='darkblue')
        self.plot_lines()
        plt.xlim([0, 1500])
        plt.ylim([0, 150])
        plt.xlabel('Protein sequence length', size=12)
        plt.ylabel('LC Motifs - No Missing Residues', size=12)

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