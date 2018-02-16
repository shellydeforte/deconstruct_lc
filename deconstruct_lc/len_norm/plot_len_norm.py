import configparser
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from deconstruct_lc import tools_lc
from deconstruct_lc import read_config
from deconstruct_lc.len_norm.len_norm import LenNorm

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class PlotLenNorm(object):
    def __init__(self):
        config = read_config.read_config()
        self.data_dp = os.path.join(config['fps']['data_dp'])
        self.pdb_dp = os.path.join(self.data_dp, 'pdb_prep')
        self.norm_fpi = os.path.join(self.pdb_dp, 'pdb_norm_cd100.tsv')
        self.k = 6
        self.lca = 'SGEQAPDTNKR'
        self.lce = 1.6
        self.lc_m = 0.066213297264721263
        self.lc_b = 1.7520712972708843
        self.lc_b_up = 16.5
        self.grey_b = 36.5

    def demonstrate_mb(self):
        """
        LinregressResult(slope=0.066213297264721263, intercept=1.7520712972708843,
        rvalue=0.7809848592680475, pvalue=0.0, stderr=0.0002720887990701889)
        """
        ln = LenNorm()
        lr = ln.mb_lc(self.k, self.lca, self.lce)
        print(lr)

    def plot_all(self):
        plt.subplot(1, 2, 1)
        self.plot_scatter()
        plt.subplot(1, 2, 2)
        self.plot_nomiss()
        plt.subplots_adjust(hspace=0.5)
        plt.show()

    def plot_a_line(self):
        x = np.arange(0, 1500, 0.01)
        y1 = self.plot_line(self.lc_m, self.lc_b, x)
        y2 = self.plot_line(self.lc_m, self.lc_b_up, x)
        y3 = self.plot_line(self.lc_m, self.grey_b, x)
        plt.plot(x, y1, color='black', lw=2)
        plt.plot(x, y2, color='black', lw=2, linestyle='--')
        plt.plot(x, y3, color='grey', lw=2)
        plt.xlim([0, 1500])
        plt.ylim([0, 150])

    def plot_scatter(self):
        df = pd.read_csv(self.norm_fpi, sep='\t', index_col=0)
        seqs = df['Sequence']
        lens = [len(seq) for seq in seqs]
        raw_scores = tools_lc.calc_lc_motifs(seqs, self.k, self.lca, self.lce)
        plt.scatter(lens, raw_scores, alpha=0.1, color='darkblue')
        x = np.arange(0, 1500, 0.01)
        y1 = self.plot_line(self.lc_m, self.lc_b, x)
        y2 = self.plot_line(self.lc_m, self.lc_b_up, x)
        y3 = self.plot_line(self.lc_m, self.grey_b, x)
        plt.plot(x, y1, color='black', lw=2)
        plt.plot(x, y2, color='black', lw=2, linestyle='--')
        plt.plot(x, y3, color='grey', lw=2)
        plt.xlim([0, 1500])
        plt.ylim([0, 150])
        plt.xlabel('Protein sequence length', size=12)
        plt.ylabel('Raw LC score', size=12)
        #plt.show()

    def plot_nomiss(self):
        """Show that the PDB norm dataset moves below the trendline when you
        don't count missing residues"""
        df = pd.read_csv(self.norm_fpi, sep='\t', index_col=0)
        seqs = df['Sequence']
        miss_seqs = df['Missing']
        lens = [len(seq) for seq in seqs]
        raw_scores = tools_lc.calc_lc_motifs_nomiss(seqs, miss_seqs, self.k,
                                                    self.lca, self.lce)
        plt.scatter(lens, raw_scores, alpha=0.1, color='darkblue')
        x = np.arange(0, 1500, 0.01)
        y1 = self.plot_line(self.lc_m, self.lc_b, x)
        y2 = self.plot_line(self.lc_m, self.lc_b_up, x)
        y3 = self.plot_line(self.lc_m, self.grey_b, x)
        plt.plot(x, y1, color='black', lw=2)
        plt.plot(x, y2, color='black', lw=2, linestyle='--')
        plt.plot(x, y3, color='grey', lw=2)
        plt.xlim([0, 1500])
        plt.ylim([0, 150])
        plt.xlabel('Protein sequence length', size=12)
        plt.ylabel('Raw LC score', size=12)
        #plt.show()

    def plot_line(self, m, b, x):
        y = m*x + b
        return y


def main():
    pln = PlotLenNorm()
    pln.plot_all()


if __name__ == '__main__':
    main()