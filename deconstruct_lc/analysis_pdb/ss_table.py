"""Results: Neither the table nor the plots by bin are particularly compelling.
X goes up, but most other things are pretty similar without a lot of movement
Structure is about a 15 point difference"""

import matplotlib.pyplot as plt
from scipy.interpolate import spline
import os
import numpy as np
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import motif_seq


class SsTable(object):
    """Calculate the secondary structure inside and outside of LC motifs"""
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.pdb_dp = os.path.join(self.data_dp, 'pdb_prep')
        self.pdb_an_dp = os.path.join(self.data_dp,
                                      'pdb_analysis')
        self.an_fpi = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.miss_fp = os.path.join(self.pdb_an_dp, 'miss_in_out.tsv')
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'

    def read_ss(self):
        all_ss_in = ''
        all_ss_out = ''
        df = pd.read_csv(self.an_fpi, sep='\t', index_col=0)
        for i, row in df.iterrows():
            seq = row['Sequence']
            ind_in, ind_out = self.get_inds(seq)
            ss = row['Secondary Structure']
            miss = row['Missing']
            xss = self.add_x(ss, miss)
            ss_in, ss_out = self.get_ss(xss, ind_in, ind_out)
            all_ss_in += ss_in
            all_ss_out += ss_out
        all_ss = set(all_ss_in)
        ss_in_dict = {}
        ss_out_dict = {}
        for an_ss in all_ss:
            ss_in_dict[an_ss] = (all_ss_in.count(an_ss))/len(all_ss_in)
            ss_out_dict[an_ss] = (all_ss_out.count(an_ss))/len(all_ss_out)
        print(ss_in_dict)
        print(ss_out_dict)

    def get_ss(self, ss, ind_in, ind_out):
        ss_in = ''
        ss_out = ''
        for ii in ind_in:
            ss_in += ss[ii]
        for io in ind_out:
            ss_out += ss[io]
        return ss_in, ss_out

    def add_x(self, ss, miss):
        nss = ''
        for s, m in zip(ss, miss):
            if m == 'X':
                nss += m
            else:
                nss += s
        return nss

    def get_inds(self, seq):
        lcas = motif_seq.LcSeq(seq, self.k, self.lca, 'lca')
        lces = motif_seq.LcSeq(seq, self.k, self.lce, 'lce')
        lca_in, lca_out = lcas._get_motif_indexes()
        lce_in, lce_out = lces._get_motif_indexes()
        ind_in = lca_in.union(lce_in)
        ind_out = lca_out.union(lce_out)
        return ind_in, ind_out


class PlotSs(object):
    """For each bin, take the in/out regions as whole sequences and then
    calculate the fraction"""
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.pdb_dp = os.path.join(self.data_dp, 'pdb_prep')
        self.pdb_an_dp = os.path.join(self.data_dp,
                                      'pdb_analysis')
        self.an_fpi = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.ss_out_fp = os.path.join(self.pdb_an_dp, 'ss_out.tsv')
        self.ss_in_fp = os.path.join(self.pdb_an_dp, 'ss_in.tsv')
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'

    def plot_bar(self):
        data = [[0.2, 0.4, 0.4], [0.1, 0.5, 0.6]]
        data1 = [0.2, 0.1]
        data2 = [0.4, 0.5]
        data3 = [0.4, 0.6]
        plt.bar([0, 1], data1, color='darkblue')
        plt.bar([0, 1], data2, bottom=data1, color='darkorange')
        plt.show()

    def bar_plot(self):
        df_out = pd.read_csv(self.ss_out_fp, sep='\t', index_col=0)
        df_in = pd.read_csv(self.ss_in_fp, sep='\t', index_col=0)
        x = list(range(0, 10))
        x2 = [i+0.25 for i in x]
        self.abar(df_in, x2, 1)
        self.abar(df_out, x, 1)
        plt.show()

    def abar(self, df, x, a):
        missing = []
        noss = []
        turns = []
        struct = []
        for i, row in df.iterrows():
            # each row is a bin
            missing.append(row['X'])
            noss.append(row['P'])
            turns.append((row['S'] + row['T']))
            struct.append((row['E'] + row['H'] + row['B'] + row['G'] + row['I']))
        plt.bar(x, struct, color='black', width=0.2, alpha=a)
        plt.bar(x, turns, color='grey', bottom=struct, width=0.2, alpha=a)
        plt.bar(x, noss, color='darkgrey', bottom=np.array(turns)+np.array(struct), width=0.2, alpha=a)
        plt.bar(x, missing, color='darkred', bottom=np.array(turns)+np.array(struct)+np.array(noss), width=0.2, alpha=a)
        plt.ylim([0, 1.1])


    def read_plot(self):
        df = pd.read_csv(self.ss_out_fp, sep='\t', index_col=0)
        all_ss = ['P', 'X', 'T', 'S', 'H', 'E', 'B', 'G', 'I']
        ss_in_dict = {}
        x = list(range(0, 10))
        for ss in all_ss:
            ss_in_dict[ss] = list(df[ss])
            print(len(ss_in_dict[ss]))
            xnew = np.linspace(0, 10, 300)
            power_smooth = spline(x, ss_in_dict[ss], xnew)
            plt.plot(xnew, power_smooth, label=ss)
            # plt.plot(ss_out_dict[ss], linestyle='--')
        plt.ylim([0, 0.35])
        plt.legend()
        plt.show()

    def get_bins(self):
        df = pd.read_csv(self.an_fpi, sep='\t', index_col=0)
        bins = range(0, 50, 5)
        ss_in_dict = {'P': [], 'X': [], 'T': [], 'S': [], 'H': [], 'E': [],
                      'B': [], 'G': [], 'I': []}
        ss_out_dict = {'P': [], 'X': [], 'T': [], 'S': [], 'H': [], 'E': [],
                      'B': [], 'G': [], 'I': []}
        for i in list(bins):
            ndf = df[(df['LC Raw'] >= i) & (df['LC Raw'] < i + 5)]
            self.read_ss(ndf, ss_in_dict, ss_out_dict)

        print(ss_in_dict)
        print(ss_out_dict)
        df_ssin = pd.DataFrame(ss_in_dict)
        df_ssout = pd.DataFrame(ss_out_dict)
        df_ssin.to_csv(self.ss_in_fp, sep='\t')
        df_ssout.to_csv(self.ss_out_fp, sep='\t')

    def read_ss(self, df, ss_in_dict, ss_out_dict):
        all_ss_in = ''
        all_ss_out = ''
        for i, row in df.iterrows():
            seq = row['Sequence']
            ind_in, ind_out = self.get_inds(seq)
            ss = row['Secondary Structure']
            miss = row['Missing']
            xss = self.add_x(ss, miss)
            ss_in, ss_out = self.get_ss(xss, ind_in, ind_out)
            all_ss_in += ss_in
            all_ss_out += ss_out
        all_ss = ['P', 'X', 'T', 'S', 'H', 'E', 'B', 'G', 'I']
        for an_ss in all_ss:
            ss_in_dict[an_ss].append((all_ss_in.count(an_ss))/len(all_ss_in))
            ss_out_dict[an_ss].append((all_ss_out.count(an_ss))/len(all_ss_out))

    def get_ss(self, ss, ind_in, ind_out):
        ss_in = ''
        ss_out = ''
        for ii in ind_in:
            ss_in += ss[ii]
        for io in ind_out:
            ss_out += ss[io]
        return ss_in, ss_out

    def add_x(self, ss, miss):
        nss = ''
        for s, m in zip(ss, miss):
            if m == 'X':
                nss += m
            else:
                nss += s
        return nss

    def get_inds(self, seq):
        lcas = motif_seq.LcSeq(seq, self.k, self.lca, 'lca')
        lces = motif_seq.LcSeq(seq, self.k, self.lce, 'lce')
        lca_in, lca_out = lcas._get_motif_indexes()
        lce_in, lce_out = lces._get_motif_indexes()
        ind_in = lca_in.union(lce_in)
        ind_out = lca_out.union(lce_out)
        return ind_in, ind_out


def main():
    st = PlotSs()
    st.bar_plot()


if __name__ == '__main__':
    main()