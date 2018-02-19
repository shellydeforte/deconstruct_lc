"""Results: Neither the table nor the plots by bin are particularly compelling.
X goes up, but most other things are pretty similar without a lot of movement
Structure is about a 15 point difference"""

import matplotlib.pyplot as plt
from scipy.interpolate import spline
import os
import numpy as np
import pandas as pd
from collections import OrderedDict

from deconstruct_lc import read_config
from deconstruct_lc import motif_seq
from deconstruct_lc import tools_lc


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
        self.ss_one_in_fp = os.path.join(self.pdb_an_dp, 'ss_one_in.tsv')
        self.ss_one_out_fp = os.path.join(self.pdb_an_dp, 'ss_one_out.tsv')
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'

    def one_bar_plot(self):
        """Plot the average of all, but also mention that the bars go down by bins"""
        df_out = pd.read_csv(self.ss_one_out_fp, sep='\t', index_col=0)
        df_in = pd.read_csv(self.ss_one_in_fp, sep='\t', index_col=0)

        x1 = [0]
        x2 = [0.3]
        h = [0.2]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_yticks([0.1, 0.4])

        self.data_one_bar(df_out, x1, 1, h)
        plt.legend(fontsize=12)
        self.data_one_bar(df_in, x2, 1, h)
        labels = ['Outside Motifs', 'Inside Motifs']

        ax.set_yticklabels(labels, size=12)
        plt.ylim([0, 1])
        plt.xlim([0, 1.0])
        plt.xlabel('Fraction Secondary Structure')
        plt.tight_layout()
        plt.show()

    def data_one_bar(self, df, x, a, w):
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
        bot1 = (np.array(turns) + np.array(struct) + np.array(noss))[0]
        bot2 = (np.array(turns) + np.array(struct))[0]
        plt.barh(x, missing, color='white', left=bot1, height=w, alpha=a, label='Missing')
        plt.barh(x, noss, color='darkgrey', left=bot2, height=w, alpha=a, label='Coils')
        plt.barh(x, turns, color='grey', left=struct, height=w, alpha=a, label='Turns and Bends')
        plt.barh(x, struct, color='black', height=w, alpha=a, label='Alpha Helix and Beta Sheet')

    def bar_plot(self):
        df_out = pd.read_csv(self.ss_out_fp, sep='\t', index_col=0)
        df_in = pd.read_csv(self.ss_in_fp, sep='\t', index_col=0)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x = list(range(0, 10))
        x2 = [i+0.45 for i in x]
        self.abar(df_in, x2, 1)
        ax.legend()
        #ax.legend(bbox_to_anchor=(1, 1.2))
        self.abar(df_out, x, 1)
        labels = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30',
                  '30-35', '35-40', '40-45', '45-50', '50+']
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, size=12)
        #plt.tight_layout()

        plt.ylim([0, 1.5])
        plt.xlim([-1, len(x)+1])
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
        plt.bar(x, struct, color='black', width=0.4, alpha=a, label='Alpha Helix and Beta Sheet')
        plt.bar(x, turns, color='grey', bottom=struct, width=0.4, alpha=a, label='Turns and Bends')
        plt.bar(x, noss, color='darkgrey', bottom=np.array(turns)+np.array(struct), width=0.4, alpha=a, label='Coils')
        plt.bar(x, missing, color='darkred', bottom=np.array(turns)+np.array(struct)+np.array(noss), width=0.4, alpha=a, label='Missing')
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

    def one_bin(self):
        df = pd.read_csv(self.an_fpi, sep='\t', index_col=0)
        ss_in_dict = {'P': [], 'X': [], 'T': [], 'S': [], 'H': [], 'E': [],
                      'B': [], 'G': [], 'I': []}
        ss_out_dict = {'P': [], 'X': [], 'T': [], 'S': [], 'H': [], 'E': [],
                      'B': [], 'G': [], 'I': []}
        self.read_ss(df, ss_in_dict, ss_out_dict)
        print(ss_in_dict)
        print(ss_out_dict)
        df_ssin = pd.DataFrame(ss_in_dict)
        df_ssout = pd.DataFrame(ss_out_dict)
        df_ssin.to_csv(self.ss_one_in_fp, sep='\t')
        df_ssout.to_csv(self.ss_one_out_fp, sep='\t')

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


class SsComp(object):
    """If it is both in motif and S/T/P/X - what do the motifs look like?"""
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.pdb_dp = os.path.join(self.data_dp, 'pdb_prep')
        self.pdb_an_dp = os.path.join(self.data_dp,
                                      'pdb_analysis')
        self.an_fpi = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.ss_out_fp = os.path.join(self.pdb_an_dp, 'ss_out.tsv')
        self.ss_in_fp = os.path.join(self.pdb_an_dp, 'ss_in.tsv')
        self.ss_one_in_fp = os.path.join(self.pdb_an_dp, 'ss_one_in.tsv')
        self.ss_one_out_fp = os.path.join(self.pdb_an_dp, 'ss_one_out.tsv')
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'

    def comp(self):
        df = pd.read_csv(self.an_fpi, sep='\t', index_col=0)
        all_kmers = {}
        for i, row in df.iterrows():
            print(i)
            seq = row['Sequence']
            ss = row['Secondary Structure']
            miss = row['Missing']
            xss = self.add_x(ss, miss)

            seq_kmers = tools_lc.seq_to_kmers(seq, self.k)
            ss_kmers = tools_lc.seq_to_kmers(xss, self.k)
            for seq_kmer, ss_kmer in zip(seq_kmers, ss_kmers):
                if tools_lc.lca_motif(seq_kmer, self.lca) or tools_lc.lce_motif(seq_kmer, self.lce):
                    if set(ss_kmer) <= {'S', 'T', 'P', 'X'}:
                        if seq_kmer in all_kmers:
                            all_kmers[seq_kmer] += 1
                        else:
                            all_kmers[seq_kmer] = 1
        for item in all_kmers:
            if all_kmers[item] > 200:
                print(item)
                print(all_kmers[item])


    def add_x(self, ss, miss):
        nss = ''
        for s, m in zip(ss, miss):
            if m == 'X':
                nss += m
            else:
                nss += s
        return nss

def main():
    pss = PlotSs()
    pss.one_bar_plot()


if __name__ == '__main__':
    main()