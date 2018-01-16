import configparser
import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from deconstruct_lc import motif_seq
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class MissMotif(object):
    """
    Do missing residues occur more frequently within blobs?
    """
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.score_dp = os.path.join(config['filepaths']['data_dp'], 'scores')
        self.pdb_an_fp = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.k_lca = 6
        self.k_lce = 6
        self.alph_lca = 'SGEQAPDTNKR'
        self.thresh_lce = 1.6
        self.lca_label = '{}_{}'.format(self.k_lca, self.alph_lca)
        self.lce_label = '{}_{}'.format(self.k_lce, self.thresh_lce)

    def read_df(self):
        df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        bins = range(0, 20, 5)
        mean_mm = []
        std_mm = []
        mean_mp = []
        std_mp = []
        for i in bins:
            print(i)
            ndf = df[(df[self.lca_label] >= i) & (df[self.lca_label] < i+5)]
            miss_in_motifs, motif_percs = self.lc_blobs(ndf)
            mean_mm.append(np.mean(miss_in_motifs))
            std_mm.append(np.std(miss_in_motifs))
            mean_mp.append(np.mean(motif_percs))
            std_mp.append(np.std(motif_percs))
        print(mean_mm)
        print(std_mm)
        print(mean_mp)
        print(std_mp)
        plt.errorbar(bins, mean_mm, std_mm, linestyle='None', marker='o')
        plt.errorbar(bins, mean_mp, std_mp, linestyle='None', marker='o')
        plt.show()

    def lc_blobs(self, df):
        miss_in_motifs = []
        motif_percs = []
        for i, row in df.iterrows():
            miss = row['Missing']
            seq = row['Sequence']
            ind_miss = set([i for i, c in enumerate(miss) if c == 'X'])
            if len(ind_miss) > 0:
                ind_in = self.get_inds(seq)
                miss_in_motifs.append(len(ind_in & ind_miss) / len(ind_miss))
                motif_percs.append(len(ind_in)/len(seq))
        return miss_in_motifs, motif_percs

    def get_inds(self, seq):
        lcas = motif_seq.LcSeq(seq, self.k_lca, self.alph_lca, 'lca')
        lces = motif_seq.LcSeq(seq, self.k_lce, self.thresh_lce, 'lce')
        lca_in, lca_out = lcas._get_motif_indexes()
        lce_in, lce_out = lces._get_motif_indexes()
        ind_in = lca_in.union(lce_in)
        return ind_in


def main():
    mm = MissMotif()
    mm.read_df()


if __name__ == '__main__':
    main()