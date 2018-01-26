import configparser
import os
from Bio import SeqIO
import random
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import numpy as np
import pandas as pd
from deconstruct_lc import motif_seq
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class LcLens(object):
    """
    What is the relationship between the lc score and the longest continuous
    length?
    """
    def __init__(self):
        self.dp = os.path.join(config['filepaths']['data_dp'])
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.pdb_an_dp = os.path.join(config['filepaths']['data_dp'],
                                      'pdb_analysis')
        self.pdb_an_fp = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.train_fp = os.path.join(self.dp, 'train.tsv')
        self.k_lca = 6
        self.k_lce = 6
        self.alph_lca = 'SGEQAPDTNKR'
        self.thresh_lce = 1.6
        self.lca_label = '{}_{}'.format(self.k_lca, self.alph_lca)
        self.lce_label = '{}_{}'.format(self.k_lce, self.thresh_lce)

    def lens_vs_lc(self):
        df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        ndf = df[(df['LC'] >= 30)]
        seqs = ndf['Sequence']
        all_lens = []
        for seq in seqs:
            lens = tools_lc.lc_to_lens(seq, self.k_lca, self.alph_lca,
                                       self.thresh_lce)
            if len(lens) > 0:
                all_lens.append(max(lens))
        print(np.mean(all_lens))
        plt.hist(all_lens, bins=50)
        plt.show()

    def bc_pdb_lens(self):
        df = pd.read_csv(self.train_fp, sep='\t', index_col=0)
        bc_df = df[df['y'] == 0]
        pdb_df = df[df['y'] == 1]
        pdb_df = pdb_df[pdb_df['LC'] > 30]
        bc_seqs = bc_df['Sequence']
        pdb_seqs = pdb_df['Sequence']
        all_bc_lens = []
        all_pdb_lens = []
        for bc_seq in bc_seqs:
            bc_lens = tools_lc.lc_to_lens(bc_seq, self.k_lca, self.alph_lca,
                                       self.thresh_lce)
            if len(bc_lens) > 0:
                all_bc_lens.append(max(bc_lens))
        for pdb_seq in pdb_seqs:
            pdb_lens = tools_lc.lc_to_lens(pdb_seq, self.k_lca,
                                           self.alph_lca, self.thresh_lce)
            if len(pdb_lens) > 0:
                all_pdb_lens.append(max(pdb_lens))
        #plt.hist(all_bc_lens, bins=50)
        plt.hist(all_pdb_lens, bins=50)
        plt.show()


def main():
    ll = LcLens()
    ll.lens_vs_lc()


if __name__ == '__main__':
    main()