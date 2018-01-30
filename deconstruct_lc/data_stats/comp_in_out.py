import configparser
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc
from deconstruct_lc.motif_seq import LcSeq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class CompStats(object):
    def __init__(self):
        self.train_fpi = os.path.join(config['filepaths']['train_fp'])
        self.k = 6
        self.lca = 'SGEQAPDTNKR'
        self.lce = 1.6

    def comp_lc(self):
        """What is the composition inside LCE motifs?
        Put all LCE motifs into a single string, and do fractions"""
        bc_seqs = self.get_seqs(0)
        pdb_seqs = self.get_seqs(1)
        all_lca_seqs, all_lce_seqs, all_lc_seqs = self.all_lc_seqs(bc_seqs)
        aas = 'SGEQAPDTNKRLHVYFIMCW'
        aas_list = [aa for aa in aas]
        ind = range(len(aas))
        lca_bins = self.get_aa_bins(all_lca_seqs)
        lce_bins = self.get_aa_bins(all_lce_seqs)
        lc_bins = self.get_aa_bins(all_lc_seqs)
        plt.bar(ind, lca_bins, color='darkblue', alpha=0.7, label='LCA')
        plt.bar(ind, lce_bins, color='orange', alpha=0.7, label='LCE')
        plt.bar(ind, lc_bins, color='black', alpha=0.7, label='LC')
        plt.xticks(ind, aas_list)
        plt.legend()
        plt.xlabel('Amino Acids')
        plt.ylabel('Relative Fraction in full dataset')
        plt.show()

    def get_seqs(self, y):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        df = df[df['y'] == y]
        seqs = list(df['Sequence'])
        return seqs

    def get_aa_bins(self, seq):
        aas = 'SGEQAPDTNKRLHVYFIMCW'
        pa = ProteinAnalysis(seq)
        bc_dict = pa.get_amino_acids_percent()
        aa_bins = []
        for aa in aas:
            aa_bins.append(bc_dict[aa])
        return aa_bins

    def all_lc_seqs(self, seqs):
        all_lca_seqs = ''
        all_lce_seqs = ''
        all_lc_seqs = ''
        for seq in seqs:
            lca_seq, lce_seq, lc_seq = self.lc_seqs(seq)
            all_lca_seqs += lca_seq
            all_lce_seqs += lce_seq
            all_lc_seqs += lc_seq
        return all_lca_seqs, all_lce_seqs, all_lc_seqs

    def lc_seqs(self, seq):
        lca, lce, lc = self.get_lc_inds(seq)
        lca_seq = ''
        lce_seq = ''
        lc_seq = ''
        for i, aa in enumerate(seq):
            if i in lca:
                lca_seq += aa
            if i in lce:
                lce_seq += aa
            if i in lc:
                lc_seq += aa
        return lca_seq, lce_seq, lc_seq



    def get_lc_inds(self, seq):
        lce_inds = tools_lc.lce_to_indexes(seq, self.k, self.lce)
        lca_inds = tools_lc.lca_to_indexes(seq, self.k, self.lca)
        lca = lca_inds - lce_inds
        lce = lce_inds - lca_inds
        lc = lca_inds & lce_inds
        return lca, lce, lc


def main():
    cs = CompStats()
    cs.comp_lc()


if __name__ == '__main__':
    main()