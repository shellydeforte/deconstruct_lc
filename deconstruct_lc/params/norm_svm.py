import os
import numpy as np
import pandas as pd
from scipy.stats import linregress

from deconstruct_lc import read_config
from deconstruct_lc.svm import svms
from deconstruct_lc import tools_lc


class NormSvm(object):
    def __init__(self, config):
        self.config = config
        data_dp = config['fps']['data_dp']
        pdb_dp = os.path.join(config['fps']['data_dp'], 'pdb_prep')
        pdb_fp = os.path.join(pdb_dp, 'pdb_norm_cd100.tsv')
        pdb_df = pd.read_csv(pdb_fp, sep='\t', index_col=0)
        self.pdb_seqs = list(pdb_df['Sequence'])
        self.pdb_lens = [len(pdb_seq) for pdb_seq in self.pdb_seqs]

        self.param_dp = os.path.join(data_dp, 'params')
        self.lce_lab_fp = os.path.join(self.param_dp, 'top_svm_lce.tsv')
        self.lca_lab_fp = os.path.join(self.param_dp, 'rep_lca.txt')

        train_fp = os.path.join(data_dp, 'train.tsv')
        train_df = pd.read_csv(train_fp, sep='\t', index_col=0)
        self.y = list(train_df['y'])
        self.seqs = list(train_df['Sequence'])

    def read_labels(self):
        lce_df = pd.read_csv(self.lce_lab_fp, sep='\t')
        lce_labs = list(lce_df['Label'])
        lca_labs = []
        with open(self.lca_lab_fp, 'r') as fpi:
            for line in fpi:
                lca_labs.append(line.strip())
        combos = []
        for lce in lce_labs:
            for lca in lca_labs:
                combos.append((lce, lca))
        return combos

    def run_svm(self):
        """
        if k is not the same
        """
        lce_lca = self.read_labels()
        for lce_lab, lca_lab in lce_lca:
            k_lce = int(lce_lab.split('_')[0])
            k_lca = int(lca_lab.split('_')[0])
            lce = float(lce_lab.split('_')[1])
            lca = str(lca_lab.split('_')[1])

    def lca_lce_vec(self, k_lca, k_lce, lca, lce):
        """for all check lca, lce, and [lca, lce]"""
        lca_raw = tools_lc.calc_lca_motifs(self.seqs, k_lca, lca)

    def norm_function(self, m, b, raw_scores, lengths):
        norm_scores = []
        for raw_score, length in zip(raw_scores, lengths):
            norm_score = raw_score - ((m * length) + b)
            norm_scores.append(norm_score)
        return norm_scores

    def count_lca_not_lce(self, seqs, k, lca, lce):
        all_counts = []
        for seq in seqs:
            count = 0
            kmers = tools_lc.seq_to_kmers(seq, k)
            for kmer in kmers:
                if tools_lc.lca_motif(kmer, lca):
                    if not tools_lc.lce_motif(kmer, lce):
                        count += 1
            all_counts.append(count)
        return all_counts

    def count_not_lca_lce(self, seqs, k, lca, lce):
        all_counts = []
        for seq in seqs:
            count = 0
            kmers = tools_lc.seq_to_kmers(seq, k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, lce):
                    if not tools_lc.lca_motif(kmer, lca):
                        count += 1
            all_counts.append(count)
        return all_counts

    def count_lca_and_lce(self, seqs, k, lca, lce):
        all_counts = []
        for seq in seqs:
            count = 0
            kmers = tools_lc.seq_to_kmers(seq, k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, lce):
                    if tools_lc.lca_motif(kmer, lca):
                        count += 1
            all_counts.append(count)
        return all_counts

    def count_lca_or_lce(self, seqs, k, lca, lce):
        lc_motifs = tools_lc.calc_lc_motifs(seqs, k, lca, lce)
        return lc_motifs

    def get_mb(self, ln_scores):
        lr = linregress(self.pdb_lens, ln_scores)
        m = lr[0]
        b = lr[1]
        return m, b


def main():
    config = read_config.read_config()
    ns = NormSvm(config)
    ns.run_svm()


if __name__ == '__main__':
    main()