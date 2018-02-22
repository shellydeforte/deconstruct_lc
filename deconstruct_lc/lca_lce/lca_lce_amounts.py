import os
import numpy as np
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc.svm import svms
from deconstruct_lc import tools_lc


class AmountsTrain(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.train_fpi = os.path.join(data_dp, 'train.tsv')
        self.k = int(config['score']['k'])
        self.lca = str(config['score']['lca'])
        self.lce = float(config['score']['lce'])

    def venn_lc(self):
        bc_seqs = self.get_seqs(0)
        bc_kmers = self.get_counts(bc_seqs)
        pdb_seqs = self.get_seqs(1)
        pdb_kmers = self.get_counts(pdb_seqs)
        bc_tot_kmers = self.count_kmers(bc_seqs)
        pdb_tot_kmers = self.count_kmers(pdb_seqs)
        ind = ['LCA', 'LCA & ~LCE', 'LCE', '~LCA & LCE', 'LCA & LCE', 'LCA || LCE']
        bc_fracs = [item/bc_tot_kmers for item in bc_kmers]
        pdb_fracs = [item/pdb_tot_kmers for item in pdb_kmers]
        df_dict = {'LCA/LCE': ind, 'BC 6-mers': bc_kmers, 'BC fracs': bc_fracs,
                   'PDB 6-mers': pdb_kmers, 'PDB fracs': pdb_fracs}
        cols = ['LCA/LCE', 'BC 6-mers', 'BC fracs', 'PDB 6-mers', 'PDB fracs']
        df = pd.DataFrame(df_dict, columns=cols)
        print(df)

    def get_counts(self, seqs):
        lca_count = sum(tools_lc.calc_lca_motifs(seqs, self.k, self.lca))
        lca_not_lce_count = sum(self.count_lca_not_lce(seqs))
        lce_count = sum(tools_lc.calc_lce_motifs(seqs, self.k, self.lce))
        lce_not_lca_count = sum(self.count_not_lca_lce(seqs))
        lca_lce_count = sum(self.count_lca_and_lce(seqs))
        lc_count = sum(self.count_lca_or_lce(seqs))
        return [lca_count, lca_not_lce_count, lce_count, lce_not_lca_count, lca_lce_count, lc_count]

    def count_kmers(self, seqs):
        total_kmers = 0
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            total_kmers += len(kmers)
        return total_kmers

    def count_lca_not_lce(self, seqs):
        all_counts = []
        for seq in seqs:
            count = 0
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lca_motif(kmer, self.lca):
                    if not tools_lc.lce_motif(kmer, self.lce):
                        count += 1
            all_counts.append(count)
        return all_counts

    def count_not_lca_lce(self, seqs):
        all_counts = []
        for seq in seqs:
            count = 0
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, self.lce):
                    if not tools_lc.lca_motif(kmer, self.lca):
                        count += 1
            all_counts.append(count)
        return all_counts

    def count_lca_and_lce(self, seqs):
        all_counts = []
        for seq in seqs:
            count = 0
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, self.lce):
                    if tools_lc.lca_motif(kmer, self.lca):
                        count += 1
            all_counts.append(count)
        return all_counts

    def count_lca_or_lce(self, seqs):
        lc_motifs = tools_lc.calc_lc_motifs(seqs, self.k, self.lca, self.lce)
        return lc_motifs

    def get_seqs(self, y):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        df = df[df['y'] == y]
        seqs = list(df['Sequence'])
        return seqs

    def sum_to_svm(self):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        seqs = list(df['Sequence'])
        y = np.array(df['y']).T
        lca_count = tools_lc.calc_lca_motifs(seqs, self.k, self.lca)
        lca_not_lce_count = self.count_lca_not_lce(seqs)
        lce_count = tools_lc.calc_lce_motifs(seqs, self.k, self.lce)
        lce_not_lca_count = self.count_not_lca_lce(seqs)
        lca_lce_count = self.count_lca_and_lce(seqs)
        lc_count = self.count_lca_or_lce(seqs)
        ind = ['LCA', 'LCA & ~LCE', 'LCE', '~LCA & LCE', 'LCA & LCE',
               'LCA || LCE']
        accuracy = [self.run_svm(lca_count, y), self.run_svm(lca_not_lce_count, y),
                    self.run_svm(lce_count, y), self.run_svm(lce_not_lca_count, y),
                    self.run_svm(lca_lce_count, y), self.run_svm(lc_count, y)]
        df_dict = {'LCA/LCE': ind, 'Accuracy': accuracy}
        cols = ['LCA/LCE', 'Accuracy']
        df = pd.DataFrame(df_dict, columns=cols)
        print(df)

    def run_svm(self, motif_sum, y):
        X = np.array([motif_sum]).T
        clf = svms.linear_svc(X, y)
        return clf.score(X, y)


def main():
    at = AmountsTrain()
    at.venn_lc()


if __name__ == '__main__':
    main()