import os
import pandas as pd

from deconstruct_lc import read_config
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
        ind = ['LCA & ~LCE', '~LCA & LCE', 'LCA & LCE', 'LCA || LCE', 'Total 6-mers']
        df_dict = {'LCA/LCE': ind, 'BC 6-mers': bc_kmers, 'PDB 6-mers': pdb_kmers}
        cols = ['LCA/LCE', 'BC 6-mers', 'PDB 6-mers']
        df = pd.DataFrame(df_dict, columns=cols)
        print(df)

    def get_counts(self, seqs):
        lca_count = self.count_lca_not_lce(seqs)
        lce_count = self.count_not_lca_lce(seqs)
        lca_lce_count = self.count_lca_and_lce(seqs)
        lc_count = self.count_lca_or_lce(seqs)
        total_kmers = self.count_kmers(seqs)
        return [lca_count, lce_count, lca_lce_count, lc_count, total_kmers]

    def count_kmers(self, seqs):
        total_kmers = 0
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            total_kmers += len(kmers)
        return total_kmers

    def count_lca_not_lce(self, seqs):
        total_count = 0
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lca_motif(kmer, self.lca):
                    if not tools_lc.lce_motif(kmer, self.lce):
                        total_count += 1
        return total_count

    def count_not_lca_lce(self, seqs):
        total_count = 0
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, self.lce):
                    if not tools_lc.lca_motif(kmer, self.lca):
                        total_count += 1
        return total_count

    def count_lca_and_lce(self, seqs):
        total_count = 0
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, self.lce):
                    if tools_lc.lca_motif(kmer, self.lca):
                        total_count += 1
        return total_count

    def count_lca_or_lce(self, seqs):
        lc_motifs = tools_lc.calc_lc_motifs(seqs, self.k, self.lca, self.lce)
        return sum(lc_motifs)

    def get_seqs(self, y):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        df = df[df['y'] == y]
        seqs = list(df['Sequence'])
        return seqs


def main():
    at = AmountsTrain()
    at.venn_lc()


if __name__ == '__main__':
    main()