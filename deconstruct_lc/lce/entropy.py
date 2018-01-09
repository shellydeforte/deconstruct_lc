import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from predict_llps import tools_lc
from predict_llps import tools_fasta


class CalcEntropy(object):
    def __init__(self):
        self.fd = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.cb = os.path.join(self.fd, 'quickgo', 'quickgo_cb.fasta')
        self.pdb = os.path.join(self.fd, 'pdb_nomiss_cd50.fasta')

    def find_entropy_values(self, k):
        """
        Return all the possible shannon entropies in my data set for the given
        k-mer length, rounded up to the nearest 0.1.
        """
        all_seqs = tools_fasta.fasta_to_seq(self.cb)
        all_shannon = set()
        for seq in all_seqs:
            kmers = tools_lc.seq_to_kmers(seq, k)
            for kmer in kmers:
                s = tools_lc.shannon(kmer)
                all_shannon.add(s)
        new_scores = []
        all_shannon = sorted(list(all_shannon), reverse=True)
        print(sorted(all_shannon))
        for score in all_shannon[1:]:
            new_score = self.round_up(score)
            new_scores.append(new_score)
            new_scores = sorted(list(set(new_scores)), reverse=True)
        return new_scores

    def round_up(self, score):
        """
        Given a float, round up to the nearest 10th, and then format into
        string.
        """
        new_score = round(score, 1)
        if new_score < score:
            new_score += 0.1
            new_score = round(new_score, 1)  # necessary because trailing 0.999..
        return new_score

    def plot_entropy(self):
        k = 6
        cb, pdb = self.load_entropy()
        plt.hist(pdb, 10, normed=1, edgecolor='darkblue', alpha=0.5,
                                    cumulative=True, histtype='step', lw=3)
        plt.hist(cb, 10, normed=1, edgecolor='orangered', alpha=0.5,
                                    cumulative=True, histtype='step', lw=3)
        ent_values = [2.6, 2.3, 2.0, 1.8, 1.6, 1.5, 1.3, 1.0, 0.7, 0.0]
        pdb_am = [2016, 2016, 2016, 2003, 1801, 1724, 929, 230, 118, 0]
        cb_am = [3436, 3436, 3436, 3428, 3390, 3378, 2950, 1935, 1416, 0]
        pdb_n = []
        cb_n = []
        for item in pdb_am:
            pdb_n.append(item/pdb_am[0])
        for item in cb_am:
            cb_n.append(item/cb_am[0])
        plt.scatter(ent_values, cb_n, c='orangered', alpha=0.5, marker='s',
                    label='BC proteins w/ motif')
        plt.scatter(ent_values, pdb_n, c='darkblue', alpha=0.5, marker='s',
                    label='PDB proteins w/ motif')
        plt.xlabel('Shannon Information Entropy')
        plt.ylabel('Cumulative Fraction of 6-mers below LCE threshold')
        plt.title('LCE Thresholds for k = {}'.format(str(k)))
        plt.legend()
        plt.show()

    def write_entropy(self):
        k = 6
        cb = self.count_entropy(k, self.cb)
        pdb = self.count_entropy(k, self.pdb)
        df_cb = pd.DataFrame({'cb': cb})
        df_pdb = pd.DataFrame({'pdb': pdb})
        df_cb.to_csv('cb.tsv', sep='\t')
        df_pdb.to_csv('pdb.tsv', sep='\t')

    def count_entropy(self, k, fasta_fp):
        """
        Place entropy values into bins
        """
        all_seqs = tools_fasta.fasta_to_seq(fasta_fp)
        all_shannon = []
        for seq in all_seqs:
            kmers = tools_lc.seq_to_kmers(seq, k)
            for kmer in kmers:
                s = tools_lc.shannon(kmer)
                if s == 2.584962500721156: # This is how I filled in the
                    # entropy table
                    print(kmer)
                all_shannon.append(s)
        return all_shannon

    def load_entropy(self):
        df_cb = pd.read_csv('cb.tsv', sep='\t', index_col=0)
        cb = list(df_cb['cb'])
        df_pdb = pd.read_csv('pdb.tsv', sep='\t', index_col=0)
        pdb = list(df_pdb['pdb'])
        return cb, pdb

    def count_perc_w_le(self):
        """
        pdb = [2016, 2016, 2016, 2003, 1801, 1724, 929, 230, 118, 0]
        cb = [3436, 3436, 3436, 3428, 3390, 3378, 2950, 1935, 1416, 0]
        """
        k = 6
        ent_values = [2.6, 2.3, 2.0, 1.8, 1.6, 1.5, 1.3, 1.0, 0.7, 0.0]
        all_seqs = tools_fasta.fasta_to_seq(self.cb)
        all_ents = []
        for ent in ent_values:
            has_ent = 0
            for seq in all_seqs:
                kmers = tools_lc.seq_to_kmers(seq, k)
                for kmer in kmers:
                    s = tools_lc.shannon(kmer)
                    if s < ent:
                        has_ent += 1
                        break
            all_ents.append(has_ent)
        print(all_ents)


def main():
    ce = CalcEntropy()
    ce.plot_entropy()


if __name__ == '__main__':
    main()
