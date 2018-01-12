"""
Plot omega, sliding net charge
plot sliding fraction of Q/N, S/T, A/G, P
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from deconstruct_lc import tools_lc
from deconstruct_lc.complementarity import motif_seq

class Fraction(object):

    def __init__(self):
        self.k = 6
        self.lca = 'SGEQAPDTNKR'
        self.lce = 1.6


    def process_seq(self, seq, k):
        kmers = tools_lc.seq_to_kmers(seq, k)
        qn = self.alph_fracs(kmers, 'QN')
        st = self.alph_fracs(kmers, 'ST')
        ag = self.alph_fracs(kmers, 'AG')
        p = self.alph_fracs(kmers, 'P')
        ed = self.alph_fracs(kmers, 'ED')
        kr = self.alph_fracs(kmers, 'KR')
        f = self.alph_fracs(kmers, 'F')
        r = self.alph_fracs(kmers, 'R')
        #plt.plot(qn, label='QN')
        plt.plot(st, label='ST')
        #plt.plot(ag, label='AG')
        #plt.plot(r, label='R')
        #plt.plot(f, label='F')
        lca_x, lca_y, lce_x, lce_y = self.get_motif_index(seq)
        plt.scatter(lca_x, lca_y, color='black', s=2)
        plt.scatter(lce_x, lce_y, color='red', s=2)
        #plt.plot(ed, label='ED')
        #plt.plot(kr, label='KR')
        plt.plot(p, label='P')
        plt.legend()
        plt.show()

    def alph_fracs(self, kmers, alph):
        fracs = []
        for kmer in kmers:
            frac = self.get_frac(kmer, alph)
            fracs.append(frac)
        return fracs

    def get_frac(self, kmer, alph):
        tot_count = 0
        for aa in alph:
            tot_count += kmer.count(aa)
        assert tot_count <= len(kmer)
        frac = float(tot_count)/float(len(kmer))
        return frac

    def get_motif_index(self, sequence):
        mot = motif_seq.LcSeq(sequence, self.k, self.lca, 'lca')
        ind_in, ind_out = mot._get_motif_indexes()
        lca_x = list(ind_in)
        lca_y = [1]*(len(lca_x))
        mot = motif_seq.LcSeq(sequence, self.k, self.lce, 'lce')
        ind_in, ind_out = mot._get_motif_indexes()
        lce_x = list(ind_in)
        lce_y = [1.1]*(len(lce_x))
        return lca_x, lca_y, lce_x, lce_y



class Pipeline(object):

    def __init__(self):
        self.base_fp = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.nmo_fpi = os.path.join(self.base_fp, 'scores',
                                    'nmo_6_SGEQAPDTNKR_6_1.6_seq_scores.tsv')
        self.pdb_fpi = os.path.join(self.base_fp, 'scores',
                                    'pdb_nomiss_cd50_6_SGEQAPDTNKR_6_1.6_seq_scores.tsv')
        self.lca_label = '6_SGEQAPDTNKR'
        self.lce_label = '6_1.6'

    def sandbox(self):
        label = self.lca_label
        df = pd.read_csv(self.nmo_fpi, sep='\t', index_col=0)
        df = df[(df[label] > 30)]
        df = df.sort_values(by=[label])
        df = df.reset_index(drop=True)
        for i, row in df.iterrows():
            sequence = row['Sequence']
            print(len(sequence))
            print(row[label])
            print(row['Protein ID'])
            frac = Fraction()
            frac.process_seq(sequence, 6)



def main():
    pipe = Pipeline()
    pipe.sandbox()


if __name__ == '__main__':
    main()