from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_lc

class AaComp(object):
    """
    Write a continuous string from each type of LC motif to a single fasta
    Use overlapping motifs
    """
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.fdo = os.path.join(data_dp, 'lca_lce')
        self.train_fpi = os.path.join(data_dp, 'train.tsv')
        self.k = int(config['score']['k'])
        self.lca = str(config['score']['lca'])
        self.lce = float(config['score']['lce'])

    def plot_charge(self):
        """
        Result: This is nto working, too much variety
        Plot, in the BC dataset, the KRE fraction in in motif vs. out motif.
        No, wait. No matter how I fucking do this, the fraction will be higher
        because of the reduced alphabet.
        """
        bc_seqs = self.get_seqs(1)
        lca_counts, seq_kmers = self.seq_lca(bc_seqs)
        lca_fracs = self.charge_frac(seq_kmers)
        seq_fracs = self.charge_frac(bc_seqs)
        diff_fracs = []
        for lca_frac, seq_frac in zip(lca_fracs, seq_fracs):
            diff_fracs.append(lca_frac - seq_frac)
        plt.scatter(seq_fracs, lca_fracs, alpha=0.5)
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.show()

    def charge_frac(self, seqs):
        cfracs = []
        for seq in seqs:
            if len(seq) == 0:
                cfracs.append(0)
            else:
                ccount = seq.count('K') + seq.count('E') + seq.count('R')
                cfrac = ccount/len(seq)
                cfracs.append(cfrac)
        return cfracs

    def seq_lca(self, seqs):
        seq_kmers = []
        lca_counts = []
        for seq in seqs:
            lca_motifs = 0
            kmer_str = ''
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lca_motif(kmer, self.lca):
                    kmer_str += kmer
                    lca_motifs += 1
            lca_counts.append(lca_motifs)
            seq_kmers.append(kmer_str)
        return lca_counts, seq_kmers

    def comp_lc(self):
        """What is the composition inside LCE motifs?
        Put all LCE motifs into a single string, and do fractions"""
        bc_seqs = self.get_seqs(0)
        pdb_seqs = self.get_seqs(1)
        bc_seqs = self.get_seqs(0)
        pdb_seqs = self.get_seqs(1)
        bc_lca = self.seq_lca(bc_seqs)
        pdb_lca = self.seq_lca(pdb_seqs)
        #aas = 'SGEQAPDTNKR'
        #aas = 'ERKPSQTGAND'
        aas = 'KRESQPANDGT'
        aas_list = [aa for aa in aas]
        ind = list(range(len(aas)))

        bc_lca_bins = self.get_aa_bins(bc_lca)
        print(bc_lca_bins)
        pdb_lca_bins = self.get_aa_bins(pdb_lca)
        print(pdb_lca_bins)
        plt.bar(ind, pdb_lca_bins, color='darkblue', alpha=0.7, label='PDB')
        plt.bar(ind, bc_lca_bins, color='darkorange', alpha=0.7, label='BC')
        ind = [i + 0.4 for i in ind]
        plt.xticks(ind, aas_list)
        plt.legend()
        plt.xlabel('Amino Acids')
        plt.ylabel('Relative Fraction in LCA motifs')
        plt.ylim([0, 0.16])
        plt.xlim([-1, 12])
        plt.show()

    def get_aa_bins(self, seq):
        #aas = 'SGEQAPDTNKRLHVYFIMCW'
        #aas = 'ERKPSQTGAND'
        aas = 'KRESQPANDGT'
        pa = ProteinAnalysis(seq)
        bc_dict = pa.get_amino_acids_percent()
        aa_bins = []
        for aa in aas:
            aa_bins.append(bc_dict[aa])
        return aa_bins

    def run_seqs(self):
        ind = ['LCA', 'LCA & ~LCE', 'LCE', '~LCA & LCE', 'LCA & LCE',
               'LCA || LCE']
        bc_seqs = self.get_seqs(0)
        pdb_seqs = self.get_seqs(1)
        fpo = os.path.join(self.fdo, 'lca_bc.fasta')
        bc_lca = self.seq_lca(bc_seqs)
        self.write_seqs(bc_lca, fpo)
        pdb_lca = self.seq_lca(pdb_seqs)
        fpo = os.path.join(self.fdo, 'lca_pdb.fasta')
        self.write_seqs(pdb_lca, fpo)

    def write_seqs(self, seq_str, fpo):
        with open(fpo, 'w') as fo:
            fo.write('>\n')
            fo.write(seq_str)

    def get_seqs(self, y):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        df = df[df['y'] == y]
        seqs = list(df['Sequence'])
        return seqs

    def seq_lca_not_lce(self, seqs):
        all_kmers = ''
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lca_motif(kmer, self.lca):
                    if not tools_lc.lce_motif(kmer, self.lce):
                        all_kmers += kmer
        return all_kmers

    def seq_not_lca_lce(self, seqs):
        all_kmers = ''
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, self.lce):
                    if not tools_lc.lca_motif(kmer, self.lca):
                        all_kmers += kmer
        return all_kmers

    def seq_lca_and_lce(self, seqs):
        all_kmers = ''
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, self.lce):
                    if tools_lc.lca_motif(kmer, self.lca):
                        all_kmers += kmer
        return all_kmers

    def seq_lca_or_lce(self, seqs):
        all_kmers = ''
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, self.lce):
                    all_kmers += kmer
                elif tools_lc.lca_motif(kmer, self.lca):
                    all_kmers += kmer
        return all_kmers

    def seq_lca2(self, seqs):
        all_kmers = ''
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lca_motif(kmer, self.lca):
                    all_kmers += kmer
        return all_kmers

    def seq_lce(self, seqs):
        all_kmers = ''
        for seq in seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lce_motif(kmer, self.lce):
                    all_kmers += kmer
        return all_kmers


def main():
    ac = AaComp()
    ac.plot_charge()


if __name__ == '__main__':
    main()