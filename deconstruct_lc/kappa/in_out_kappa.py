from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from deconstruct_lc.kappa import kappa
from deconstruct_lc import read_config
from deconstruct_lc import tools_lc
from deconstruct_lc.svm import svms
from deconstruct_lc import motif_seq

class InOutKappa(object):
    """
    Results: the composition is not enough to tell apart these regions
    I have made the observation that certain residues, particularly charged
    residues are more highly represented in LCA motifs in BC vs. PDB.

    There is a certain number of BC proteins that are below the score lines of
    20, and 0. Here are my questions:
    Of these proteins, do any have 0 motifs?
    For those that have > 0 motifs, can we compare the amino acid composition
    within the motifs to the amino acid composition within PDB motifs?
    Does that help us classify within this scoring region?
    """
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.fdo = os.path.join(data_dp, 'lca_lce')
        self.train_fpi = os.path.join(data_dp, 'train.tsv')
        self.k = int(config['score']['k'])
        self.lca = str(config['score']['lca'])
        self.lce = float(config['score']['lce'])

    def in_out_kappa(self):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        df = df[df['y'] == 0]
        seqs = list(df['Sequence'])
        all_deltas = []
        net_charges = []
        frac_charges = []
        for seq in seqs:
            ms = motif_seq.LcSeq(seq, self.k, self.lca, 'lca')
            in_seq, out_seq = ms.seq_in_motif()
            in_kmer, out_kmer = ms.overlapping_kmer_in_motif()
            if len(in_kmer) > 20:
                ka = kappa.KappaKmers(out_kmer, out_seq)
                if ka.FCR() > 0.1:
                    delta = ka.deltaForm()
                    net_charges.append(ka.NCPR())
                    print(out_seq)
                    print(delta)
                    all_deltas.append(delta)
                    frac_charges.append(ka.FCR())
        #plt.hist(net_charges)
        plt.scatter(net_charges, all_deltas, alpha=0.5, color='grey')
        #plt.ylim([0, 0.35])
        plt.ylim([0, 0.5])
        plt.xlim([-0.8, 0.8])
        #plt.xlim([0, 0.4])
        plt.xlabel('Net charge per residue', size=14)
        plt.ylabel('Charge Asymmetry (Delta)', size=14)
        plt.title('Outside LC Motifs')
        plt.show()


def main():
    ls = InOutKappa()
    ls.in_out_kappa()


if __name__ == '__main__':
    main()