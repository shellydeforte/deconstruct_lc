from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import localcider
from localcider.sequenceParameters import SequenceParameters

from deconstruct_lc import read_config
from deconstruct_lc import tools_lc
from deconstruct_lc.svm import svms
from deconstruct_lc import motif_seq

class LcaSvmComp(object):
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
        for seq in seqs:
            ms = motif_seq.LcSeq(seq, self.k, self.lca, 'lca')
            in_seq, out_seq = ms.seq_in_motif()
            SeqOb = SequenceParameters(in_seq)
            print(SeqOb.get_kappa())
            seqOb = SequenceParameters(out_seq)
            print(seqOb.get_kappa())
            print('')

    def check_one_charge(self):
        """
        Result. If you remove K, R, E, your classification accuracy goes to 0.71
        Hypothesis: it is the LCAs with K/R/E that matter the most for
        classification. So what if we only count LCAs with a charged residue?
        """
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        #df = df[df['y'] == 0]
        seqs = list(df['Sequence'])
        lca_counts = self.count_lca_charge(seqs)
        #plt.hist(lca_counts, bins=20, range=(0, 70))
        #plt.ylim([0, 900])
        #plt.show()
        X = np.array([lca_counts]).T
        y = np.array(df['y']).T
        clf = svms.linear_svc(X, y)
        print(clf.score(X, y))

    def count_lca_charge(self, seqs):
        lca_counts = []
        for seq in seqs:
            lca_motifs = 0
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                if tools_lc.lca_motif(kmer, self.lca):
                    if not tools_lc.lce_motif(kmer, self.lce):
                        if ('K' in kmer) and ('R' in kmer) and ('E' in kmer):
                            lca_motifs += 1
            lca_counts.append(lca_motifs)
        return lca_counts

    def create_feature_vecs(self):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        seqs = list(df['Sequence'])
        y = list(df['y'])
        self.feat_vec(seqs, y)

    def feat_vec(self, seqs, y):
        """
        For each sequence, create a feature vector that is # motifs,
        and fraction for each amino acid of the LCA
        """
        lca_counts, seq_kmers = self.seq_lca(seqs)
        df_dict = {'seq_kmer': seq_kmers, 'lca_count': lca_counts, 'y': y}
        df = pd.DataFrame(df_dict)
        print(len(df))
        ndf = df
        #ndf = df[(df['lca_count'] > 20) & (df['lca_count'] < 30)]
        print(len(ndf[ndf['y'] == 0]))
        print(len(ndf[ndf['y'] == 1]))
        y = np.array(ndf['y']).T
        xs = []
        for i, row in ndf.iterrows():
            feat_vec = self.one_feat_vec(str(row['seq_kmer']))
            #feat_vec = []
            #feat_vec.append(int(row['lca_count']))
            xs.append(feat_vec)
        X = np.array(xs)
        clf = svms.normal_rbf(X, y)
        print(clf.score(X, y))

    def one_feat_vec(self, seq_kmer):
        aas = 'KRESQPANDGT'
        feat_vec = []
        for aa in aas:
            if seq_kmer.count(aa) > 0:
                feat_vec.append(seq_kmer.count(aa)/len(seq_kmer))
            else:
                feat_vec.append(0)
        return feat_vec

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


def main():
    ls = LcaSvmComp()
    ls.in_out_kappa()


if __name__ == '__main__':
    main()