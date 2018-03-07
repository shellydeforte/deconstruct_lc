import os
import pandas as pd
from scipy.stats import linregress

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc


class LenNorm(object):
    """
    Note that the linear regression algorithm returns the following:
    m, b, pearsons correlation coefficient, p-value (assuming no slope),
    standard error of the estimated gradient
    """
    def __init__(self, config):
        pdb_dp = os.path.join(config['fps']['data_dp'], 'pdb_prep')
        self.pdb_fp = os.path.join(pdb_dp, 'pdb_norm_cd100.tsv')
        self.seqs = self.read_seqs()
        self.lens = tools_fasta.get_lengths(self.seqs)

    def read_seqs(self):
        df = pd.read_csv(self.pdb_fp, sep='\t', index_col=0)
        seqs = list(df['Sequence'])
        return seqs

    def mb_lca(self, k, lca):
        """LCA"""
        scores = tools_lc.calc_lca_motifs(self.seqs, k, lca)
        lr = linregress(self.lens, scores)
        return lr

    def mb_lce(self, k, lce):
        """LCE"""
        scores = tools_lc.calc_lce_motifs(self.seqs, k, lce)
        lr = linregress(self.lens, scores)
        return lr

    def mb_lc(self, k, lca, lce):
        """LCA || LCE"""
        scores = tools_lc.calc_lc_motifs(self.seqs, k, lca, lce)
        lr = linregress(self.lens, scores)
        return lr

    def mb_lca_and_lce(self, k, lca, lce):
        """LCA & LCE"""
        scores = []
        for seq in self.seqs:
            scores.append(tools_lc.count_lca_and_lce(seq, k, lca, lce))
        lr = linregress(self.lens, scores)
        return lr

    def mb_lca_not_lce(self, k, lca, lce):
        """LCA & ~LCE"""
        scores = []
        for seq in self.seqs:
            scores.append(tools_lc.count_lca_not_lce(seq, k, lca, lce))
        lr = linregress(self.lens, scores)
        return lr

    def mb_not_lca_lce(self, k, lca, lce):
        """~LCA & LCE"""
        scores = []
        for seq in self.seqs:
            scores.append(tools_lc.count_not_lca_lce(seq, k, lca, lce))
        lr = linregress(self.lens, scores)
        return lr


def main():
    config = read_config.read_config()
    k = int(config['score']['k'])
    lca = str(config['score']['lca'])
    lce = float(config['score']['lce'])
    ln = LenNorm(config)
    lr = ln.mb_lc(k, lca, lce)
    print("The slope is {} and the intercept is {}".format(lr[0], lr[1]))


if __name__ == '__main__':
    main()