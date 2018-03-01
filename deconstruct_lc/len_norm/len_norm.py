import os
import pandas as pd
from scipy.stats import linregress

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc


class LenNorm(object):
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
        scores = tools_lc.calc_lca_motifs(self.seqs, k, lca)
        lr = linregress(self.lens, scores)
        m = lr[0]
        b = lr[1]
        return m, b

    def mb_lce(self, k, lce):
        scores = tools_lc.calc_lce_motifs(self.seqs, k, lce)
        lr = linregress(self.lens, scores)
        m = lr[0]
        b = lr[1]
        return m, b

    def mb_lc(self, k, lca, lce):
        scores = tools_lc.calc_lc_motifs(self.seqs, k, lca, lce)
        lr = linregress(self.lens, scores)
        m = lr[0]
        b = lr[1]
        return m, b


def main():
    config = read_config.read_config()
    k = int(config['score']['k'])
    lca = str(config['score']['lca'])
    lce = float(config['score']['lce'])
    ln = LenNorm(config)
    m, b = ln.mb_lc(k, lca, lce)
    print("The slope is {} and the intercept is {}".format(m, b))


if __name__ == '__main__':
    main()