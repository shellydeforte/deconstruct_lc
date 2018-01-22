import configparser
import os
from scipy.stats import linregress
import pandas as pd
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class LenNorm(object):
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.pdb_fp = os.path.join(self.pdb_dp, 'pdb_norm_cd100.tsv')
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
    ln = LenNorm()
    k = 6
    lca = 'SGEQAPDTNKR'
    m, b = ln.mb_lca(k, lca)
    print(m)
    print(b)


if __name__ == '__main__':
    main()