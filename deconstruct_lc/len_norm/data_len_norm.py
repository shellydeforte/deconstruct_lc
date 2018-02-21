import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_lc


class DataLenNorm(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        pdb_dp = os.path.join(data_dp, 'pdb_prep')
        pdb_an_dp = os.path.join(data_dp, 'pdb_analysis')
        self.norm_fpi = os.path.join(pdb_dp, 'pdb_norm_cd100.tsv')
        self.fpo = os.path.join(pdb_an_dp, 'pdb_len_norm.tsv')
        self.k = int(config['score']['k'])
        self.lca = str(config['score']['lca'])
        self.lce = float(config['score']['lce'])

    def write_tsv(self):
        """Write a tsv file that is score, nomiss_score, length"""
        df = pd.read_csv(self.norm_fpi, sep='\t', index_col=0)
        seqs = df['Sequence']
        miss_seqs = df['Missing']
        lens = [len(seq) for seq in seqs]
        raw_scores = tools_lc.calc_lc_motifs(seqs, self.k, self.lca, self.lce)
        nomiss_scores = tools_lc.calc_lc_motifs_nomiss(seqs, miss_seqs, self.k,
                                                       self.lca, self.lce)
        df_dict = {'score': raw_scores, 'nomiss_score': nomiss_scores, 'Length': lens}
        cols = ['score', 'nomiss_score', 'Length']
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.fpo, sep='\t')


def main():
    dl = DataLenNorm()
    dl.write_tsv()


if __name__ == '__main__':
    main()