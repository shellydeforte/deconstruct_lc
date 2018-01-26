import configparser
import os
import pandas as pd

from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class PdbAnalysis(object):
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.all_fpi = os.path.join(self.pdb_dp, 'pdb_all.tsv')
        self.an_fpo = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.k_lca = 6
        self.k_lce = 6
        self.alph_lca = 'SGEQAPDTNKR'
        self.thresh_lce = 1.6
        self.lca_label = '{}_{}'.format(self.k_lca, self.alph_lca)
        self.lce_label = '{}_{}'.format(self.k_lce, self.thresh_lce)

    def write_analysis(self):
        df = pd.read_csv(self.all_fpi, sep='\t')
        print("Size of dataframe before filtering is {}".format(len(df)))
        df = df.drop_duplicates(subset=['Sequence', 'Missing'], keep=False)
        print("Size of dataframe after filtering is {}".format(len(df)))
        df = df.reset_index(drop=True)
        df = self.add_scores(df)
        df.to_csv(self.an_fpo, sep='\t')

    def add_scores(self, df):
        seqs = list(df['Sequence'])
        miss_seqs = list(df['Missing'])
        lcas = tools_lc.calc_lca_motifs(seqs, self.k_lca, self.alph_lca)
        lces = tools_lc.calc_lce_motifs(seqs, self.k_lce, self.thresh_lce)
        lcs = tools_lc.calc_lc_motifs(seqs, self.k_lca, self.alph_lca,
                                      self.thresh_lce)
        lengths = tools_fasta.get_lengths(seqs)
        miss_count = self.get_missing(miss_seqs)
        df['Length'] = lengths
        df['Miss Count'] = miss_count
        df[self.lca_label] = lcas
        df[self.lce_label] = lces
        df['LC'] = lcs
        return df

    def get_missing(self, miss_seqs):
        miss = []
        for seq in miss_seqs:
            miss.append(seq.count('X'))
        return miss


def main():
    pa = PdbAnalysis()
    pa.write_analysis()


if __name__ == '__main__':
    main()