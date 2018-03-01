import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc
from deconstruct_lc.scores.norm_score import NormScore


class PdbAnalysis(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.pdb_dp = os.path.join(data_dp, 'pdb_prep')
        self.all_fpi = os.path.join(self.pdb_dp, 'pdb_all.tsv')
        self.an_fpo = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.k_lca = 6
        self.k_lce = 6
        self.alph_lca = 'SGEQAPDTNKR'
        self.thresh_lce = 1.6

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
        ns = NormScore(seqs)
        lc_raw = tools_lc.calc_lc_motifs(seqs, self.k_lca, self.alph_lca,
                                         self.thresh_lce)
        lc_norms = ns.lc_norm_score()
        lengths = tools_fasta.get_lengths(seqs)
        miss_count = self.get_missing(miss_seqs)
        df['Length'] = lengths
        df['Miss Count'] = miss_count
        df['LC Norm'] = lc_norms
        df['LC Raw'] = lc_raw
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