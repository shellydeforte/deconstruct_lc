import os
import pandas as pd

from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc
from deconstruct_lc.scores.norm_score import NormScore


class PdbAnalysis(object):
    def __init__(self, config):
        data_dp = config['fps']['data_dp']
        pdb_dp = os.path.join(data_dp, 'data_pdb')
        self.all_fpi = os.path.join(pdb_dp, 'pdb_all.tsv')
        self.an_fpo = os.path.join(pdb_dp, 'pdb_analysis.tsv')
        self.k = config['score'].getint('k')
        self.alph_lca = config['score'].get('lca')
        self.thresh_lce = config['score'].getfloat('lce')

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
        ns = NormScore()
        lc_raw = tools_lc.calc_lc_motifs(seqs, self.k, self.alph_lca, self.thresh_lce)
        lc_norms = ns.lc_norm_score(seqs)
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