import os
import pandas as pd
from Bio import SeqIO
from deconstruct_lc import read_config
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc.analysis_bc.write_bc_score import BcScore
from deconstruct_lc import tools_fasta


class ScoreProfile(object):
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.bc_dp = os.path.join(self.data_dp, 'bc_prep')
        self.bc_an_dp = os.path.join(self.data_dp, 'bc_analysis')
        # Use fasta file with all bc sequences
        self.fasta = os.path.join(self.bc_dp, 'quickgo_bc.fasta')
        self.bc_ss = os.path.join(self.bc_dp, 'quickgo_bc.xlsx')
        self.bc_score_fp = os.path.join(self.bc_an_dp, 'bc_all_score.tsv')

    def open_files(self):
        bc = BcScore()
        bc_names = bc.get_sheets()
        for name in bc_names:
            fn = os.path.join(self.bc_an_dp, '{}_score.tsv'.format(name))
            df_in = pd.read_csv(fn, sep='\t', index_col=0)
            ndf = df_in[df_in['Organism'] == 'YEAST']
            if len(ndf) > 0:
                print(len(ndf))
                print(name)
                sdf = ndf[ndf['LC Score'] < 0]
                print(len(sdf)/len(ndf))
                sdf = ndf[(ndf['LC Score'] >= 0) & (ndf['LC Score'] < 20)]
                print(len(sdf) / len(ndf))
                sdf = ndf[(ndf['LC Score'] >= 20)]
                print(len(sdf) / len(ndf))


def main():
    sp = ScoreProfile()
    sp.open_files()


if __name__ == '__main__':
    main()