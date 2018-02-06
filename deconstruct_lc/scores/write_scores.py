import configparser
import os
import pandas as pd
from deconstruct_lc.scores.norm_score import NormScore
import numpy as np


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class WriteNorm(object):
    """Write tsv file with column header proteome: yeast, human, BC PDB
    LC score"""
    def __init__(self):
        data = config['filepaths']['data_dp']
        self.train_fpi = os.path.join(data, 'train.tsv')
        self.train_scores = os.path.join(data, 'scores', 'train_scores.tsv')

    def write_train(self):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        seqs = df['Sequence']
        y = df['y']
        ns = NormScore(seqs)
        scores = ns.lc_norm_score()
        cols = ['y', 'LC score']
        df_dict = {'y': y, 'LC score': scores}
        df_out = pd.DataFrame(df_dict, columns=cols)
        df_out.to_csv(self.train_scores, sep='\t')


def main():
    wn = WriteNorm()
    wn.write_train()


if __name__ == '__main__':
    main()