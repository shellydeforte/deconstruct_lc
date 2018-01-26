import configparser
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class MissScore(object):
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.pdb_an_dp = os.path.join(config['filepaths']['data_dp'],
                                      'pdb_analysis')
        self.an_fpi = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.lc_vs_miss_fp = os.path.join(self.pdb_an_dp, 'lc_vs_miss.tsv')

    def plot_miss_score(self):
        """Plot binned normed scores with avg missing and frac missing"""
        pass

    def write_lc_vs_miss(self):
        df = pd.read_csv(self.an_fpi, sep='\t', index_col=0)
        labels = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30',
                  '30-35', '35-40', '40-45', '45-50', '50+']
        bins = range(0, 50, 5)
        frac_w_miss = []
        num_miss = []
        std_num_miss = []
        for i in bins:
            print(i)
            ndf = df[(df['LC Raw'] >= i) & (df['LCA+LCE'] < i + 5)]
            nm_ndf = ndf[ndf['Miss Count'] > 0]
            frac_w_miss.append(len(nm_ndf)/len(ndf))
            num_miss.append(np.mean(list(nm_ndf['Miss Count'])))
            std_num_miss.append(np.std(list(nm_ndf['Miss Count'])))
        ndf = df[(df['LC Raw'] >= 50)]
        nm_ndf = ndf[ndf['Miss Count'] > 0]
        frac_w_miss.append(len(nm_ndf) / len(ndf))
        num_miss.append(np.mean(list(nm_ndf['Miss Count'])))
        std_num_miss.append(np.std(list(nm_ndf['Miss Count'])))
        df_dict = {'Fraction Missing': frac_w_miss,
                   'Average Missing Residues': num_miss,
                   'STD Missing Residues': std_num_miss,
                   'Labels': labels}
        df_out = pd.DataFrame(df_dict)
        df_out.to_csv(self.lc_vs_miss_fp, sep='\t')


def main():
    ms = MissScore()
    ms.write_lc_vs_miss()


if __name__ == '__main__':
    main()