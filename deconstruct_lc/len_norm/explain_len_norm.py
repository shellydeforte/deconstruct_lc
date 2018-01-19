import configparser
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
from deconstruct_lc.params import len_norm


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class ExplainLenNorm(object):
    def __init__(self):
        self.data_dp = config['filepaths']['data_dp']
        self.pdb_norm_fp = os.path.join(self.data_dp, 'scores',
                                      'pdbmiss_6_SGEQAPDTNKR_6_1.6_raw.tsv')
        self.pdb_train_fp = os.path.join(self.data_dp, 'scores',
                                         'pdb_nomiss_cd90_6_SGEQAPDTNKR_6_1.6_norm.tsv')


    def hist_lens(self):
        pdb_train_df = pd.read_csv(self.pdb_train_fp, sep='\t', index_col=0)
        pdb_train_lens = list(pdb_train_df['Length'])
        pdb_train_heights, pdb_train_bins = np.histogram(pdb_train_lens,
                                                       bins=19, range=(100, 2000))
        print("The number of proteins in this pdb train dataset over 700 is "
              "{}".format(sum(pdb_train_heights[5:])))
        pdb_norm_df = pd.read_csv(self.pdb_norm_fp, sep='\t', index_col=0)
        pdb_norm_lens = list(pdb_norm_df['Length'])
        print(max(pdb_norm_lens))
        print(min(pdb_norm_lens))
        pdb_norm_heights, pdb_norm_bins = np.histogram(pdb_norm_lens,
                                                       bins=38, range=(100,
                                                                      2000))
        plt.bar(pdb_norm_bins[:-1],
                pdb_norm_heights,
                width=45,
                color="darkblue",
                alpha=0.7,
                label='PDB')
        print("The number of proteins in this pdb norm dataset over 700 is "
              "{}".format(sum(pdb_norm_heights[5:])))
        print(pdb_norm_heights)

        plt.xlabel('Protein Length')
        plt.ylabel('Relative Fraction')
        plt.legend()
        plt.show()




    def scatter_plots(self):
        pass


def main():
    el = ExplainLenNorm()
    el.hist_lens()


if __name__ == '__main__':
    main()