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

class MissLen(object):
    def __init__(self):
        self.score_dp = os.path.join(config['filepaths']['data_dp'], 'scores')
        self.pdb_score = os.path.join(self.score_dp,
                                      'pdbnorm_6_SGEQAPDTNKR_6_1.6_raw_notunique.tsv')
        self.pdb_df = pd.read_csv(self.pdb_score, sep='\t', index_col=0)
        self.lca_label = '6_SGEQAPDTNKR'
        self.lce_label = '6_1.6'
        self.bin_range = [100, 1900]
        self.bin_size = 50

    def create_bins(self):
        """Display the number of proteins in each length bin"""
        pdb_df = pd.read_csv(self.pdb_score, sep='\t', index_col=0)
        means = []
        stds = []
        x = []
        lengths = list(pdb_df['Length'])
        print(max(lengths))
        print(min(lengths))
        totals = []
        miss_means = []
        lca_means = []
        lce_means = []
        for i in range(self.bin_range[0], self.bin_range[1], self.bin_size):
            df = pdb_df.loc[(pdb_df['Length'] <= i+self.bin_size)
                             & (pdb_df['Length'] > i)]
            totals.append(len(df))
            miss_means.append(np.mean(df['Missing']))
            lca_means.append(np.mean(df['6_SGEQAPDTNKR']))
            lce_means.append(np.mean(df['6_1.6']))
        plt.plot(totals)
        #plt.plot(miss_means)
        #plt.plot(lca_means)
        #plt.plot(lce_means)
        plt.show()

    def len_hist(self):
        df = self.pdb_df[(self.pdb_df['Length'] <= 2000)
                             & (self.pdb_df['Length'] >= 100)]
        df_nomiss = df[df['Missing'] == 0]
        nm_miss_lens = list(df_nomiss['Length'])
        miss_lens = list(df['Length'])
        pdb_heights, pdb_bins = np.histogram(miss_lens, range=(100,
                                                               1000), bins=18)
        npdb_heights, npdb_bins = np.histogram(nm_miss_lens, range=(100,
                                                                    1000),
                                               bins=18)
        print(pdb_heights)
        print(pdb_bins)
        plt.bar(pdb_bins[:-1], pdb_heights, width=50, color='darkblue',
                label="PDB with missing regions")
        plt.bar(npdb_bins[:-1], npdb_heights, width=50, color='orangered',
                label="PDB without missing regions")
        plt.xlabel("Protein Length")
        plt.ylabel("Number of Proteins")
        plt.legend()
        plt.show()

    def avg_missing(self):
        df, df_nm = self.get_df()
        miss_means = []
        stds = []
        x = []
        for i in range(100, 1300, 50):
            x.append(i)
            ndf = df[(df['Length'] >= i) & (df['Length'] < i+50) & df[
                'Missing'] != 0]
            miss_means.append(np.mean(ndf['Missing']))
            stds.append(np.std(ndf['Missing']))
        #plt.scatter(x, miss_means)
        #regr = linear_model.LinearRegression()
        #X = np.array(x).reshape(-1, 1)
        #y = np.array(miss_means).reshape(-1, 1)
        #regr.fit(x, miss_means)
        #m = regr.coef_
        #b = regr.intercept_
        #y_line = self._line_equation(x, m, b)
        #plt.plot(x, y_line)
        plt.errorbar(x, miss_means, stds, linestyle='None', marker='o',
                     color='darkblue')
        plt.ylim([0, 300])
        plt.xlabel("Protein Length")
        plt.ylabel("Average missing residues")
        plt.show()

    def len_vs_miss(self, df):
        miss_means = []
        stds = []
        x = []
        totals = []
        for i in range(100, 1000, 50):
            x.append(i)
            ndf = df[(df['Length'] >= i) & (df['Length'] < i+50) & df[
                'Missing'] != 0]
            totals.append(len(ndf))
            miss_means.append(np.mean(ndf['Missing']))
            stds.append(np.std(ndf['Missing']))
        print(totals)
        plt.errorbar(x, miss_means, stds, linestyle='None', marker='o',
                     color='darkblue')
        plt.ylim([0, 300])
        plt.xlabel("Protein Length")
        plt.ylabel("Average missing residues")
        plt.show()

    def get_df(self):
        df = self.pdb_df[(self.pdb_df['Length'] <= 2000)
                             & (self.pdb_df['Length'] >= 100)]
        df_nm = df[df['Missing'] == 0]
        return df, df_nm

    def _line_equation(self, x, m, b):
        y = m*x + b
        return y

    def score_miss_len(self):
        df = self.pdb_df[(self.pdb_df['Length'] <= 2000)
                         & (self.pdb_df['Length'] >= 100)]
        df = df[(df['6_SGEQAPDTNKR'] <= 40) & (df['6_SGEQAPDTNKR'] >= 30)]
        return df


def main():
    ml = MissLen()
    df = ml.score_miss_len()
    ml.len_vs_miss(df)


if __name__ == '__main__':
    main()