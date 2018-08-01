import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

from deconstruct_lc import read_config

class MarcotteAnalysis(object):
    """
    Chi square analysis for marcotte data against Huh
    """
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.puncta = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')
        self.nopuncta = os.path.join(data_dp, 'experiment', 'marcotte_nopuncta_scores.tsv')

    def read_files(self):
        """
        All marcotte proteins are >= 100 and <= 2000
        862/886 meet this condition for Huh proteins
        """
        puncta_df = pd.read_csv(self.puncta, sep='\t')
        puncta_short = puncta_df[puncta_df['Length'] <= 250]
        puncta_df = puncta_df[(puncta_df['Length'] >= 100) & (puncta_df['Length'] <= 2000)]
        nopuncta_df = pd.read_csv(self.nopuncta, sep='\t')
        nopuncta_df = nopuncta_df[(nopuncta_df['Length'] >= 100) & (nopuncta_df['Length'] <= 2000)]
        nopuncta_short = nopuncta_df[nopuncta_df['Length'] < 250]
        nopuncta_df['LC Score'].hist(bins=30, range=(-60, 250), normed=True)
        puncta_df['LC Score'].hist(bins=30, range=(-60, 250), normed=True, alpha=0.5)
        plt.show()
        nopuncta_df['Length'].hist(bins=30, range=(100, 2000), normed=True)
        puncta_df['Length'].hist(bins=30, range=(100, 2000), normed=True, alpha=0.5)
        plt.show()
        mlt, mm, mgt = self.get_bins(puncta_df)
        hlt, hm, hgt = self.get_bins(nopuncta_df)
        cont = np.array([[mlt, mm, mgt], [hlt, hm, hgt]])
        print(cont)
        p = chi2_contingency(cont)[1]
        print(p)


    def get_bins(self, df):
        ndf = df[df['LC Score'] < 0]
        lt = len(ndf)
        ndf = df[(df['LC Score'] >= 0) & (df['LC Score'] <= 20)]
        m = len(ndf)
        ndf = df[df['LC Score'] > 20]
        gt = len(ndf)
        return lt, m, gt

def main():
    ma = MarcotteAnalysis()
    ma.read_files()


if __name__ == '__main__':
    main()