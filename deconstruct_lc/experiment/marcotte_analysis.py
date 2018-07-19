import os
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

from deconstruct_lc import read_config

class MarcotteAnalysis(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.marcotte = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')
        self.huh = os.path.join(data_dp, 'experiment', 'huh_scores.tsv')

    def read_files(self):
        """
        All marcotte proteins are >= 100 and <= 2000
        862/886 meet this condition for Huh proteins
        """
        marc_df = pd.read_csv(self.marcotte, sep='\t')
        marc_df = marc_df[(marc_df['Length'] >= 100) & (marc_df['Length'] <= 2000)]
        huh_df = pd.read_csv(self.huh, sep='\t')
        huh_df = huh_df[(huh_df['Length'] >= 100) & (huh_df['Length'] <= 2000)]
        huh_df = huh_df[huh_df['ORF']]
        mlt, mm, mgt = self.get_bins(marc_df)
        hlt, hm, hgt = self.get_bins(huh_df)
        cont = np.array([[mlt, mm, mgt], [hlt, hm, hgt]]).T
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