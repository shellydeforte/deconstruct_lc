import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

from deconstruct_lc import read_config

class Expression(object):
    """
    Chi square analysis for marcotte data against Huh
    """
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.exp_dp = os.path.join(data_dp, 'expression_files_for_S3', 'Gasch_2000_PMID_11102521')
        self.exp_fp = os.path.join(self.exp_dp, '2010.Gasch00_stationaryPhase(y14).flt.knn.avg.pcl')
        self.puncta_fpi = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')

    def read_file(self):
        print(self.exp_fp)
        hex = self.hex_low()
        df = pd.read_csv(self.exp_fp, sep='\t')
        df = df[['YORF', 'YPD_2_d_30C; src: t=0<->2_d']]
        puncta_df = pd.read_csv(self.puncta_fpi, sep='\t')
        puncta_low = puncta_df[puncta_df['LC Score'] < -10]
        puncta_hi = puncta_df[puncta_df['LC Score'] >= 20]
        low_orf = list(puncta_low['ORF'])
        hi_orf = list(puncta_hi['ORF'])
        low_df = df[df['YORF'].isin(hex)]
        #hi_df = df[df['YORF'].isin(hi_orf)]
        print(low_df['YPD_2_d_30C; src: t=0<->2_d'].mean())
        #print(hi_df['YPD_2_d_30C; src: t=0<->2_d'].mean())
        print(low_df)

    def hex_low(self):
        hex = ['YDR539W', 'YGR117C', 'YKL035W', 'YER081W', 'YLR028C', 'YDR127W', 'YBL055C',
        'YMR120C', 'YBL039C', 'YDR450W', 'YER175C', 'YJR103W', 'YCL030C', 'YJR057W',
        'YGR210C', 'YER052C', 'YMR303C', 'YLR343W', 'YNL220W', 'YKL001C', 'YGR185C',
        'YMR169C', 'YKL127W', 'YLR344W']
        return hex




def main():
    ex = Expression()
    ex.read_file()


if __name__ == '__main__':
    main()