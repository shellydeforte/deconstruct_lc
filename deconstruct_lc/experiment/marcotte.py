"""
Write puncta yes/no, write scores,
make sure all puncta present in set
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import statsmodels.stats.power as smp
from scipy.stats import chi2_contingency

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc.display import display_lc

class Puncta(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.puncta_fp = os.path.join(data_dp, 'experiment', 'marcotte_puncta_proteins.xlsx')
        self.allproteins_fp = os.path.join(data_dp, 'experiment', 'marcotte_proteins.xlsx')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')

    def check_puncta(self):
        puncta = pd.read_excel(self.puncta_fp, sheetname='ST1')
        all = pd.read_excel(self.allproteins_fp, sheetname='Sheet1')
        puncta_orf = list(puncta['ORF'])
        all_orf = list(all['Gene Systematic Name'])
        puncta = []
        all = []
        for item in puncta_orf:
            puncta.append(item[0:7])
        for item in all_orf:
            all.append(item[0:7])
        puncta = set(puncta)
        all = set(all)
        print(len(puncta-all))


def main():
    p = Puncta()
    p.check_puncta()


if __name__ == '__main__':
    main()