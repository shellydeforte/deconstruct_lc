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

    def plot_lc_vs_miss(self):
        df = pd.read_csv(self.lc_vs_miss_fp, sep='\t', index_col=0)
        frac_w_miss = list(df['Fraction Missing'])
        num_miss = list(df['Average Missing Residues'])
        std_num_miss = list(df['STD Missing Residues'])
        labels = list(df['Labels'])
        x = list(range(len(frac_w_miss)))
        plt.xticks(x, labels, rotation=45)
        plt.xlim([-1, len(x)+1])
        #plt.errorbar(x, num_miss, std_num_miss, linestyle='None', marker='o',
        #             capsize=3, label='Average missing residues')
        #plt.ylabel('Average Missing Residues')
        plt.xlabel('LC motifs')
        #plt.ylim([0,200])
        plt.ylabel('Fraction of proteins with missing residues')
        plt.scatter(x, frac_w_miss)
        plt.plot(x, frac_w_miss)
        plt.show()

    def plot_lc_miss(self):
        df = pd.read_csv(self.lc_vs_miss_fp, sep='\t', index_col=0)
        frac_w_miss = list(df['Fraction Missing'])
        num_miss = list(df['Average Missing Residues'])
        std_num_miss = list(df['STD Missing Residues'])
        labels = list(df['Labels'])
        labels = labels
        x = list(range(len(frac_w_miss)))
        fig, ax1 = plt.subplots(sharex=True)
        ax1.plot(x, frac_w_miss, marker='o')
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel('Fraction of proteins with missing residues', color='b')
        ax1.tick_params('y', colors='b')
        plt.xticks(x, labels, rotation=45)
        ax1.set_ylim([0.75, 1.0])
        ax1.set_xlim([-1, len(x)+1])
        ax1.set_xlabel('LC Motifs')

        ax2 = ax1.twinx()
        ax2.errorbar(x, num_miss, std_num_miss, linestyle='None', marker='o',
                     capsize=3, label='Average missing residues', color='r')
        ax2.set_ylabel('Average Missing Residues', color='r')
        ax2.tick_params('y', colors='r')
        ax2.set_ylim([0, 200])
        ax2.set_xlim([-1, len(x)+1])

        fig.tight_layout()
        plt.show()

    def plot_mean(self):
        mean_mm, std_mm, mean_mp, std_mp = self.mean_data()
        x = list(range(len(mean_mm)))
        labels = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30',
                  '30-35', '35-40', '40-45', '45-50', '50+']
        plt.errorbar(x, mean_mm, std_mm, linestyle='None', marker='o',
                     capsize=3, label='Fraction missing residues in LC motif')
        plt.errorbar(x, mean_mp, std_mp, linestyle='None', marker='o',
                     capsize=3, label='Fraction residues in LC motif')
        plt.xticks(x, labels, rotation=45)
        plt.xlim([-1, len(x)+1])
        #plt.ylim([0, 0.8])
        plt.xlabel('LC motifs')
        plt.legend(loc=4)
        plt.show()

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
            ndf = df[(df['LC Raw'] >= i) & (df['LC Raw'] < i + 5)]
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

    def mean_data(self):
        mean_mm = [0.15119716529756219, 0.2758867067395091,
                   0.33919911651251144,
         0.38925749618984801, 0.4596892469792353, 0.45675615911402828,
         0.4864237185593116, 0.47843336509996348, 0.47722958598203197,
         0.52296341132184865, 0.53371100558725326]
        std_mm = [0.267896467804773, 0.31001593805679722,
                    0.29755128257322389,
         0.29214897153214725, 0.29618672624311254, 0.28878338867998538,
         0.27766447616029249, 0.26516401342522217, 0.24012679453077757,
         0.23249365650538631, 0.23073066874878609]
        mean_mp = [0.14288089382642194, 0.19447891989162036,
                   0.2171816720664799,
         0.23594776589707467, 0.25346468713519443, 0.26288893104698952,
         0.27484725570710161, 0.27239470296870616, 0.26238778404020702,
         0.27150317759143594, 0.26612460664234783]
        std_mp = [0.14335880427343892, 0.11564355104930381,
                 0.099416983023802502,
         0.090527165333543019, 0.082859300918348588, 0.083315470100230646,
         0.08419892402540298, 0.077321014349445147, 0.074297419859518155,
         0.064961335129703535, 0.067440855726631221]
        return mean_mm, std_mm, mean_mp, std_mp


def main():
    ms = MissScore()
    ms.plot_lc_miss()


if __name__ == '__main__':
    main()