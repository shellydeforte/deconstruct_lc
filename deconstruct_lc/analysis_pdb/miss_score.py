import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import motif_seq

class MissScore(object):
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.pdb_dp = os.path.join(self.data_dp, 'pdb_prep')
        self.pdb_an_dp = os.path.join(self.data_dp,
                                      'pdb_analysis')
        self.an_fpi = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.miss_fp = os.path.join(self.pdb_an_dp, 'miss_in_out.tsv')
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'

    def plot_all(self):
        """
        subplot(nrows, ncolumns, index)
        """
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = ax1.twinx()
        self.frac_miss_box(ax1, ax2)
        #ax3 = fig.add_subplot(212, sharex=ax1)
        ax3 = fig.add_subplot(212)
        ax1.set_xlim([0, 12])
        ax2.set_xlim([0, 12])
        ax3.set_xlim([0, 12])
        self.plot_inout_box(ax3)
        #plt.tight_layout()
        plt.show()

    def frac_miss_box(self, ax1, ax2):
        # note labels in boxplots
        df = pd.read_csv(self.an_fpi, sep='\t', index_col=0)
        bins = range(0, 50, 5)
        miss_counts = []
        frac_miss = []
        bp = {'color': 'black'}
        wp = {'color': 'black', 'linestyle':'-'}
        medianprops = dict(linestyle='-', color='black')
        meanpointprops = dict(marker='D', markeredgecolor='black',
                              markerfacecolor='black', markersize=2)
        for i in bins:
            ndf = df[(df['LC Raw'] >= i) & (df['LC Raw'] < i + 5)]
            nm_ndf = ndf[ndf['Miss Count'] > 0]
            miss_counts.append(list(nm_ndf['Miss Count']))
            frac_miss.append(len(nm_ndf)/len(ndf))
        ndf = df[(df['LC Raw'] >= 50)]
        nm_ndf = ndf[ndf['Miss Count'] > 0]
        miss_counts.append(list(nm_ndf['Miss Count']))
        frac_miss.append(len(nm_ndf) / len(ndf))
        x = list(range(1, len(frac_miss)+1))
        #x = list(range(0, 22, 2))
        print(len(x))
        print(len(frac_miss))
        ax1.plot(x, frac_miss, color='black')
        ax1.scatter(x, frac_miss, marker='o', color='green', s=70)
        ax1.set_ylabel('Fraction w/ missing', color='darkgreen', size=12)
        ax1.tick_params('y', colors='black')
        ax1.set_ylim([0.75, 1.0])
        ax1.tick_params(axis='x', which='both', labelbottom='off')
        ax2.tick_params(axis='x', which='both', labelbottom='off')
        ax2.boxplot(miss_counts, vert=True, whis=[5, 95], widths=0.5,
                   boxprops=bp, whiskerprops=wp, showfliers=False,
                    showmeans=True, medianprops=medianprops,
                    meanprops=meanpointprops)
        ax2.set_ylim([0, 260])
        ax2.set_ylabel('Missing residues', color='black')

    def plot_mean(self, ax3):
        mean_mm, std_mm, mean_mp, std_mp = self.mean_data()
        x = list(range(1, len(mean_mm)+1))
        labels = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30',
                  '30-35', '35-40', '40-45', '45-50', '50+']
        ax3.errorbar(x, mean_mm, std_mm, linestyle='None', marker='o',
                     capsize=4, lw=2, label='Fraction missing residues in LC '
                                      'motif', color='black', markersize=8)
        ax3.errorbar(x, mean_mp, std_mp, linestyle='None', marker='o',
                     capsize=5, lw=3, label='Fraction residues in LC motif',
                     color='grey', markersize=8)
        #ax3.set_xticks(x, labels, rotation=45)
        #ax3.set_xlim([-1, len(x)+1])
        ax3.set_ylim([0, 0.8])
        ax3.set_xlabel('LC motifs', size=12)
        #plt.legend(loc=4)
        #plt.legend(bbox_to_anchor=(1.017, 1.14))
        #plt.tight_layout()
        #plt.show()

    def plot_inout_box(self, ax3):
        labels = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30',
                  '30-35', '35-40', '40-45', '45-50', '50+']
        df = pd.read_csv(self.miss_fp, sep='\t', index_col=0)
        bins = range(0, 50, 5)
        miss_in = []
        motif_perc = []
        all_box = []
        for i in bins:
            ndf = df[(df['LC Raw'] >= i) & (df['LC Raw'] < i + 5)]
            all_box.append(list(ndf['Motif perc']))
            all_box.append(list(ndf['Miss in motif']))

            miss_in.append(list(ndf['Miss in motif']))
            motif_perc.append(list(ndf['Motif perc']))
        ndf = df[(df['LC Raw'] >= 50)]
        miss_in.append(list(ndf['Miss in motif']))
        motif_perc.append(list(ndf['Motif perc']))
        bp1 = {'color': 'grey', 'facecolor': 'white'}
        wp1 = {'color': 'black', 'linestyle':'-'}
        bp2 = {'color': 'black'}
        wp2 = {'color': 'black', 'linestyle':'-'}
        medianprops = dict(linestyle='-', color='black')
        bp = ax3.boxplot(all_box, vert=True, whis=[5, 95],
                         widths=0.5, boxprops=bp1, whiskerprops=wp1,
                         showfliers=False, patch_artist=True,
                         medianprops=medianprops)
        #plt.setp(bp['boxes'], color='grey')
        #ax3.boxplot(miss_in, vert=True, whis=[5, 95], widths=0.5,
        #            boxprops=bp2, whiskerprops=wp2, showfliers=False)
        ax3.set_ylim([-0.5, 1.5])
        x = list(range(2, 22, 2))
        plt.xticks(x, labels, rotation=45)
        plt.show()

    def write_in_motif(self):
        df = pd.read_csv(self.an_fpi, sep='\t', index_col=0)
        miss_in_motifs, motif_percs, lc_raw = self.lc_blobs(df)
        df_dict = {'Miss in motif': miss_in_motifs, 'Motif perc':
            motif_percs, 'LC Raw': lc_raw}
        df_out = pd.DataFrame(df_dict)
        df_out.to_csv(self.miss_fp, sep='\t')

    def lc_blobs(self, df):
        miss_in_motifs = []
        motif_percs = []
        lc_raw = []
        for i, row in df.iterrows():
            print(i)
            miss = row['Missing']
            seq = row['Sequence']
            ind_miss = set([i for i, c in enumerate(miss) if c == 'X'])
            if len(ind_miss) > 0:
                ind_in = self.get_inds(seq)
                miss_in_motifs.append(len(ind_in & ind_miss) / len(ind_miss))
                motif_percs.append(len(ind_in)/len(seq))
                lc_raw.append(row['LC Raw'])
        return miss_in_motifs, motif_percs, lc_raw

    def get_inds(self, seq):
        lcas = motif_seq.LcSeq(seq, self.k, self.lca, 'lca')
        lces = motif_seq.LcSeq(seq, self.k, self.lce, 'lce')
        lca_in, lca_out = lcas._get_motif_indexes()
        lce_in, lce_out = lces._get_motif_indexes()
        ind_in = lca_in.union(lce_in)
        return ind_in

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
    ms.plot_all()


if __name__ == '__main__':
    main()