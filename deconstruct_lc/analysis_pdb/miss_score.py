import matplotlib.pyplot as plt
import os
import numpy as np
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
        pos2 = list(range(1, 33, 3))
        x = [i + 0.5 for i in pos2]
        ax1 = fig.add_subplot(211)
        ax2 = ax1.twinx()
        self.frac_miss_box(ax1, ax2, x)
        ax3 = fig.add_subplot(212)
        ax1.set_xlim([0, 33])
        ax2.set_xlim([0, 33])
        ax3.set_xlim([0, 33])
        self.plot_inout_box(ax3, x)
        plt.tight_layout()
        plt.show()

    def frac_miss_box(self, ax1, ax2, x):
        # note labels in boxplots
        frac_miss, miss_counts = self.frac_count_data()
        self.plot_frac(ax1, x, frac_miss)
        self.plot_count(ax2, x, miss_counts)
        ax1.tick_params('y', colors='black')
        ax1.set_ylim([0.75, 1.0])
        ax1.tick_params(axis='x', which='both', labelbottom='off')
        ax2.tick_params(axis='x', which='both', labelbottom='off')
        ax2.set_ylim([0, 260])

    def plot_frac(self, ax, x, frac_miss):
        top_color = 'peru'

        ax.plot(x, frac_miss,
                color='grey')
        ax.scatter(x, frac_miss,
                   marker='o',
                   color=top_color,
                   s=20)
        ax.set_ylabel('Fraction w/ Missing',
                      color=top_color,
                      size=12)

    def plot_count(self, ax, x, miss_counts):
        bot_color = 'maroon'
        bp = {'color': 'grey'}
        wp = {'color': 'grey',
              'linestyle':'-'}
        medianprops = dict(linestyle='-',
                           color='black')
        meanpointprops = dict(marker='D',
                              markeredgecolor=bot_color,
                              markerfacecolor=bot_color,
                              markersize=3)
        ax.boxplot(miss_counts,
                   vert=True,
                   positions=x,
                   whis=[5, 95],
                   widths=1,
                   showfliers=False,
                   showmeans=True,
                   boxprops=bp,
                   whiskerprops=wp,
                   medianprops=medianprops,
                   meanprops=meanpointprops)
        ax.set_ylabel('Mean Missing',
                      color=bot_color,
                      size=12)

    def frac_count_data(self):
        df = pd.read_csv(self.an_fpi, sep='\t', index_col=0)
        bins = range(0, 50, 5)
        miss_counts = []
        frac_miss = []
        for i in bins:
            ndf = df[(df['LC Raw'] >= i) & (df['LC Raw'] < i + 5)]
            nm_ndf = ndf[ndf['Miss Count'] > 0]
            miss_counts.append(list(nm_ndf['Miss Count']))
            frac_miss.append(len(nm_ndf) / len(ndf))
        ndf = df[(df['LC Raw'] >= 50)]
        nm_ndf = ndf[ndf['Miss Count'] > 0]
        miss_counts.append(list(nm_ndf['Miss Count']))
        frac_miss.append(len(nm_ndf) / len(ndf))
        return frac_miss, miss_counts

    def plot_inout_box(self, ax3, x):
        top_color = 'peru'
        bot_color = 'maroon'
        labels = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30',
                  '30-35', '35-40', '40-45', '45-50', '50+']
        motif_perc, miss_in, all_box = self.inout_data()
        pos1 = list(range(2, 34, 3))
        pos2 = list(range(1, 33, 3))
        miss_in_means = []
        for item in miss_in:
            miss_in_means.append(np.mean(item))
        motif_perc_means = []
        for item in motif_perc:
            motif_perc_means.append(np.mean(item))
        bp1 = {'color': 'grey', 'facecolor': 'white'}
        wp1 = {'color': 'grey', 'linestyle':'-', 'alpha': 0.5}
        bp2 = {'color': 'grey', 'alpha': 0.8, 'facecolor': 'grey'}
        wp2 = {'color': 'grey', 'linestyle':'-', 'alpha': 0.5}
        medianprops = dict(linestyle='-', color='black')
        meanpointprops1 = dict(marker='o',
                              markeredgecolor=top_color,
                              markerfacecolor=top_color,
                              markersize=5)
        meanpointprops2 = dict(marker='o',
                              markeredgecolor=bot_color,
                              markerfacecolor=bot_color,
                              markersize=3)

        ax3.plot(pos1, miss_in_means, color=top_color, alpha=0.5)
        ax3.plot(pos2, motif_perc_means, color=bot_color, alpha=0.5)
        bp = ax3.boxplot(miss_in,
                         vert=True,
                         positions=pos1,
                         whis=[5, 95],
                         widths=0.75,
                         showfliers=False,
                         showmeans=True,
                         patch_artist=True,
                         boxprops=bp1,
                         whiskerprops=wp1,
                         medianprops=medianprops,
                         meanprops=meanpointprops1)

        bp2 = ax3.boxplot(motif_perc,
                          vert=True,
                          positions=pos2,
                          whis=[5, 95],
                          widths=1,
                          showfliers=False,
                          showmeans=True,
                          patch_artist=True,
                          boxprops=bp2,
                          whiskerprops=wp2,
                          medianprops=medianprops,
                          meanprops=meanpointprops2)

        ax3.set_ylim([-0.1, 1.1])
        ax3.set_xlim([0, 33])
        ax3.set_ylabel('Missing in LC', color=top_color, size=12)
        ax2 = ax3.twinx()
        ax2.set_ylabel('LC Fraction', color=bot_color, size=12)
        ax3.set_xlabel('LC Motifs')

        ax3.set_xticks(x)
        ax3.set_xticklabels(labels, rotation=45)

    def inout_data(self):
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
        return motif_perc, miss_in, all_box

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