import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import pandas as pd
from deconstruct_lc import read_config
import numpy as np


class PlotScores(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.puncta = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')
        self.nopuncta = os.path.join(data_dp, 'experiment', 'marcotte_nopuncta_scores.tsv')

    def plot_bg(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.add_patch(patches.Rectangle((-30, 0), 30, 5, facecolor='grey'))
        ax.add_patch(patches.Rectangle((0, 0), 20, 5, facecolor='darkgrey'))
        ax.add_patch(patches.Rectangle((20, 0), 100, 5, facecolor='white'))
        ax.set_xlim([-30 ,120])
        ax.set_ylim([0, 4])
        plt.show()

    def matplot_box_plots(self):
        """
        For doing background:
        https://stackoverflow.com/questions/18215276/how-to-fill-rainbow-color-under-a-curve-in-python-matplotlib
        """
        puncta_df = pd.read_csv(self.puncta, sep='\t', index_col=0)
        nopuncta_df = pd.read_csv(self.nopuncta, sep='\t', index_col=0)
        puncta_scores = list(puncta_df['LC Score'])
        nopuncta_scores = list(nopuncta_df['LC Score'])
        fig = plt.figure(figsize=(7.5, 3))
        ax = fig.add_subplot(111)
        #fig.set_facecolor('white')
        #ax.grid(False)
        ax.add_patch(patches.Rectangle((-30, 0), 30, 6, facecolor='grey'))
        ax.add_patch(patches.Rectangle((0, 0), 20, 6, facecolor='darkgrey'))
        ax.add_patch(patches.Rectangle((20, 0), 100, 6, facecolor='white'))
        ax.set_xlim([-30 ,110])
        ax.set_ylim([0, 4])
        labs = ['Does not Form Puncta', 'Forms Puncta']
        bp = {'color': 'black'}
        wp = {'color': 'black', 'linestyle':'-'}
        meanprops = dict(marker='o',
                              markeredgecolor='black',
                              markerfacecolor='black',
                              markersize=3)
        medianprops = dict(linestyle='-', color='black')
        all_scores = [nopuncta_scores, puncta_scores]
        ax.boxplot(all_scores,
                   vert=False,
                   whis=[5, 95],
                   labels=labs,
                   widths=0.5,
                   showmeans=True,
                   showfliers=False,
                   boxprops=bp,
                   whiskerprops=wp,
                   meanprops=meanprops,
                   medianprops=medianprops)
        #plt.xlim([-30, 120])
        plt.xticks(np.arange(-30, 111, 10))
        plt.xlabel('LC score')
        plt.tick_params(axis='both', left='on', top='on', right='on',
                        bottom='on', labelleft='off', labeltop='off',
                        labelright='on', labelbottom='on')
        plt.tight_layout()
        plt.show()

def main():
    ps = PlotScores()
    ps.matplot_box_plots()


if __name__ == '__main__':
    main()