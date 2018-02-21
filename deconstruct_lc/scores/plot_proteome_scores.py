import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import pandas as pd
from deconstruct_lc import read_config
import numpy as np


class PlotScores(object):
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.fpi = os.path.join(self.data_dp, 'scores', 'pdb_bc_scores.tsv')

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
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        #df = df[(df['Proteome'] == 'BC') | (df['Proteome'] == 'PDB')]
        bc_scores = list(df[df['Proteome'] == 'BC']['LC Score'])
        pdb_scores = list(df[df['Proteome'] == 'PDB']['LC Score'])
        #data = np.concatenate((bc_scores, pdb_scores), 0)
        fig = plt.figure(figsize=(7.5, 3))
        ax = fig.add_subplot(111)
        #fig.set_facecolor('white')
        #ax.grid(False)
        ax.add_patch(patches.Rectangle((-30, 0), 30, 6, facecolor='grey'))
        ax.add_patch(patches.Rectangle((0, 0), 20, 6, facecolor='darkgrey'))
        ax.add_patch(patches.Rectangle((20, 0), 100, 6, facecolor='white'))
        ax.set_xlim([-30 ,110])
        ax.set_ylim([0, 4])
        labs = ['PDB', 'Yeast', 'Yeast Nucleolus', 'Yeast Stress Granule', 'Yeast P Body']
        bp = {'color': 'black'}
        wp = {'color': 'black', 'linestyle':'-'}
        meanprops = dict(marker='o',
                              markeredgecolor='black',
                              markerfacecolor='black',
                              markersize=3)
        medianprops = dict(linestyle='-', color='black')
        all_scores = self.get_scores()
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

    def get_scores(self):
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        yeast = list(df[df['Proteome'] == 'Yeast']['LC Score'])
        yeast_sg = list(df[(df['Proteome'] == 'Cytoplasmic_Stress_Granule') &
                       (df['Organism'] == 'YEAST')]['LC Score'])
        yeast_pb = list(df[(df['Proteome'] == 'P_Body') &
                       (df['Organism'] == 'YEAST')]['LC Score'])
        yeast_nuc = list(df[(df['Proteome'] == 'Nucleolus') &
                        (df['Organism'] == 'YEAST')]['LC Score'])
        pdb = list(df[df['Proteome'] == 'PDB']['LC Score'])
        return [pdb, yeast, yeast_nuc, yeast_sg, yeast_pb]


def main():
    ps = PlotScores()
    ps.matplot_box_plots()


if __name__ == '__main__':
    main()