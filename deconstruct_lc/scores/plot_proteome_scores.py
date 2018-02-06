import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import configparser
import os
import pandas as pd
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
import numpy as np


class PlotScores(object):
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.fpi = os.path.join(self.data_dp, 'scores', 'proteomes.tsv')

    def matplot_box_plots(self):
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        #df = df[(df['Proteome'] == 'BC') | (df['Proteome'] == 'PDB')]
        bc_scores = list(df[df['Proteome'] == 'BC']['LC Score'])
        pdb_scores = list(df[df['Proteome'] == 'PDB']['LC Score'])
        #data = np.concatenate((bc_scores, pdb_scores), 0)
        fig, ax = plt.subplots(figsize=(7.5, 1))
        fig.set_facecolor('white')
        ax.grid(False)
        bp = {'color': 'black'}
        wp = {'color': 'black', 'linestyle':'-'}
        ax.boxplot([pdb_scores, bc_scores], vert=False, whis=[5, 95],
                    labels=['BC', 'PDB'], widths=0.5, boxprops=bp,
                   whiskerprops=wp)
        plt.xlim([-30, 150])
        plt.xticks(np.arange(-30, 151, 10))
        plt.tight_layout()
        plt.show()


class WriteScores(object):
    """
    Write the normed scores for yeast, human, pdb, bc for plotting
    """
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.train_fp = os.path.join(self.data_dp, 'train.tsv')
        self.yeast_fp = os.path.join(self.data_dp, 'proteomes',
                                     'UP000002311_559292_Yeast.fasta')
        self.human_fp = os.path.join(self.data_dp, 'proteomes',
                                     'UP000005640_9606_Human.fasta')
        self.fpo = os.path.join(self.data_dp, 'scores', 'proteomes.tsv')

    def write_scores(self):
        hum_seqs = tools_fasta.fasta_to_seq(self.human_fp, minlen=100,
                                            maxlen=2000)
        yeast_seqs = tools_fasta.fasta_to_seq(self.yeast_fp, minlen=100,
                                              maxlen=2000)
        df = pd.read_csv(self.train_fp, sep='\t', index_col=0)
        bc_seqs = list(df[df['y'] == 0]['Sequence'])
        pdb_seqs = list(df[df['y'] == 1]['Sequence'])
        ns = NormScore()
        hum_scores = ns.lc_norm_score(hum_seqs)
        yeast_scores = ns.lc_norm_score(yeast_seqs)
        bc_scores = ns.lc_norm_score(bc_seqs)
        pdb_scores = ns.lc_norm_score(pdb_seqs)
        lc_score = hum_scores + yeast_scores + bc_scores + pdb_scores
        proteome = ['Human']*len(hum_scores) + ['Yeast']*len(yeast_scores) + \
                   ['BC']*len(bc_scores) + ['PDB']*len(pdb_scores)
        df_dict = {'Proteome': proteome, 'LC Score': lc_score}
        cols = ['Proteome', 'LC Score']
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.fpo, sep='\t')


def main():
    #ws = WriteScores()
    #ws.write_scores()
    ps = PlotScores()
    ps.matplot_box_plots()


if __name__ == '__main__':
    main()