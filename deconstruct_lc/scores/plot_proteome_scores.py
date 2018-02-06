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

class CreateTable(object):
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.fpi = os.path.join(self.data_dp, 'scores', 'proteomes.tsv')

    def get_stats(self):
        """Get numbers for < 0, 0-20, > 20"""
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        bc_df = df[df['Proteome'] == 'BC']
        pdb_df = df[df['Proteome'] == 'PDB']
        human_df = df[df['Proteome'] == 'Human']
        yeast_df = df[df['Proteome'] == 'Yeast']
        print("BC")
        self.get_bins(bc_df)
        print("PDB")
        self.get_bins(pdb_df)
        print("Human")
        self.get_bins(human_df)
        print("Yeast")
        self.get_bins(yeast_df)

    def get_bins(self, df):
        ndf = df[df['LC Score'] <= 0]
        lz = len(ndf)/len(df)
        print("The fraction < 0 is {}".format(lz))
        ndf = df[(df['LC Score'] > 0) & (df['LC Score'] <= 20)]
        gz = len(ndf)/len(df)
        print("The fraction > 0 and <= 20 is {}".format(gz))
        ndf = df[df['LC Score'] > 20]
        gt = len(ndf)/len(df)
        print("The fraction > 20 is {}".format(gt))


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
        fig, ax = plt.subplots(figsize=(7.5, 1.5))
        fig.set_facecolor('white')
        ax.grid(False)
        bp = {'color': 'black'}
        wp = {'color': 'black', 'linestyle':'-'}
        ax.boxplot([pdb_scores, bc_scores], vert=False, whis=[5, 95],
                    labels=['PDB', 'BC'], widths=0.5, boxprops=bp,
                   whiskerprops=wp)
        plt.xlim([-30, 150])
        plt.xticks(np.arange(-30, 151, 10))
        plt.xlabel('LC score')
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
    #ps = PlotScores()
    #ps.matplot_box_plots()
    ct = CreateTable()
    ct.get_stats()


if __name__ == '__main__':
    main()