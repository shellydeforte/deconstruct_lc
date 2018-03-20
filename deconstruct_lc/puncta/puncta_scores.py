import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc.display import display_lc


class PunctaScores(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.puncta_fp = os.path.join(data_dp, 'puncta', 'puncta_proteins.xlsx')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')

    def run_plot(self):
        #st3_scores = self.insol_remove_puncta()
        st1_scores = self.get_scores('ST1')
        st2_scores = self.get_scores('ST2')
        st3_scores = self.get_scores('ST3')
        all_scores = [st1_scores, st2_scores, st3_scores]
        labs = ['ST1 (puncta)', 'ST2 (no puncta)', 'ST3 (insoluble->soluble)']
        self.matplot_box_plots(all_scores, labs)

    def run_display(self):
        st1_ids = self.get_ids('ST1')
        st2_ids = self.get_ids('ST2')
        st3_ids = self.get_ids('ST3')
        st1_seqs, st1_genes = tools_fasta.get_yeast_seq_gene_from_ids(self.orf_trans, st1_ids)
        st2_seqs, st2_genes = tools_fasta.get_yeast_seq_gene_from_ids(self.orf_trans, st2_ids)
        st3_seqs, st2_genes = tools_fasta.get_yeast_seq_gene_from_ids(self.orf_trans, st3_ids)
        disp = display_lc.Display(st1_seqs, 'st1.html')
        disp.write_body()
        disp = display_lc.Display(st2_seqs, 'st2.html')
        disp.write_body()
        disp = display_lc.Display(st3_seqs, 'st3.html')
        disp.write_body()
        disp = display_lc.Display(st1_seqs, 'st1_color.html', color=True)
        disp.write_body()
        disp = display_lc.Display(st2_seqs, 'st2_color.html', color=True)
        disp.write_body()
        disp = display_lc.Display(st3_seqs, 'st3_color.html', color=True)
        disp.write_body()


    def insol_remove_puncta(self):
        st3_ids = self.get_ids('ST3')
        print(len(st3_ids))
        st1_ids = self.get_ids('ST1')
        no_puncta = set(st3_ids) - set(st1_ids)
        print(len(no_puncta))
        seqs = self.get_seqs(no_puncta)
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        return scores

    def matplot_box_plots(self, scores, labs):
        """
        For doing background:
        https://stackoverflow.com/questions/18215276/how-to-fill-rainbow-color-under-a-curve-in-python-matplotlib
        """
        fig = plt.figure(figsize=(7.5, 3))
        ax = fig.add_subplot(111)
        ax.add_patch(patches.Rectangle((-30, 0), 30, 6, facecolor='grey'))
        ax.add_patch(patches.Rectangle((0, 0), 20, 6, facecolor='darkgrey'))
        ax.add_patch(patches.Rectangle((20, 0), 100, 6, facecolor='white'))
        ax.set_xlim([-30 ,110])
        ax.set_ylim([0, 4])
        bp = {'color': 'black'}
        wp = {'color': 'black', 'linestyle':'-'}
        meanprops = dict(marker='o',
                              markeredgecolor='black',
                              markerfacecolor='black',
                              markersize=3)
        medianprops = dict(linestyle='-', color='black')
        ax.boxplot(scores,
                   vert=False,
                   whis=[5, 95],
                   widths=0.5,
                   labels=labs,
                   showmeans=True,
                   showfliers=False,
                   boxprops=bp,
                   whiskerprops=wp,
                   meanprops=meanprops,
                   medianprops=medianprops)
        plt.xticks(np.arange(-30, 111, 10))
        plt.xlabel('LC score')
        plt.tick_params(axis='both', left='on', top='on', right='on',
                        bottom='on', labelleft='off', labeltop='off',
                        labelright='on', labelbottom='on')
        plt.tight_layout()
        plt.show()

    def get_scores(self, sn):
        orf_ids = self.get_ids(sn)
        seqs = self.get_seqs(orf_ids)
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        return scores

    def get_ids(self, sn):
        df = pd.read_excel(self.puncta_fp, sheetname=sn)
        orf_ids = list(df['ORF'])
        return orf_ids

    def get_seqs(self, orf_ids):
        seqs = tools_fasta.get_yeast_seq_from_ids(self.orf_trans, orf_ids)
        return seqs


def main():
    ps = PunctaScores()
    ps.run_display()


if __name__ == '__main__':
    main()