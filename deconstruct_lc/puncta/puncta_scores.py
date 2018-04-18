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


class PunctaScores(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.puncta_fp = os.path.join(data_dp, 'puncta', 'marcotte', 'puncta_proteins.xlsx')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.agg_fp = os.path.join(data_dp, 'puncta', 'oconnel_agg_list')

    def agg_vs_puncta(self):
        ns = NormScore()
        labels = ['Foci (180)', 'Aggregates (117)', 'Foci-Aggregates (158)',
                  'Aggregates-Foci (95)', 'Aggregates&Foci (22)']
        agg_orfs = set(self.read_agg())
        puncta_orfs = set(self.get_ids('ST1'))
        agg = ns.lc_norm_score(self.get_seqs(agg_orfs))
        puncta = ns.lc_norm_score(self.get_seqs(puncta_orfs))
        puncta_agg = ns.lc_norm_score(self.get_seqs(puncta_orfs - agg_orfs))
        agg_puncta = ns.lc_norm_score(self.get_seqs(agg_orfs - puncta_orfs))
        puncta_and_agg = ns.lc_norm_score(self.get_seqs(puncta_orfs&agg_orfs))
        all_scores = [puncta, agg, puncta_agg, agg_puncta, puncta_and_agg]
        self.matplot_box_plots(all_scores, labels)

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

    def read_agg(self):
        orfs = []
        with open(self.agg_fp, 'r') as fi:
            for line in fi:
                orf = line.strip()
                orfs.append(orf)
        return orfs

    def run_plot(self):
        #st3_scores = self.insol_remove_puncta()
        st1_scores = self.get_scores('ST1')
        st2_scores = self.get_scores('ST2')
        st3_scores = self.get_scores('ST3')
        all_scores = [st1_scores, st2_scores, st3_scores]
        labs = ['ST1 (puncta)', 'ST2 (no puncta)', 'ST3 (insoluble->soluble)']
        self.matplot_box_plots(all_scores, labs)

    def cont_table_power(self):
        rows = 5
        cols = 2
        df = (rows - 1) * (cols - 1)
        nbins = df + 1
        alpha = 0.05
        power = 0.8
        st1_scores = self.get_scores('ST1')
        st2_scores = self.get_scores('ST2')
        col1 = self.bin_three(st1_scores)
        col2 = self.bin_three(st2_scores)
        n = sum(col1) + sum(col2)
        print(n)
        ct = np.array([col1, col2]).T
        print(ct)
        chi2, p, dof, ex = chi2_contingency(ct, correction=False)
        es = np.sqrt(chi2 / n * df)  # cramer's v
        print(es) # medium effect
        sample_size = smp.GofChisquarePower().solve_power(es, n_bins=nbins, alpha=alpha,
                                                power=power)
        print(sample_size)

    def bin_scores(self, scores):
        bins = [0, 0, 0]
        for score in scores:
            if score <= 0:
                bins[0] += 1
            elif 20 >= score > 0:
                bins[1] += 1
            else:
                bins[2] += 1
        return bins

    def bin_two(self, scores):
        bins = [0, 0]
        for score in scores:
            if score <= 0:
                bins[0] += 1
            else:
                bins[1] += 1
        return bins

    def bin_three(self, scores):
        bins = [0, 0, 0, 0, 0]
        for score in scores:
            if score <= -20:
                bins[0] += 1
            elif -20 < score <= -10:
                bins[1] += 1
            elif -10 < score <= 0:
                bins[2] += 1
            elif 0 < score <= 20:
                bins[3] += 1
            elif score > 20:
                bins[4] += 1
            else:
                print("Binning Problem")
        return bins

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

    def calc_scores(self, seqs):
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        return scores

    def get_seqs(self, orf_ids):
        seqs = tools_fasta.get_yeast_seq_from_ids(self.orf_trans, orf_ids)
        return seqs


def main():
    ps = PunctaScores()
    ps.agg_vs_puncta()


if __name__ == '__main__':
    main()