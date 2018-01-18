# look at https://stackoverflow.com/questions/3949226/calculating-pearson-correlation-and-significance-in-python
import configparser
import os
from Bio import SeqIO
import random
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import numpy as np
import pandas as pd
from deconstruct_lc import motif_seq
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class MissMotif(object):
    """
    Plot 1: LC bins vs. fraction of protein with miss residue and avg
    missing residue when present  (plot_lc_vs_miss)
    Plot 2: LC bins vs. missing residues with fixed length
    Plot 3: LC bins vs. LC motif in missing residues (done)
    Do missing residues occur more frequently within blobs?
    """
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.pdb_an_dp = os.path.join(config['filepaths']['data_dp'],
                                      'pdb_analysis')
        self.pdb_an_fp = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')
        self.k_lca = 6
        self.k_lce = 6
        self.alph_lca = 'SGEQAPDTNKR'
        self.thresh_lce = 1.6
        self.lca_label = '{}_{}'.format(self.k_lca, self.alph_lca)
        self.lce_label = '{}_{}'.format(self.k_lce, self.thresh_lce)
        self.lc_vs_miss_fp = os.path.join(self.pdb_an_dp, 'lc_vs_miss.tsv')

    def motif_vs_coverage(self):
        df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        bin = range(100, 1000, 100)
        all_motif_percs = []
        all_motif_stds = []
        x = []
        for i in bin:
            motif_percs = []
            print(i)
            ndf = df[(df['Length'] >= i) & (df['Length'] < i + 100)]
            for i, row in ndf.iterrows():
                seq = row['Sequence']
                ind_in = self.get_inds(seq)
                motif_percs.append(len(ind_in) / len(seq))
            x.append(i)
            all_motif_percs.append(np.mean(motif_percs))
            all_motif_stds.append(np.std(motif_percs))
            print(all_motif_percs)
            print(all_motif_stds)

    def coverage_random(self):
        df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        bin = range(100, 1000, 100)
        all_motif_percs = []
        all_motif_stds = []
        x = []
        for i in bin:
            motif_percs = []
            print(i)
            ndf = df[(df['Length'] >= i) & (df['Length'] < i + 100)]
            aseq = list(ndf['Sequence'])
            seqs = self.create_random(aseq)
            for seq in seqs:
                ind_in = self.get_inds(seq)
                motif_percs.append(len(ind_in) / len(seq))
            x.append(i)
            all_motif_percs.append(np.mean(motif_percs))
            all_motif_stds.append(np.std(motif_percs))
            print(all_motif_percs)
            print(all_motif_stds)

    def create_random(self, seqs):
        nseqs = []
        for seq in seqs:
            lseq = [a for a in seq]
            random.shuffle(lseq)
            nseqs.append(''.join(lseq))
        return nseqs

    def plot_coverage(self):
        tmean = [0.23613392480397224, 0.24129299067670479, 0.20363240784003156,
         0.21521984605747985, 0.21560380306075025, 0.2126832223655015,
         0.20074931437836224, 0.19808298265652774, 0.20585607288722238]
        tstd = [0.098707938782962051, 0.093531289182195776,
               0.065869324533671433,
         0.059857030693938475, 0.055764428149622389, 0.052605567548994127,
         0.056823970944672043, 0.0423112359041719, 0.033929398744321125]
        x = [100, 200, 300, 400, 500, 600, 700, 800, 900]
        plt.xlim([0, 1000])
        plt.errorbar(x, tmean, tstd, linestyle='None', marker='o',
                     capsize=3)
        plt.show()

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

    def plot_fix_len(self):
        df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        df = df[(df['Length'] >= 400) & (df['Length'] < 600)]
        labels = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30',
                  '30-35', '35-40', '40-45', '45-50', '50+']
        bins = range(0, 50, 5)
        frac_w_miss = []
        num_miss = []
        std_num_miss = []
        for i in bins:
            print(i)
            ndf = df[(df['LCA+LCE'] >= i) & (df['LCA+LCE'] < i + 5)]
            nm_ndf = ndf[ndf['Miss Count'] > 0]
            frac_w_miss.append(len(nm_ndf)/len(ndf))
            num_miss.append(np.mean(list(nm_ndf['Miss Count'])))
            std_num_miss.append(np.std(list(nm_ndf['Miss Count'])))
        ndf = df[(df['LCA+LCE'] >= 50)]
        nm_ndf = ndf[ndf['Miss Count'] > 0]
        frac_w_miss.append(len(nm_ndf) / len(ndf))
        num_miss.append(np.mean(list(nm_ndf['Miss Count'])))
        std_num_miss.append(np.std(list(nm_ndf['Miss Count'])))
        x = list(range(len(frac_w_miss)))
        plt.xticks(x, labels, rotation=45)
        plt.xlim([-1, len(x)+1])
        plt.errorbar(x, num_miss, std_num_miss, linestyle='None', marker='o',
                     capsize=3, label='Average missing residues')
        plt.ylabel('Average Missing Residues')
        plt.xlabel('LC motifs')
        plt.ylim([0,200])
        plt.show()
        plt.ylabel('Fraction of proteins with missing residues')
        plt.scatter(x, frac_w_miss)
        plt.plot(x, frac_w_miss)
        plt.show()

    def write_lc_vs_miss(self):
        df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        labels = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30',
                  '30-35', '35-40', '40-45', '45-50', '50+']
        bins = range(0, 50, 5)
        frac_w_miss = []
        num_miss = []
        std_num_miss = []
        for i in bins:
            print(i)
            ndf = df[(df['LCA+LCE'] >= i) & (df['LCA+LCE'] < i + 5)]
            nm_ndf = ndf[ndf['Miss Count'] > 0]
            frac_w_miss.append(len(nm_ndf)/len(ndf))
            num_miss.append(np.mean(list(nm_ndf['Miss Count'])))
            std_num_miss.append(np.std(list(nm_ndf['Miss Count'])))
        ndf = df[(df['LCA+LCE'] >= 50)]
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

    def plot_box_whisker(self):
        """I'm not sure if a boxplot is better"""
        df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        bins = range(45, 50, 5)
        mm = []
        mp = []
        for i in bins:
            print(i)
            ndf = df[(df[self.lca_label] >= i) & (df[self.lca_label] < i+5)]
            print(len(ndf))
            miss_in_motifs, motif_percs = self.lc_blobs(ndf)
            mm.append(miss_in_motifs)
            mp.append(motif_percs)
        plt.boxplot([mm, mp])
        #plt.boxplot(mp)
        plt.ylim([-0.1, 1.1])
        plt.show()

    def read_df(self):
        df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        bins = range(0, 50, 5)
        mean_mm = []
        std_mm = []
        mean_mp = []
        std_mp = []
        for i in bins:
            print(i)
            ndf = df[(df['LCA+LCE'] >= i) & (df['LCA+LCE'] < i+5)]
            print(len(ndf))
            miss_in_motifs, motif_percs = self.lc_blobs(ndf)
            mean_mm.append(np.mean(miss_in_motifs))
            std_mm.append(np.std(miss_in_motifs))
            mean_mp.append(np.mean(motif_percs))
            std_mp.append(np.std(motif_percs))
        ndf = df[(df['LCA+LCE'] >= 50)]
        miss_in_motifs, motif_percs = self.lc_blobs(ndf)
        mean_mm.append(np.mean(miss_in_motifs))
        std_mm.append(np.std(miss_in_motifs))
        mean_mp.append(np.mean(motif_percs))
        std_mp.append(np.std(motif_percs))
        print(mean_mm)
        print(std_mm)
        print(mean_mp)
        print(std_mp)
        plt.errorbar(bins, mean_mm, std_mm, linestyle='None', marker='o')
        plt.errorbar(bins, mean_mp, std_mp, linestyle='None', marker='o')
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

    def lc_blobs(self, df):
        miss_in_motifs = []
        motif_percs = []
        for i, row in df.iterrows():
            miss = row['Missing']
            seq = row['Sequence']
            ind_miss = set([i for i, c in enumerate(miss) if c == 'X'])
            if len(ind_miss) > 0:
                ind_in = self.get_inds(seq)
                miss_in_motifs.append(len(ind_in & ind_miss) / len(ind_miss))
                motif_percs.append(len(ind_in)/len(seq))
        return miss_in_motifs, motif_percs

    def get_inds(self, seq):
        lcas = motif_seq.LcSeq(seq, self.k_lca, self.alph_lca, 'lca')
        lces = motif_seq.LcSeq(seq, self.k_lce, self.thresh_lce, 'lce')
        lca_in, lca_out = lcas._get_motif_indexes()
        lce_in, lce_out = lces._get_motif_indexes()
        ind_in = lca_in.union(lce_in)
        return ind_in


def main():
    mm = MissMotif()
    mm.coverage_random()


if __name__ == '__main__':
    main()