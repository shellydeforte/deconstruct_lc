import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from deconstruct_lc import read_config
from deconstruct_lc.svm import svms


class LenComp(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.aas = 'SGEQAPDTNKRLHVYFIMCW'
        self.lca = str(config['score']['lca'])
        self.train_fp = os.path.join(data_dp, 'train.tsv')
        self.comp_fp = os.path.join(data_dp, 'len_comp', 'train_comp.tsv')

    def plot_lencomp(self):
        plt.subplot(2, 1, 1)
        self.plot_len()
        plt.subplot(2, 1, 2)
        self.plot_comp()
        plt.subplots_adjust(hspace=0.5)
        plt.show()

    def plot_len(self):
        df_train = pd.read_csv(self.train_fp, sep='\t', index_col=0)
        bc_lens = list(df_train[df_train['y'] == 0]['Length'])
        pdb_lens = list(df_train[df_train['y'] == 1]['Length'])
        cb_heights, cb_bins = np.histogram(bc_lens, bins=20, range=(0,2000))
        cbn_heights = cb_heights / sum(cb_heights)

        pdb_heights, pdb_bins = np.histogram(pdb_lens, bins=20, range=(0,
                                                                       2000))
        pdbn_heights = pdb_heights / sum(pdb_heights)
        plt.bar(pdb_bins[:-1], pdbn_heights, width=(max(pdb_bins) - min(
            pdb_bins)) / len(pdb_bins), color="darkblue", alpha=0.7,
                label='PDB')

        plt.bar(cb_bins[:-1], cbn_heights, width=(max(cb_bins) - min(
            cb_bins)) / len(cb_bins), color="orangered", alpha=0.7,
                label='BC')

        plt.xlabel('Protein Length', size=12)
        plt.ylabel('Relative Fraction', size=12)
        plt.legend()

    def plot_comp(self):
        df_train = pd.read_csv(self.train_fp, sep='\t', index_col=0)
        bc_seqs = list(df_train[df_train['y'] == 0]['Sequence'])
        pdb_seqs = list(df_train[df_train['y'] == 1]['Sequence'])
        aas_list = [aa for aa in self.aas]
        ind = range(len(self.aas))
        pdb_seq = ''
        for seq in pdb_seqs:
            pdb_seq += seq
        cb_seq = ''
        for seq in bc_seqs:
            cb_seq += seq
        an_pdb_seq = ProteinAnalysis(pdb_seq)
        pdb_dict = an_pdb_seq.get_amino_acids_percent()
        an_cb_seq = ProteinAnalysis(cb_seq)
        cb_dict = an_cb_seq.get_amino_acids_percent()
        pdb_bins = []
        cb_bins = []
        for aa in aas_list:
            pdb_bins.append(pdb_dict[aa])
            cb_bins.append(cb_dict[aa])
        plt.bar(ind, pdb_bins, color='darkblue', alpha=0.7, label='PDB',
                align='center')
        plt.bar(ind, cb_bins, color='orangered', alpha=0.7,
                label='BC', align='center')
        plt.xticks(ind, aas_list)
        plt.xlim([-1, len(self.aas)])
        plt.legend()
        plt.xlabel('Amino Acids', size=12)
        plt.ylabel('Relative Fraction', size=12)

    def svm_comp(self):
        df_train = pd.read_csv(self.comp_fp, sep='\t', index_col=0)
        y = np.array(df_train['y']).T
        scores = []
        for aa in self.aas:
            cols = [aa]
            X = np.array(df_train[cols])
            lin_clf = svms.linear_svc(X, y)
            score = lin_clf.score(X, y)
            scores.append(score)
        print("The mean accuracy is {} and standard deviation is {} for the "
              "fraction of each amino acid used separately to "
              "classify.".format(np.mean(scores), np.std(scores)))
        all_cols = [aa for aa in self.aas]
        X = np.array(df_train[all_cols])
        rbf_clf = svms.smooth_rbf(X, y)
        score = rbf_clf.score(X, y)
        print("The accuracy score for the fraction of all amino acids used "
              "to classify is {}".format(score))

    def svm_len(self):
        df_train = pd.read_csv(self.train_fp, sep='\t', index_col=0)
        y = np.array(df_train['y']).T
        X = np.array(df_train[['Length']])
        lin_clf = svms.linear_svc(X, y)
        print("The accuracy score for the length is {}".format(lin_clf.score(X, y)))

    def write_aa_comp(self):
        cols = ['Protein ID', 'y'] + [aa for aa in self.aas]
        df_train = pd.read_csv(self.train_fp, sep='\t', index_col=0)
        bc_seqs = list(df_train[df_train['y'] == 0]['Sequence'])
        pdb_seqs = list(df_train[df_train['y'] == 1]['Sequence'])
        df_dict = dict()
        for aa in self.aas:
            df_dict[aa] = []
        df_dict['y'] = list(df_train['y'])
        df_dict['Protein ID'] = list(df_train['Protein ID'])
        for bc_seq in bc_seqs:
            a_bc_seq = ProteinAnalysis(bc_seq)
            bc_aas = a_bc_seq.get_amino_acids_percent()
            for aa in self.aas:
                df_dict[aa].append(bc_aas[aa])
        for pdb_seq in pdb_seqs:
            a_pdb_seq = ProteinAnalysis(pdb_seq)
            pdb_aas = a_pdb_seq.get_amino_acids_percent()
            for aa in self.aas:
                df_dict[aa].append(pdb_aas[aa])
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.comp_fp, sep='\t')


def main():
    lc = LenComp()
    lc.plot_lencomp()
    lc.write_aa_comp()
    lc.svm_comp()
    lc.svm_len()


if __name__ == '__main__':
    main()