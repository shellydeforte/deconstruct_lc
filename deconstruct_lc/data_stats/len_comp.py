import configparser
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from deconstruct_lc import tools_fasta
from deconstruct_lc.svm import svms

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class LenComp(object):
    """Run an SVM on the length and composition together and separately"""
    def __init__(self):
        self.aas = 'SGEQAPDTNKRLHVYFIMCW'
        self.lca = 'SGEQAPDTNKR'
        self.fd = os.path.join(config['filepaths']['data_fp'])
        self.cb90 = os.path.join(self.fd, 'bc_train_cd90.fasta')
        self.pdb90 = os.path.join(self.fd, 'pdb_train_cd90.fasta')
        self.pdb_ids, self.pdb_seqs = tools_fasta.fasta_to_id_seq(
            self.pdb90, minlen=100, maxlen=2000)
        self.bc_ids, self.bc_seqs = tools_fasta.fasta_to_id_seq(self.cb90,
                                                                minlen=100,
                                                                maxlen=2000)
        self.comp_len_fp = os.path.join(self.fd, 'len_comp', 'len_comp.tsv')

    def individual_aas(self):
        df = pd.read_csv(self.comp_len_fp, sep='\t', index_col=0)
        y = np.array(df['y']).T
        scores = []
        for aa in self.aas:
            cols = [aa]
            X = np.array(df[cols])
            lin_clf = svms.linear_svc(X, y)
            score = lin_clf.score(X, y)
            scores.append(score)
        print("The mean accuracy is {} and standard deviation is {} for the "
              "fraction of each amino acid used separately to "
              "classify.".format(np.mean(scores), np.std(scores)))
        all_cols = [aa for aa in self.aas]
        X = np.array(df[all_cols])
        rbf_clf = svms.normal_rbf(X, y)
        score = rbf_clf.score(X, y)
        print("The accuracy score for the fraction of all amino acids used "
              "to classify is {}".format(score))

    def lca_comp(self):
        df = pd.read_csv(self.comp_len_fp, sep='\t', index_col=0)
        y = np.array(df['y']).T
        lca_col = [aa for aa in self.lca]
        X = np.array(df[lca_col].sum(axis=1)).reshape(len(y), 1)
        lin_clf = svms.linear_svc(X, y)
        score = lin_clf.score(X, y)
        print("The accuracy using the LCA composition is {}".format(score))

    def lca_comp_plot(self):
        df = pd.read_csv(self.comp_len_fp, sep='\t', index_col=0)
        lca_col = [aa for aa in self.lca]
        dfz = df[df['y']==0]
        bc_lca_sum = list(dfz[lca_col].sum(axis=1))
        dfo = df[df['y']==1]
        pdb_lca_sum = list(dfo[lca_col].sum(axis=1))
        print("The BC mean and std are {} and {}".format(np.mean(
            bc_lca_sum), np.std(bc_lca_sum)))
        print("The PDB mean and std are {} and {}".format(np.mean(
            pdb_lca_sum), np.std(pdb_lca_sum)))
        labels = ['BC', 'PDB']
        plt.boxplot([bc_lca_sum, pdb_lca_sum], labels=labels)
        plt.ylabel("Fraction of LCA residues")
        plt.show()

    def write_aa_comp(self):
        cols = ['Protein ID', 'Length', 'y'] + [aa for aa in self.aas]
        df_dict = self.init_df(cols)
        for bc_seq in self.bc_seqs:
            a_bc_seq = ProteinAnalysis(bc_seq)
            bc_aas = a_bc_seq.get_amino_acids_percent()
            for aa in self.aas:
                df_dict[aa].append(bc_aas[aa])
        for pdb_seq in self.pdb_seqs:
            a_pdb_seq = ProteinAnalysis(pdb_seq)
            pdb_aas = a_pdb_seq.get_amino_acids_percent()
            for aa in self.aas:
                df_dict[aa].append(pdb_aas[aa])
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.comp_len_fp, sep='\t')

    def init_df(self, cols):
        df_dict = {}
        for col in cols:
            df_dict[col] = []
        pdb_lens = tools_fasta.get_lengths(self.pdb_seqs)
        bc_lens = tools_fasta.get_lengths(self.bc_seqs)
        pids = self.bc_ids + self.pdb_ids
        lens = bc_lens + pdb_lens
        y = [0]*len(self.bc_ids) + [1]*len(self.pdb_ids)
        df_dict['Protein ID'] = pids
        df_dict['Length'] = lens
        df_dict['y'] = y
        return df_dict


def main():
    lc = LenComp()
    lc.individual_aas()
    lc.lca_comp()
    lc.lca_comp_plot()



if __name__ == '__main__':
    main()