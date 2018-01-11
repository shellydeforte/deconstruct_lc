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
        self.fd = os.path.join(config['filepaths']['data_fp'])
        self.cb90 = os.path.join(self.fd, 'bc_train_cd90.fasta')
        self.pdb90 = os.path.join(self.fd, 'pdb_train_cd90.fasta')
        self.pdb_ids, self.pdb_seqs = tools_fasta.fasta_to_id_seq(
            self.pdb90, minlen=100, maxlen=2000)
        self.bc_ids, self.bc_seqs = tools_fasta.fasta_to_id_seq(self.cb90,
                                                                minlen=100,
                                                                maxlen=2000)
        self.comp_len_fp = os.path.join(self.fd, 'len_comp', 'len_comp.tsv')

    def run_svm(self):
        X, y = self.construct_Xy()
        #clf = svms.smooth_rbf(X, y)
        #svms.plot_contour(clf, X, y)
        clf2 = svms.linear_svc(X, y)
        svms.plot_contour(clf2, X, y)

    def construct_Xy(self):
        df = pd.read_csv(self.comp_len_fp, sep='\t', index_col=0)
        y = np.array(df['y']).T
        aas = 'SGEQAPDTNKRLHVYFIMCW'
        #cols = ['Length'] + [aa for aa in aas]
        cols = ['Length']
        X = np.array(df[cols])
        return X, y

    def get_aa_comp(self):
        aas = 'SGEQAPDTNKRLHVYFIMCW'
        cols = ['Protein ID', 'Length', 'y'] + [aa for aa in aas]
        df_dict = self.init_df(cols)
        for pdb_seq in self.pdb_seqs:
            a_pdb_seq = ProteinAnalysis(pdb_seq)
            pdb_aas = a_pdb_seq.get_amino_acids_percent()
            for aa in aas:
                df_dict[aa].append(pdb_aas[aa])
        for bc_seq in self.bc_seqs:
            a_bc_seq = ProteinAnalysis(bc_seq)
            bc_aas = a_bc_seq.get_amino_acids_percent()
            for aa in aas:
                df_dict[aa].append(bc_aas[aa])
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
    lc.run_svm()



if __name__ == '__main__':
    main()