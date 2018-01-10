import configparser
import os
import matplotlib.pyplot as plt
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

    def get_seqs(self):
        pdb_seqs = tools_fasta.fasta_to_seq(self.pdb90)
        pdb_lens = tools_fasta.get_lengths(pdb_seqs)
        bc_seqs = tools_fasta.fasta_to_seq(self.cb90)
        bc_lens = tools_fasta.get_lengths(bc_seqs)
        return pdb_seqs, pdb_lens, bc_seqs, bc_lens

    def construct_Xy(self):
        lca_nmo, lce_nmo, lca_pdb, lce_pdb = self.norm_scores()
        nmo_arr = np.array([lca_nmo, lce_nmo]).T
        pdb_arr = np.array([lca_pdb, lce_pdb]).T
        X = np.concatenate((nmo_arr, pdb_arr))
        y = np.array([0]*len(lca_nmo) + [1]*len(lca_pdb)).T
        return X, y

    def get_aa_comp(self):
        aas = 'SGEQAPDTNKRLHVYFIMCW'
        df_dict = {'Protein ID', 'Length'}
        cols = ['Protein ID', 'Length'] + [aa for aa in aas]
        print(cols)
        #pdb_seqs, pdb_lens, bc_seqs, bc_lens = self.get_seqs()
        #an_pdb_seq = ProteinAnalysis(pdb_seq)
        #pdb_dict = an_pdb_seq.get_amino_acids_percent()
        #an_cb_seq = ProteinAnalysis(cb_seq)
        #cb_dict = an_cb_seq.get_amino_acids_percent()
        #pdb_bins = []
        #cb_bins = []
        #for aa in aas:
        #    pdb_bins.append(pdb_dict[aa])
        #    cb_bins.append(cb_dict[aa])


def main():
    lc = LenComp()
    lc.get_aa_comp()



if __name__ == '__main__':
    main()