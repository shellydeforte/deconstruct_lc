import os
import numpy as np
import pandas as pd
from scipy.stats import linregress

from deconstruct_lc import read_config
from deconstruct_lc.svm import svms
from deconstruct_lc import tools_lc


class NormSvm(object):
    def __init__(self, config):
        self.config = config
        data_dp = config['fps']['data_dp']
        pdb_dp = os.path.join(config['fps']['data_dp'], 'pdb_prep')
        pdb_fp = os.path.join(pdb_dp, 'pdb_norm_cd100.tsv')
        pdb_df = pd.read_csv(pdb_fp, sep='\t', index_col=0)
        self.pdb_seqs = list(pdb_df['Sequence'])
        self.pdb_lens = [len(pdb_seq) for pdb_seq in self.pdb_seqs]

        self.param_dp = os.path.join(data_dp, 'params')
        self.lce_lab_fp = os.path.join(self.param_dp, 'top_svm_lce.tsv')
        self.lca_lab_fp = os.path.join(self.param_dp, 'rep_lca.txt')

        train_fp = os.path.join(data_dp, 'train.tsv')
        train_df = pd.read_csv(train_fp, sep='\t', index_col=0)
        self.y = list(train_df['y'])
        self.seqs = list(train_df['Sequence'])

    def read_labels(self):
        lce_df = pd.read_csv(self.lce_lab_fp, sep='\t')
        lce_labs = list(lce_df['Label'])
        lca_labs = []
        with open(self.lca_lab_fp, 'r') as fpi:
            for line in fpi:
                lca_labs.append(line.strip())
        combos = []
        for lce in lce_labs:
            for lca in lca_labs:
                combos.append((lce, lca))
        return combos

    def run_svm(self):
        """
        if k is not the same
        """
        lce_lca = self.read_labels()
        for lce_lab, lca_lab in lce_lca:
            k_lce = int(lce_lab.split('_')[0])
            k_lca = int(lca_lab.split('_')[0])
            lce = float(lce_lab.split('_')[1])
            lca = str(lca_lab.split('_')[1])


def main():
    config = read_config.read_config()
    ns = NormSvm(config)
    ns.run_svm()


if __name__ == '__main__':
    main()