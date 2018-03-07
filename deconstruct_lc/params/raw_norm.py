import os
import pandas as pd

from deconstruct_lc.len_norm import mb_len_norm


class RawNorm(object):
    """
    Read the m/b values for each representative label and write the normalized
    score
    """
    def __init__(self, config):
        self.config = config
        data_dp = self.config['fps']['data_dp']
        self.param_dp = os.path.join(data_dp, 'params')
        self.lca_fpi = os.path.join(self.param_dp, 'rep_lca.txt')
        self.lce_fpi = os.path.join(self.param_dp, 'top_svm_lce.tsv')
        self.mb_solo_fp = os.path.join(self.param_dp, 'mb_solo.tsv')
        self.mb_combo_fp = os.path.join(self.param_dp, 'mb_combo.tsv')
        self.k1 = self.config.getint('params', 'k1')
        self.k2 = self.config.getint('params', 'k2')

    def norm_alone(self):
        """
        Calculate LCA or LCE only
        """
        lca_labs, lce_labs = self.read_labels()
        raw_scores, lengths, cy, cpids = self.get_raw_score(lca_labs[0], 'lca')
        for lca_lab in lca_labs:
            raw_scores, lengths, y, pids = self.get_raw_score(lca_lab, 'lca')
            assert y == cy
            assert pids == cpids
        for lce_lab in lce_labs:
            raw_scores, lengths, y, pids = self.get_raw_score(lce_lab, 'lce')
            assert y == cy
            assert pids == cpids

    def get_raw_score(self, label, type):
        """Given an lc label, get the corresponding raw scores"""
        lab_sp = label.split('_')
        k = lab_sp[0]
        fpi = os.path.join(self.param_dp, 'raw_{}_{}.tsv'.format(k, type))
        df_in = pd.read_csv(fpi, sep='\t', index_col=0)
        raw_scores = list(df_in[label])
        lengths = list(df_in['Length'])
        y = list(df_in['y'])
        pids = list(df_in['Protein ID'])
        return raw_scores, lengths, y, pids

    def read_labels(self):
        lce_df = pd.read_csv(self.lce_fpi, sep='\t')
        lce_labs = list(lce_df['Label'])
        lca_labs = []
        with open(self.lca_fpi, 'r') as fpi:
            for line in fpi:
                lca_labs.append(line.strip())
        return lca_labs, lce_labs

    def label_combos(self):
        lca_labs, lce_labs = self.read_labels()
        combos = []
        for lce in lce_labs:
            for lca in lca_labs:
                combos.append((lce, lca))
        return combos

    def read_files(self):
        lca_fpi = os.path.join(self.param_dp, 'top_svm_lca.tsv')
        lce_fpi = os.path.join(self.param_dp, 'top_svm_lce.tsv')
        lca_fpo = os.path.join(self.param_dp, 'top_norm_lca.tsv')
        lce_fpo = os.path.join(self.param_dp, 'top_norm_lce.tsv')
        lca_df = pd.read_csv(lca_fpi, sep='\t', index_col=0)
        lca_dict = self.calc_norm_score(lca_df, 'lca')
        lca_df_out = pd.DataFrame(lca_dict)
        lca_df_out.to_csv(lca_fpo, sep='\t')
        lce_df = pd.read_csv(lce_fpi, sep='\t', index_col=0)
        lce_dict = self.calc_norm_score(lce_df, 'lce')
        lce_df_out = pd.DataFrame(lce_dict)
        lce_df_out.to_csv(lce_fpo, sep='\t')

    def calc_norm_score(self, df_in, type):
        df_dict = {}
        for i, row in df_in.iterrows():
            m = float(row['m'])
            b = float(row['b'])
            label = str(row['Label'])
            print(label)
            raw_scores, lengths = self.get_raw_score(label, type)
            norm_scores = self.norm_function(m, b, raw_scores, lengths)
            df_dict[label] = norm_scores
        return df_dict

    def get_mb(self):
        pass

    @staticmethod
    def norm_function(m, b, raw_scores, lengths):
        norm_scores = []
        for raw_score, length in zip(raw_scores, lengths):
            norm_score = raw_score - ((m * length) + b)
            norm_scores.append(norm_score)
        return norm_scores