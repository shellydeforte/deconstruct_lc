import os
import pandas as pd

from deconstruct_lc import tools_lc


class RawNorm(object):
    """
    Read the m/b values for each representative label and write the normalized
    score
    """
    def __init__(self, config):
        self.config = config
        data_dp = self.config['fps']['data_dp']
        self.param_dp = os.path.join(data_dp, 'params')
        self.combos_dp = os.path.join(self.param_dp, 'combos')
        self.solo_dp = os.path.join(self.param_dp, 'solo')
        self.mb_solo_fp = os.path.join(self.param_dp, 'mb_solo.tsv')
        train = os.path.join(data_dp, 'train.tsv')
        train_df = pd.read_csv(train, sep='\t', index_col=0)
        self.seqs = list(train_df['Sequence'])
        self.pids = list(train_df['Protein ID'])
        self.y = list(train_df['y'])
        self.lengths = list(train_df['Length'])

    def solo_norm(self):
        df_in = pd.read_csv(self.mb_solo_fp, sep='\t', index_col=0)
        df_in = df_in[df_in['pearsons'] > 0.7]
        for i, row in df_in.iterrows():
            fno = 'norm_{}.tsv'.format(str(row['lc label']))
            fpo = os.path.join(self.solo_dp, fno)
            m = float(row['m'])
            b = float(row['b'])
            lc_label = str(row['lc label'])
            print(lc_label)
            params = lc_label.split('_')
            k = int(params[0])
            if isinstance(params[1], str):
                lca = str(params[1])
                raw_scores = tools_lc.calc_lca_motifs(self.seqs, k, lca)
            else:
                lce = float(params[1])
                raw_scores = tools_lc.calc_lce_motifs(self.seqs, k, lce)
            norm_scores = self.norm_function(m, b, raw_scores, self.lengths)
            df_dict = {'Norm Scores': norm_scores, 'Protein ID': self.pids,
                       'y': self.y}
            df_out = pd.DataFrame(df_dict)
            df_out.to_csv(fpo, sep='\t')

    def combo_norm(self):
        fns = os.listdir(self.combos_dp)
        for fn in fns:
            fpi = os.path.join(self.combos_dp, fn)
            df_in = pd.read_csv(fpi, sep='\t', index_col=0)
            df_in = df_in[df_in['pearsons'] > 0.7]
            if len(df_in) > 0:
                df_dict = {'Protein ID': self.pids, 'y': self.y}
                fpo = os.path.join(self.combos_dp, 'norm_{}'.format(fn))
                params = fn.split('_')
                k = int(params[0])
                lce = float(params[1])
                lca = str(params[3])[:-4]
                for i, row in df_in.iterrows():
                    lc_label = str(row['LC Type'])
                    m = float(row['m'])
                    b = float(row['b'])
                    norm_scores = self.get_norm_scores(m, b, lc_label, k, lca, lce)
                    df_dict[lc_label] = norm_scores
                df_out = pd.DataFrame(df_dict)
                df_out.to_csv(fpo, sep='\t')

    def get_norm_scores(self, m, b, lc_label, k, lca, lce):
        raw_scores = self.get_raw_scores(lc_label, k, lca, lce)
        norm_scores = self.norm_function(m, b, raw_scores, self.lengths)
        return norm_scores

    def get_raw_scores(self, lc_label, k, lca, lce):
        if lc_label == 'LCA || LCE':
            scores = tools_lc.calc_lc_motifs(self.seqs, k, lca, lce)
        elif lc_label == 'LCA & LCE':
            scores = []
            for seq in self.seqs:
                scores.append(tools_lc.count_lca_and_lce(seq, k, lca, lce))
        elif lc_label == 'LCA & ~LCE':
            scores = []
            for seq in self.seqs:
                scores.append(tools_lc.count_lca_not_lce(seq, k, lca, lce))
        elif lc_label == '~LCA & LCE':
            scores = []
            for seq in self.seqs:
                scores.append(tools_lc.count_not_lca_lce(seq, k, lca, lce))
        else:
            raise Exception('Unexpected logical expression')
        return scores

    @staticmethod
    def norm_function(m, b, raw_scores, lengths):
        norm_scores = []
        for raw_score, length in zip(raw_scores, lengths):
            norm_score = raw_score - ((m * length) + b)
            norm_scores.append(norm_score)
        return norm_scores