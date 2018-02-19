import os
import pandas as pd
from deconstruct_lc.params import lc_labels
from deconstruct_lc import tools_lc


class RawScores(object):
    def __init__(self, config):
        data_dp = config['fps']['data_dp']
        self.train_fp = os.path.join(data_dp, 'train.tsv')
        self.param_dp = os.path.join(data_dp, 'params')
        self.k1 = config.getint('params', 'k1')
        self.k2 = config.getint('params', 'k2')
        self.alph = config['params']['alph']
        self.all_ids, self.all_seqs, self.all_lens, self.y = self.get_seqs()

    def get_seqs(self):
        df = pd.read_csv(self.train_fp, sep='\t', index_col=0)
        all_ids = list(df['Protein ID'])
        all_seqs = list(df['Sequence'])
        all_lens = list(df['Length'])
        y = list(df['y'])
        return all_ids, all_seqs, all_lens, y

    def write_lca(self):
        for k in range(self.k1, self.k2):
            df_dict = {'Protein ID': self.all_ids, 'Length': self.all_lens,
                       'y': self.y}
            print("Now processing LCA raw scores for k = {}".format(k))
            fno = 'raw_{}_lca.tsv'.format(k)
            fpo = os.path.join(self.param_dp, fno)
            if not os.path.exists(fpo):
                df = self.create_df_lca(k, df_dict)
                df.to_csv(fpo, sep='\t')

    def create_df_lca(self, k, df_dict):
        lc_labs = lc_labels.GetLabels(k)
        k_lcas = lc_labs.create_lcas(self.alph)
        for k_lca in k_lcas:
            lca = k_lca.split('_')[1]
            scores = tools_lc.calc_lca_motifs(self.all_seqs, k, lca)
            df_dict[k_lca] = scores
        cols = ['Protein ID', 'Length', 'y']+k_lcas
        df = pd.DataFrame(df_dict, columns=cols)
        return df

    def write_lce(self):
        for k in range(self.k1, self.k2):
            df_dict = {'Protein ID': self.all_ids, 'Length': self.all_lens,
                       'y': self.y}
            print("Now processing LCE raw scores for k = {}".format(k))
            fno = 'raw_{}_lce.tsv'.format(k)
            fpo = os.path.join(self.param_dp, fno)
            if not os.path.exists(fpo):
                df = self.create_df_lce(k, df_dict)
                df.to_csv(fpo, sep='\t')

    def create_df_lce(self, k, df_dict):
        lc_labs = lc_labels.GetLabels(k)
        k_lces = lc_labs.create_lces(self.all_seqs)
        for k_lce in k_lces:
            lce = float(k_lce.split('_')[1])
            scores = tools_lc.calc_lce_motifs(self.all_seqs, k, lce)
            df_dict[k_lce] = scores
        cols = ['Protein ID', 'Length', 'y']+k_lces
        df = pd.DataFrame(df_dict, columns=cols)
        return df


def main():
    pass


if __name__ == '__main__':
    main()