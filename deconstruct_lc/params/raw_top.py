import os
import pandas as pd

from deconstruct_lc.len_norm import len_norm


class RawTop(object):
    def __init__(self, config):
        self.config = config
        data_dp = self.config['fps']['data_dp']
        self.param_dp = os.path.join(data_dp, 'params')
        self.k1 = self.config.getint('params', 'k1')
        self.k2 = self.config.getint('params', 'k2')

    def write_top(self):
        lca_fps, lce_fps = self.get_fps()
        lca_fpo = os.path.join(self.param_dp, 'top_svm_lca.tsv')
        lce_fpo = os.path.join(self.param_dp, 'top_svm_lce.tsv')
        lca_dict = self.get_top(lca_fps)
        lce_dict = self.get_top(lce_fps)
        cols = ['Label', 'SVM score']
        lca_df = pd.DataFrame(lca_dict, columns=cols)
        lca_df.to_csv(lca_fpo, sep='\t')
        lce_df = pd.DataFrame(lce_dict, columns=cols)
        lce_df.to_csv(lce_fpo, sep='\t')

    def get_top(self, all_fps):
        df_dict = {'Label': [], 'SVM score': []}
        for fp in all_fps:
            df = pd.read_csv(fp, sep='\t', index_col=0)
            ndf = df[df['SVM score'] > 0.82]
            df_dict['Label'] += list(ndf['Label'])
            df_dict['SVM score'] += list(ndf['SVM score'])
        return df_dict

    def get_fps(self):
        lca_fps = []
        lce_fps = []
        for k in range(self.k1, self.k2):
            lca_fpo = os.path.join(self.param_dp, 'svm_{}_lca.tsv'.format(k))
            lce_fpo = os.path.join(self.param_dp, 'svm_{}_lce.tsv'.format(k))
            lca_fps.append(lca_fpo)
            lce_fps.append(lce_fpo)
        return lca_fps, lce_fps

