import os
import pandas as pd

from deconstruct_lc.len_norm import len_norm


class RawMb(object):
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
        lca_m, lca_b = self.get_mb_lca(lca_dict['Label'])
        lca_dict['m'] = lca_m
        lca_dict['b'] = lca_b
        lce_m, lce_b = self.get_mb_lce(lce_dict['Label'])
        lce_dict['m'] = lce_m
        lce_dict['b'] = lce_b
        cols = ['Label', 'SVM score', 'm', 'b']
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

    def get_mb_lca(self, labels):
        ln = len_norm.LenNorm(self.config)
        ms = []
        bs = []
        for lab in labels:
            print(lab)
            lab_sp = lab.split('_')
            k = int(lab_sp[0])
            lca = lab_sp[1]
            m, b = ln.mb_lca(k, lca)
            ms.append(m)
            bs.append(b)
        return ms, bs

    def get_mb_lce(self, labels):
        ln = len_norm.LenNorm(self.config)
        ms = []
        bs = []
        for lab in labels:
            print(lab)
            lab_sp = lab.split('_')
            k = int(lab_sp[0])
            lce = float(lab_sp[1])
            m, b = ln.mb_lce(k, lce)
            ms.append(m)
            bs.append(b)
        return ms, bs


def main():
    pass


if __name__ == '__main__':
    main()