import os
import pandas as pd

from deconstruct_lc.len_norm import mb_len_norm


class WriteMb(object):
    """
    Write the m, b normalization parameters for select lca/lce labels
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

    def write_mb_solo(self):
        lca_labs, lce_labs = self.read_labels()
        ln = mb_len_norm.LenNorm(self.config)
        df_dict = {'lc label': [], 'm': [], 'b': []}
        for lca_lab in lca_labs[0:1]:
            ll = lca_lab.split('_')
            k = int(ll[0])
            lca = str(ll[1])
            m, b = ln.mb_lca(k, lca)
            df_dict['lc label'].append(lca_lab)
            df_dict['m'].append(m)
            df_dict['b'].append(b)
        for lce_lab in lce_labs[0:1]:
            ll = lce_lab.split('_')
            k = int(ll[0])
            lce = float(ll[1])
            m, b = ln.mb_lce(k, lce)
            df_dict['lc label'].append(lce_lab)
            df_dict['m'].append(m)
            df_dict['b'].append(b)
        cols = ['lc label', 'm', 'b']
        df_out = pd.DataFrame(df_dict, columns=cols)
        df_out.to_csv(self.mb_solo_fp, sep='\t')

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