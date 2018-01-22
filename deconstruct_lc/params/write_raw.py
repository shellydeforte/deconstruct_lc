import configparser
import os
import pandas as pd
from deconstruct_lc.params import lc_labels
from deconstruct_lc import tools_lc

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class WriteRaw(object):
    def __init__(self, k, seqs, init_dict):
        self.seqs = seqs
        self.init_dict = init_dict
        self.k = k

    def write_lca(self, alph):
        lc_labs = lc_labels.GetLabels(self.k)
        k_lcas = lc_labs.create_lcas(alph)
        df_dict = self.init_dict
        for k_lca in k_lcas:
            lca = k_lca.split('_')[1]
            scores = tools_lc.calc_lca_motifs(self.seqs, self.k, lca)
            df_dict[k_lca] = scores
        cols = ['Protein ID', 'Length', 'y']+k_lcas
        df = pd.DataFrame(df_dict, columns=cols)
        return df

    def write_lce(self):
        lc_labs = lc_labels.GetLabels(self.k)
        k_lces = lc_labs.create_lces(self.seqs)
        df_dict = self.init_dict
        for k_lce in k_lces:
            lce = float(k_lce.split('_')[1])
            scores = tools_lc.calc_lce_motifs(self.seqs, self.k, lce)
            df_dict[k_lce] = scores
        cols = ['Protein ID', 'Length', 'y']+k_lces
        df = pd.DataFrame(df_dict, columns=cols)
        return df


def main():
    pass


if __name__ == '__main__':
    main()