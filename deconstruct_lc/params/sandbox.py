import os
import pandas as pd
import configparser
import matplotlib.pyplot as plt
from deconstruct_lc.params import write_raw
from deconstruct_lc.params import len_norm

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

data_dp = config['filepaths']['data_dp']
bc_fp = config['filepaths']['cb_fp']


def call_len_norm():
    pdb_df_fp = os.path.join(data_dp, 'scores',
                             'pdbmiss_6_SGEQAPDTNKR_6_1.6_raw.tsv')
    pdb_df = pd.read_csv(pdb_df_fp, sep='\t', index_col=0)
    k_lca = 6
    alph_lca = 'SGEQAPDTNKR'
    k_lce = 6
    thresh_lce = 1.6
    lca_label = '{}_{}'.format(k_lca, alph_lca)
    lce_label = '{}_{}'.format(k_lce, thresh_lce)
    ln = len_norm.PdbNorm(lca_label, lce_label, pdb_df)
    print(ln.mb_from_bins(True))
    #lens = list(pdb_df['Length'])
    #print(len(lens))
    #plt.hist(lens)
    #plt.show()


def main():
    call_len_norm()


if __name__ == '__main__':
    main()