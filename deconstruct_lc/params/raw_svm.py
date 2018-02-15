import os
import numpy as np
import pandas as pd
from deconstruct_lc import read_config
from deconstruct_lc.svm import svms


class RawSvm(object):
    def __init__(self):
        self.config = read_config.read_config()
        self.param_dp = os.path.join(self.config['fps']['data_dp'], 'params')
        self.kr = (2, 21)
        #self.kr = (2, 3)

    def svm_lca_lce(self):
        for k in range(self.kr[0], self.kr[1]):
            print("{} LCE".format(k))
            lce_fpi = os.path.join(self.param_dp, 'raw_{}_lce.tsv'.format(k))
            lce_fpo = os.path.join(self.param_dp, 'svm_{}_lce.tsv'.format(k))
            self.raw_svm(lce_fpi, lce_fpo)
            print("{} LCA".format(k))
            lca_fpi = os.path.join(self.param_dp, 'raw_{}_lca.tsv'.format(k))
            lca_fpo = os.path.join(self.param_dp, 'svm_{}_lca.tsv'.format(k))
            self.raw_svm(lca_fpi, lca_fpo)

    def raw_svm(self, fpi, fpo):
        df_dict = {'SVM score': [], 'Label': []}
        cols = ['Label', 'SVM score']
        rem_cols = ['Protein ID', 'Length', 'y']
        df_in = pd.read_csv(fpi, sep='\t', index_col=0)
        k_lcs = [lab for lab in df_in.columns.values.tolist() if lab not
                    in rem_cols]
        for i, k_lc in enumerate(k_lcs):
            print(i)
            print(k_lc)
            raw_scores = df_in[k_lc]
            X = np.array([raw_scores]).T
            y = np.array(df_in['y']).T
            clf = svms.linear_svc(X, y)
            df_dict['SVM score'].append(clf.score(X, y))
            df_dict['Label'].append(k_lc)
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(fpo, sep='\t')


def main():
    pipe = RawSvm()
    pipe.svm_lca_lce()


if __name__ == '__main__':
    main()