import configparser
import os
import numpy as np
import pandas as pd
from deconstruct_lc.params import lc_labels
from deconstruct_lc.svm import svms

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class RawSvm(object):
    def __init__(self):
        self.param_dp = os.path.join(config['filepaths']['data_dp'], 'params')
        #self.kr = (3, 21)
        self.kr = (3, 4)
        #self.alph = 'SGEQAPDTNKRL'
        self.alph = 'SGE'

    def svm_lca(self):
        df_dict = {'SVM score': [], 'Label': []}
        cols = ['Label', 'SVM score']
        for k in range(self.kr[0], self.kr[1]):
            lc_labs = lc_labels.GetLabels(k)
            k_lcas = lc_labs.create_lcas(self.alph)
            fni = 'raw_{}_lca.tsv'.format(k, self.alph)
            fpi = os.path.join(self.param_dp, fni)
            df_in = pd.read_csv(fpi, sep='\t', index_col=0)
            for k_lca in k_lcas:
                scores = df_in[k_lca]
                X = np.array([scores]).T
                y = np.array(df_in['y']).T
                clf = svms.linear_svc(X, y)
                df_dict['SVM score'].append(clf.score(X, y))
                df_dict['Label'].append(k_lca)
            fno = 'svm_{}_lca.tsv'.format(k)
            fpo = os.path.join(self.param_dp, fno)
            df = pd.DataFrame(df_dict, columns=cols)
            df.to_csv(fpo, sep='\t')

    def svm_lce(self):
        df_dict = {'SVM score': [], 'Label': []}
        cols = ['Label', 'SVM score']
        rem_cols = ['Protein ID', 'Length', 'y']
        for k in range(self.kr[0], self.kr[1]):
            fni = 'raw_{}_lca.tsv'.format(k, self.alph)
            fpi = os.path.join(self.param_dp, fni)
            df_in = pd.read_csv(fpi, sep='\t', index_col=0)
            k_lces = [lab for lab in df_in.columns.values.tolist() if lab not
                    in rem_cols]
            for k_lce in k_lces:
                scores = df_in[k_lce]
                X = np.array([scores]).T
                y = np.array(df_in['y']).T
                clf = svms.linear_svc(X, y)
                df_dict['SVM score'].append(clf.score(X, y))
                df_dict['Label'].append(k_lce)
            fno = 'svm_{}_lce.tsv'.format(k)
            fpo = os.path.join(self.param_dp, fno)
            df = pd.DataFrame(df_dict, columns=cols)
            df.to_csv(fpo, sep='\t')


def main():
    pipe = RawSvm()
    pipe.svm_lca()
    pipe.svm_lce()


if __name__ == '__main__':
    main()