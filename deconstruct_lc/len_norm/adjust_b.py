import configparser
import os
import pandas as pd
import numpy as np
from deconstruct_lc.svm import svms
from deconstruct_lc.scores.norm_score import NormScore


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class AdjustB(object):
    def __init__(self):
        self.train_fpi = config['filepaths']['train_fp']

    def find_hyperplane(self):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        seqs = list(df['Sequence'])
        cs = NormScore(seqs)
        lc_norm = cs.lc_norm_score()
        X = np.array([lc_norm]).T
        y = np.array(df['y']).T
        clf = svms.linear_svc(X, y)
        xs = np.arange(-2, 2, 0.01).reshape(1, -1).T
        dists = list(clf.decision_function(xs))
        for x, dist in zip(xs, dists):
            if dist < 0:
                print(x)
                break


def main():
    ab = AdjustB()
    ab.find_hyperplane()


if __name__ == '__main__':
    main()