import os
import pandas as pd
import numpy as np

from deconstruct_lc import read_config
from deconstruct_lc.svm import svms
from deconstruct_lc.scores.norm_score import NormScore


class AdjustB(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.train_fpi = os.path.join(data_dp, 'train.tsv')

    def find_hyperplane(self):
        """Show that the hyperplane for the set intercept is close to 0"""
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        seqs = list(df['Sequence'])
        cs = NormScore()
        lc_norm = cs.lc_norm_score(seqs)
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