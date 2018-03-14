import os
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score


class BestFeatures(object):
    def __init__(self, config):
        data_dp = config['fps']['data_dp']
        self.solo_dp = os.path.join(data_dp, 'params', 'solo')
        self.combo_dp = os.path.join(data_dp, 'params', 'combos', 'norm')

    def ran_forest(self):
        X, y, all_labels = self.construct_matrix()
        for i in range(18):
            clf = RandomForestClassifier(n_estimators=100, random_state=15)
            clf = clf.fit(X, y)
            feat_imp = clf.feature_importances_
            indexes = np.argsort(feat_imp)
            X = self.remove_feature(indexes, X)
            kf = KFold(n_splits=3, random_state=0, shuffle=True)
            scores = cross_val_score(clf, X, y, cv=kf)
            print(scores.mean())
            for i in list(indexes):
                print(all_labels[i])
            print()

    def remove_feature(self, indexes, X):
        a = np.delete(X, [indexes[0]], 1)
        return a

    def construct_matrix(self):
        solo_fns = os.listdir(self.solo_dp)
        combo_fns = os.listdir(self.combo_dp)
        norm_scores = []
        all_labels = []
        for combo_fn in combo_fns:
            fpi = os.path.join(self.combo_dp, combo_fn)
            combo_df = pd.read_csv(fpi, sep='\t', index_col=0)
            labels = [label for label in list(combo_df)
                      if label != 'Protein ID' and label != 'y']
            for label in labels:
                full_label = '{} {}'.format(combo_fn[5:-4], label)
                all_labels.append(full_label)
                norm_score = list(combo_df[label])
                norm_scores.append(norm_score)
        for solo_fn in solo_fns:
            full_label = solo_fn[5:-4]
            all_labels.append(full_label)
            fpi = os.path.join(self.solo_dp, solo_fn)
            solo_df = pd.read_csv(fpi, sep='\t', index_col=0)
            norm_score = list(solo_df['Norm Scores'])
            norm_scores.append(norm_score)
        y = np.array(solo_df['y']).T
        X = np.array([norm_scores]).T.reshape(6793, 19)
        return X, y, all_labels