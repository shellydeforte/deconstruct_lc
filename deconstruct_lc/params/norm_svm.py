import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import cross_val_score

from deconstruct_lc.svm import svms


class NormSvm(object):
    """
    1. Read each normalized set separately
    2. Create a single vector and run both SVM and random forest
    3. Use random forest to identify best featuers
    """
    def __init__(self, config):
        config = config
        data_dp = config['fps']['data_dp']
        self.solo_dp = os.path.join(data_dp, 'params', 'solo')
        self.combo_dp = os.path.join(data_dp, 'params', 'combos', 'norm')
        self.oned_fpo = os.path.join(data_dp, 'params', 'oned_svm.tsv')

    def ran_forest(self):
        X, y = self.construct_matrix()
        print(X.shape)
        clf = RandomForestClassifier(n_estimators=10)
        clf = clf.fit(X, y)
        print(clf.score(X, y))
        feat_imp = clf.feature_importances_
        print(feat_imp)
        scores = cross_val_score(clf, X, y)
        print(scores.mean())

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
                print(full_label)
                all_labels.append(full_label)
                norm_score = list(combo_df[label])
                norm_scores.append(norm_score)
        for solo_fn in solo_fns:
            full_label = solo_fn[5:-4]
            print(full_label)
            all_labels.append(full_label)
            fpi = os.path.join(self.solo_dp, solo_fn)
            solo_df = pd.read_csv(fpi, sep='\t', index_col=0)
            norm_score = list(solo_df['Norm Scores'])
            norm_scores.append(norm_score)
        y = np.array(solo_df['y']).T
        X = np.array([norm_scores]).T.reshape(6793, 19)
        print(all_labels)
        print(all_labels[5])
        print(all_labels[13])
        print(all_labels[14])
        return X, y

    def oned_svm(self):
        solo_fns = os.listdir(self.solo_dp)
        combo_fns = os.listdir(self.combo_dp)
        df_dict = {'Label': [], 'Accuracy': []}
        for combo_fn in combo_fns:
            fpi = os.path.join(self.combo_dp, combo_fn)
            combo_df = pd.read_csv(fpi, sep='\t', index_col=0)
            labels = [label for label in list(combo_df)
                      if label != 'Protein ID' and label != 'y']
            for label in labels:
                full_label = '{} {}'.format(combo_fn[5:-4], label)
                df_dict['Label'].append(full_label)
                norm_scores = combo_df[label]
                X = np.array([norm_scores]).T
                y = np.array(combo_df['y']).T
                clf = svms.linear_svc(X, y)
                print(full_label)
                df_dict['Accuracy'].append(clf.score(X, y))
        for solo_fn in solo_fns:
            full_label = solo_fn[5:-4]
            print(full_label)
            df_dict['Label'].append(full_label)
            fpi = os.path.join(self.solo_dp, solo_fn)
            solo_df = pd.read_csv(fpi, sep='\t', index_col=0)
            norm_scores = solo_df['Norm Scores']
            X = np.array([norm_scores]).T
            y = np.array(solo_df['y']).T
            clf = svms.linear_svc(X, y)
            df_dict['Accuracy'].append(clf.score(X, y))
        df_out = pd.DataFrame(df_dict, columns=['Label', 'Accuracy'])
        df_out.to_csv(self.oned_fpo, sep='\t')


def main():
    pass


if __name__ == '__main__':
    main()