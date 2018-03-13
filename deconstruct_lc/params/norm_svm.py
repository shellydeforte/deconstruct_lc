import os
import numpy as np
import pandas as pd

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