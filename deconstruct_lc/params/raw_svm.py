import numpy as np
import pandas as pd
from sklearn.svm import LinearSVC

def lc_class_scores(fi):
    df = pd.read_csv(fi, sep='\t')
    col_list = df.columns[4:]
    y = np.array(df['y']).T
    df_out = pd.DataFrame(0, index=col_list, columns=['SVM_score'])
    for col in col_list:
        X = np.array(df[col]).T.reshape(-1,1)
        score = run_linear_svc(X, y)
        df_out.loc[col,'SVM_score'] = score
    return df_out


def run_linear_svc(X, y):
    clf = LinearSVC(random_state=0)
    clf.fit(X, y)
    return clf.score(X, y)


def main():
    pass


if __name__ == '__main__':
    main()