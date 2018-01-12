import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model


class PdbNorm(object):
    def __init__(self, lca_label, lce_label, pdb_df):
        self.lca_label = lca_label
        self.lce_label = lce_label
        self.pdb_df = pdb_df
        self.bin_range = [100, 700]

    def mb_from_bins(self, pl):
        lca_m, lca_b = self._score_bins(self.pdb_df, self.lca_label, pl)
        lce_m, lce_b = self._score_bins(self.pdb_df, self.lce_label, pl)
        cols = ['lca', 'lce', 'lca_m', 'lca_b', 'lce_m', 'lce_b']
        df_dict = {'lca': [self.lca_label], 'lce': [self.lce_label],
                   'lca_m': [lca_m], 'lca_b': [lca_b],
                   'lce_m': [lce_m], 'lce_b': [lce_b]}
        df = pd.DataFrame(df_dict)
        return df, cols

    def _score_bins(self, pdb_df, label, pl=False):
        means = []
        stds = []
        x = []
        for i in range(self.bin_range[0], self.bin_range[1], 50):
            df = pdb_df.loc[(pdb_df['Length'] <= i+50)
                             & (pdb_df['Length'] > i)]
            if len(df) > 0:
                dm = df[label].mean()
                x.append(i+25)
                means.append(dm)
                ds = df[label].std()
                stds.append(ds)
        regr = linear_model.LinearRegression()
        X = np.array(x).reshape(-1, 1)
        y = np.array(means).reshape(-1, 1)
        regr.fit(X, y)
        m = regr.coef_
        b = regr.intercept_
        if pl == True:
            self._plot_bins(m, b, x, means, stds)
        return m[0][0], b[0]

    def _plot_bins(self, m, b, x, means, stds):
        x_line = np.arange(self.bin_range[1])
        y_line = self._line_equation(x_line, m, b)
        plt.errorbar(x, means, stds, linestyle='None', marker='o')
        plt.plot(x_line, y_line[0])
        plt.xlabel('Protein Length')
        plt.ylabel('Mean score')
        plt.title('Mean score in bins of protein length = 50')
        plt.show()

    def _line_equation(self, x, m, b):
        y = m*x + b
        return y

def main():
    pass


if __name__ == '__main__':
    main()