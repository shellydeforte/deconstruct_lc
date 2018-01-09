import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class PlotLcProteome():
    """
    Plot the fraction and standard deviation of each amino acid in low
    complexity regions. Data comes from lc_proteome_composition.py
    """

    def __init__(self):
        self.base_dir = os.path.join(os.path.dirname(__file__), '..',
                                     'data', 'lc_composition', 'data_out')
        self.fni = os.path.join(self.base_dir, 'lc_composition.tsv')

    def plot_lc_comp(self):
        aas, mean_lcs, stds = self._read_tsv()
        indexes = list(range(20))
        indexes = np.array(indexes)
        fig, ax = plt.subplots()
        ax.set_title('Fraction of each amino acid in LCRs', size=14)
        width = 0.5
        ax.set_xticks(indexes)
        ax.set_xticklabels(aas, size=14)
        ax.bar(indexes, mean_lcs, width, color='kkkkkkkkkkkkkkkkkkkk')
        ax.set_ylim([0, 0.25])
        for pos, y, err, color in zip(indexes, mean_lcs, stds,
                                      'kkkkkkkkkkkkkkkkkkkk'):
            ax.errorbar(pos, y, err, lw=1, capsize=3, capthick=1,
                        color=color)
        plt.show()

    def _read_tsv(self):
        df = pd.read_csv(self.fni, sep='\t', index_col=0)
        aas = list(df['Amino Acid'])
        mean_lcs = list(df['Mean'])
        stds = list(df['Standard Deviation'])
        return aas, mean_lcs, stds


def main():
    lcp = PlotLcProteome()
    lcp.plot_lc_comp()


if __name__ == '__main__':
    main()