import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from deconstruct_lc import read_config

class Distance(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.fpi = os.path.join(data_dp, '..', 'lidice', 'distance locus-NE.xlsx')

    def read_file(self):
        df = pd.read_excel(self.fpi, sheetname='Hoja1')
        wt_active = df['WT active']
        wt_repressed = df['WT repressed'].dropna()
        grsIII = df['grsI,II'].dropna()
        self.plot_norm(wt_repressed, 'WT repressed')
        self.plot_norm(wt_active, 'WT active')
        self.plot_norm(grsIII, 'grsI,II$\Delta$')
        plt.xlabel('Distance locus-NE')
        plt.ylabel('P(D)')
        plt.legend()
        plt.show()

    def plot_norm(self, data, label):
        lnspc = np.linspace(-1.5, 1, len(data))
        m, s = stats.norm.fit(data)
        print(m)
        pdf_g = stats.norm.pdf(lnspc, m, s)
        plt.plot(lnspc, pdf_g, label=label, lw=2)


def main():
    d = Distance()
    d.read_file()


if __name__ == '__main__':
    main()