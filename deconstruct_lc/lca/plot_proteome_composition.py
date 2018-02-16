import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from deconstruct_lc import read_config


class PlotLcProteome():
    """
    Plot the fraction and standard deviation of each amino acid in low
    complexity regions. Data comes from lc_proteome_composition.py
    """
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.fpi = os.path.join(data_dp, 'proteomes_analysis', 'lc_composition.tsv')

    def plot_lc_comp(self):
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        medprops, meanprops, whiskerprops, boxprops = self.params()
        df.plot.box(vert=True,
                    whis=[5, 95],
                    widths=0.75,
                    showfliers=True,
                    color='black',
                    patch_artist=False,
                    showmeans=True,
                    boxprops=boxprops,
                    whiskerprops=whiskerprops,
                    medianprops=medprops,
                    meanprops=meanprops
                    )
        #plt.title('Fraction of each amino acid in LCRs', size=14)
        plt.show()

    def params(self):
        medprops = dict(linestyle='-',
                        color='black')
        meanprops = dict(marker='o',
                         markeredgecolor='darkred',
                         markerfacecolor='darkred',
                         markersize=4)
        whiskerprops = dict(color='black',
                            linestyle='-')
        boxprops = dict(color='black')
        return medprops, meanprops, whiskerprops, boxprops


def main():
    lcp = PlotLcProteome()
    lcp.plot_lc_comp()


if __name__ == '__main__':
    main()