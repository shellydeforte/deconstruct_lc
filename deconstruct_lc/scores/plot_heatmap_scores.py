import matplotlib.pyplot as plt
import matplotlib
import configparser
import os
import pandas as pd
from deconstruct_lc.scores.norm_score import NormScore
import numpy as np

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class PlotHm(object):
    """"Try this with seaborn https://seaborn.pydata.org/generated/seaborn.heatmap.html"""
    def __init__(self):
        self.train_fpi = config['filepaths']['train_fp']

    def get_scores(self):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        df = df[df['y'] == 1][0:100]
        seqs = df['Sequence']
        ns = NormScore(seqs)
        nscores = ns.lc_norm_score()
        hist, bin_edges = np.histogram(nscores)
        hist = np.array([hist])
        nascores = np.array(nscores)
        fig, ax = plt.subplots()

        im = ax.imshow(hist, cmap=matplotlib.cm.RdBu, vmin=abs(
            hist).min(),
                       vmax=abs(hist).max(), extent=[0, 1, 0, 1])
        im.set_interpolation('bilinear')

        cb = fig.colorbar(im, ax=ax)
        plt.show()


def main():
    ph = PlotHm()
    ph.get_scores()


if __name__ == '__main__':
    main()