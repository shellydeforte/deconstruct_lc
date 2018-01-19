import configparser
import os
import pandas as pd


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class ScoresNorm(object):
    """
    Explore fitted line for lca, lce, lc
    Test classification for each, normalized
    """
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.norm_fpo = os.path.join(self.pdb_dp, 'pdb_norm_cd100.tsv')
        self.train_fpo = os.path.join(self.pdb_dp, 'pdb_train_cd90.tsv')

    def regr(self):
        pass

def main():
    pass


if __name__ == '__main__':
    main()