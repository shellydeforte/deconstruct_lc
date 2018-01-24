import configparser
import os
import pandas as pd
import numpy as np
from scipy.stats import linregress
from deconstruct_lc import tools_lc
from deconstruct_lc.svm import svms


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class NormRaw(object):
    """
    Get the raw scores specifically for the norm tsv file
    """
    def __init__(self):
        self.data_dp = os.path.join(config['filepaths']['data_dp'])
        self.pdb_dp = os.path.join(self.data_dp, 'pdb_prep')
        self.norm_fpi = os.path.join(self.pdb_dp, 'pdb_norm_cd100.tsv')
        self.train_fpi = config['filepaths']['train_fp']
        self.k = 6
        self.lca = 'SGEQAPDTNKR'
        self.lce = 1.6

    def read_m_b_r(self):
        norm_df = pd.read_csv(self.norm_fpi, sep='\t', index_col=0)
        lens = norm_df['Length']
        seqs = norm_df['Sequence']
        scores = self.calc_lc(seqs)
        print(linregress(lens, scores))
        # LinregressResult(slope=0.066213297264721263, intercept=1.7520712972708843, rvalue=0.7809848592680475, pvalue=0.0, stderr=0.0002720887990701889)
        scores = self.calc_lca(seqs)
        print(linregress(lens, scores))
        # LinregressResult(slope=0.049569081348185169, intercept=1.5446026962158523, rvalue=0.71825175189073842, pvalue=0.0, stderr=0.00024674540051569595)
        scores = self.calc_lce(seqs)
        print(linregress(lens, scores))
        # LinregressResult(slope=0.022592872067308627, intercept=0.37785221922942025, rvalue=0.66148080720662494, pvalue=0.0, stderr=0.00013162378880642137)

    def get_m_b(self):
        lc = (0.066213297264721263, 1.7520712972708843)
        lca = (0.049569081348185169, 1.5446026962158523)
        lce = (0.022592872067308627, 0.37785221922942025)
        return lc, lca, lce

    def norm_function(self, m, b, raw_scores, lengths):
        norm_scores = []
        for raw_score, length in zip(raw_scores, lengths):
            norm_score = raw_score - ((m * length) + b)
            norm_scores.append(norm_score)
        return norm_scores

    def calc_lca(self, seqs):
        scores = tools_lc.calc_lca_motifs(seqs, self.k, self.lca)
        return scores

    def calc_lce(self, seqs):
        scores = tools_lc.calc_lce_motifs(seqs, self.k, self.lce)
        return scores

    def calc_lc(self, seqs):
        scores = tools_lc.calc_lc_motifs(seqs, self.k, self.lca, self.lce)
        return scores

    def classify(self):
        """
        lc = 0.814
        lca = 0.810
        lce = 0.780
        lca + lce = 0.817
        lc + lca + lce = 0.820 (linear - 0.865 normal rbf)
        :return:
        """
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        seqs = df['Sequence']
        lengths = df['Length']
        scores = self.calc_lce(seqs)
        lc, lca, lce = self.get_m_b()

        m, b = lce
        lce_norm = self.norm_function(m, b, scores, lengths)
        scores = self.calc_lca(seqs)

        m, b = lca
        lca_norm = self.norm_function(m, b, scores, lengths)
        #X = np.array([lca_norm, lce_norm]).T.reshape(-1, 1)
        scores = self.calc_lc(seqs)

        m, b = lc
        lc_norm = self.norm_function(m, b, scores, lengths)
        X = np.array([lca_norm, lce_norm, lc_norm]).T
        y = np.array(df['y']).T
        clf = svms.normal_rbf(X, y)
        print(clf.score(X, y))

    def write_raw(self):
        pass


