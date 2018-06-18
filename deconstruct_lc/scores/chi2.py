import os
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from itertools import combinations

from deconstruct_lc.analysis_bc.write_bc_score import BcScore
from deconstruct_lc import read_config


class BinChi2(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        scores_dp = os.path.join(data_dp, 'scores')
        self.prot_fpi = os.path.join(scores_dp, 'proteomes.tsv')
        self.fpi = os.path.join(scores_dp, 'pdb_bc_scores.tsv')
        self.yeast_bc = os.path.join(scores_dp, 'yeast_bc.tsv')
        self.human_bc = os.path.join(scores_dp, 'human_bc.tsv')
        self.prot = os.path.join(scores_dp, 'prot.tsv')
        self.yeast_out = os.path.join(scores_dp, 'yeast_chi2.tsv')
        self.human_out = os.path.join(scores_dp, 'human_chi2.tsv')

    def run_chi(self):
        self.load_files(['YEAST', 'Yeast'], self.yeast_out)
        self.load_files(['HUMAN', 'Human'], self.human_out)

    def load_files(self, org, fpo):
        chis = []
        df = self.create_bc(org[0])
        prot = self.create_prot(org[1])
        all = pd.concat([df, prot])
        prots = list(all['Proteome'])
        combs = list(combinations(prots, 2))
        for comb in combs:
            pval = self.get_chi2(all, comb)
            chis.append(pval)
        df_out = pd.DataFrame({'Comparison': combs, 'P-val': chis})
        df_out.to_csv(fpo, sep='\t')

    def get_chi2(self, df, comb):
        df = df[(df['Proteome'] == comb[0]) | (df['Proteome'] == comb[1])]
        df = df[['< 0', '0-20', '> 20']]
        df = np.array(df)
        return chi2_contingency(df)[1]

    def create_bc(self, org):
        names = []
        lts = []
        ms = []
        gts = []
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        bcs = BcScore()
        bc_names = bcs.get_sheets()
        for bc_name in bc_names:
            ndf = df[df['Proteome'] == bc_name]
            yndf = ndf[ndf['Organism'] == org]
            if len(yndf) > 0:
                lt, m, gt = self.get_bins(yndf)
                lts.append(lt)
                ms.append(m)
                gts.append(gt)
                names.append(bc_name)
        df_dict = {'Proteome': names, '< 0': lts, '0-20': ms, '> 20': gts}
        cols = ['Proteome', '< 0', '0-20', '> 20']
        df = pd.DataFrame(df_dict, columns=cols)
        return df

    def create_prot(self, name):
        df = pd.read_csv(self.prot_fpi, sep='\t', index_col=0)
        yndf = df[df['Proteome'] == name]
        lt, m, gt = self.get_bins(yndf)
        df_dict = {'Proteome': [name], '< 0': [lt], '0-20': [m], '> 20': [gt]}
        cols = ['Proteome', '< 0', '0-20', '> 20']
        df = pd.DataFrame(df_dict, columns=cols)
        return df

    def get_bins(self, df):
        ndf = df[df['LC Score'] < 0]
        lt = len(ndf)
        ndf = df[(df['LC Score'] >= 0) & (df['LC Score'] <= 20)]
        m = len(ndf)
        ndf = df[df['LC Score'] > 20]
        gt = len(ndf)
        return lt, m, gt


def main():
    bc = BinChi2()
    bc.run_chi()


if __name__ == '__main__':
    main()