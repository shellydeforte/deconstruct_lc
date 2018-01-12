"""
1. For proteins that are found in multiple BCs, is there anything special
about their scores?
Results. May be skewed towards higher scores, but they are spread out,
and there are still around 10 that sit in the PDB section. I would need to
look at this further to see what hypothesis I was trying to test. They do
NOT cluster though
2. Within a BC, is there anything that indicates complementarity?
Composition? Motifs? Composition within motifs?
Results: Scores are spread out - I didn't check by organism though,
there may be something present that a more detailed look could reveal
"""
import configparser
from collections import defaultdict
import os
import pandas as pd
import matplotlib.pyplot as plt

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class WithinBc(object):
    def __init__(self):
        self.minlen = 100
        self.maxlen = 2000
        self.fd = os.path.join(config['filepaths']['data_fp'])
        self.bc_fp = os.path.join(self.fd, 'scores',
                                  'quickgo_cb_cd90_6_SGEQAPDTNKR_6_1.6_norm.tsv')
        self.pdb_fp = os.path.join(self.fd, 'scores',
                                   'pdb_nomiss_cd90_6_SGEQAPDTNKR_6_1.6_norm')
        self.bc_ss = os.path.join(self.fd, 'bc_prep', 'quickgo_bc.xlsx')

    def get_within(self):
        fns = self.get_sheets()
        for sheet in fns:
            print(sheet)
            df_in = pd.read_excel(self.bc_ss, sheetname=sheet)
            pids = list(df_in['Protein ID'])
            lcas, lces = self.get_scores(pids)
            plt.scatter(lcas, lces)
            plt.xlim([-20, 100])
            plt.ylim([-20, 100])
            plt.show()

    def get_scores(self, pids):
        lcas = []
        lces = []
        df = pd.read_csv(self.bc_fp, sep='\t')
        for pid in pids:
            ndf = df[df['Protein ID'] == pid]
            if len(ndf) > 0:
                lcas.append(list(ndf['LCA Norm'])[0])
                lces.append(list(ndf['LCE Norm'])[0])
        return lcas, lces

    def find_overlap(self):
        """Find proteins that show up in more than one BC. Record the BCs,
        organism, protein"""
        pids = self.get_overlaps()
        df = pd.read_csv(self.bc_fp, sep='\t')
        lcas = []
        lces = []
        for pid in pids:
            ndf = df[df['Protein ID'] == pid]
            if len(ndf) > 0:
                lcas.append(list(ndf['LCA Norm'])[0])
                lces.append(list(ndf['LCE Norm'])[0])
        plt.scatter(lcas, lces)
        plt.xlim([-20, 100])
        plt.ylim([-20, 100])
        plt.show()

    def ss_to_dict(self):
        fns = self.get_sheets()
        pid_cb = defaultdict(set)
        for sheet in fns:
            df_in = pd.read_excel(self.bc_ss, sheetname=sheet)
            for row in df_in.itertuples():
                pid = row[1]
                pid_cb[pid].add(sheet)
        return pid_cb

    def get_overlaps(self):
        pid_cb = self.ss_to_dict()
        pids = []
        for pid in pid_cb:
            if len(pid_cb[pid]) >= 2:
                pids.append(pid)
        return pids

    def get_sheets(self):
        ex = pd.ExcelFile(self.bc_ss)
        sheet_names = ex.sheet_names
        return sorted(sheet_names)


def main():
    wb = WithinBc()
    wb.get_within()


if __name__ == '__main__':
    main()