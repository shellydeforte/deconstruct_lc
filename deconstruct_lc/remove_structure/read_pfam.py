import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from deconstruct_lc import read_config

class Pfam(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.pfam_fp = os.path.join(data_dp, 'pfam', 'Pfam-A.regions.uniprot.tsv')
        self.sgd_uni_fp = os.path.join(data_dp, 'proteomes', 'yeast_pd.xlsx')
        self.puncta = os.path.join(data_dp, 'experiment', 'marcotte_puncta_proteins.xlsx')
        self.puncta_map_excel = os.path.join(data_dp, 'experiment', 'puncta_map.xlsx')
        self.nopuncta_map = os.path.join(data_dp, 'experiment', 'nopuncta_map.tsv')
        self.nopuncta_map_excel = os.path.join(data_dp, 'experiment', 'nopuncta_map.xlsx')
        self.pfam_puncta = os.path.join(data_dp, 'experiment', 'puncta_pfam.tsv')
        self.pfam_nopuncta = os.path.join(data_dp, 'experiment', 'nopuncta_pfam.tsv')
        self.puncta_uni = os.path.join(data_dp, 'experiment', 'puncta_uni.txt')
        self.nopuncta_uni = os.path.join(data_dp, 'experiment', 'nopuncta_uni.txt')

    def read_file(self):
        """
        uniprot_acc
        pfamA_acc
        seq_start
        seq_end
        """
        unis = self.get_nopuncta_uni()
        sdf = pd.DataFrame()
        for i, chunk in enumerate(pd.read_csv(self.pfam_fp, sep='\t', chunksize=100000)):
            print(i)
            ndf = chunk[chunk['uniprot_acc'].isin(unis)]
            ndf = ndf[['uniprot_acc', 'pfamA_acc', 'seq_start', 'seq_end']]
            sdf = pd.concat([sdf, ndf])
        sdf.to_csv(self.pfam_nopuncta, sep='\t')

    def get_puncta_uni(self):
        df = pd.read_excel(self.puncta_map, sheetname='puncta_map')
        unis = list(df['Uni_ID'])
        return unis

    def get_nopuncta_uni(self):
        df = pd.read_excel(self.nopuncta_map_excel, sheetname='nopuncta_map')
        unis = list(df['Uni_ID'])
        return unis

    def write_nopuncta_map(self):
        puncta = pd.read_excel(self.puncta, sheetname='NoPuncta')
        no_puncta_orfs = list(puncta['ORF'])
        df = pd.read_excel(self.sgd_uni_fp, sep='\t')
        df = df[df['ORF'].isin(no_puncta_orfs)]
        df.to_csv(self.nopuncta_map, sep='\t')

    def write_uni(self):
        pdf = pd.read_excel(self.puncta_map_excel, sheetname='puncta_map')
        punis = set(list(pdf['Uni_ID']))
        with open(self.puncta_uni, 'w') as fo:
            for puni in punis:
                fo.write(puni+'\n')
        npdf = pd.read_excel(self.nopuncta_map_excel, sheetname='nopuncta_map')
        npunis = set(list(npdf['Uni_ID']))
        with open(self.nopuncta_uni, 'w') as fo:
            for npuni in npunis:
                fo.write(npuni+'\n')


def main():
    p = Pfam()
    p.write_uni()


if __name__ == '__main__':
    main()