import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore

class Proteins(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.hex_fpi = os.path.join(data_dp, 'experiment', '180803_HD.xls')
        self.tht_fpi = os.path.join(data_dp, 'experiment', '180803_ThT.xls')
        self.puncta_fpi = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')
        self.sg_ann = os.path.join(data_dp, 'experiment',
                                   'cytoplasmic_stress_granule_annotations.txt')
        self.sg_drop_out = os.path.join(data_dp, 'experiment', 'sg_drop_descriptions.tsv')
        self.sg_clus_out = os.path.join(data_dp, 'experiment', 'sg_clus_descriptions.tsv')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')

    def tht_stats(self):
        tht_df = pd.read_excel(self.tht_fpi, sheetname='Hoja3')
        print("There are {} proteins without data".format(
            len(tht_df[(tht_df['180708 48h'] == '-')])))
        print("There are {} proteins that did not form puncta".format(
            len(tht_df[(tht_df['180708 48h'] == 'no')])))
        print("There are {} proteins with no* (concentrated in nucleus)".format(
            len(tht_df[(tht_df['180708 48h'] == 'no*')])))

    def hex_stats(self):
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        print("There are {} proteins without data".format(
            len(hex_df[(hex_df['180708 48h'] == '-')])))
        print("There are {} proteins that did not form puncta".format(
            len(hex_df[(hex_df['180708 48h'] == 'no')])))

    def no_puncta(self):
        print("Proteins that did not form puncta at 48 hours")
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        hex_df = hex_df[(hex_df['180708 48h'] == 'no')]
        hex_ids = list(hex_df['ORF'])
        scores = self.fetch_scores(hex_ids)
        df = pd.DataFrame({'ORF': hex_ids, 'LC Score': scores})
        df = df.sort_values(by='LC Score', ascending=False)
        df = df.reset_index(drop=True)
        return df

    def yes_puncta_hex(self):
        print("Proteins that did not dissolve with hexandiol")
        tht_df = pd.read_excel(self.tht_fpi, sheetname='Hoja3')
        tht_df = tht_df[(tht_df['180803 48h HD 1h'] == 'yes') | (
        tht_df['180803 48h HD 1h'] == 'yes?') | (
                        tht_df['180803 48h HD 1h'] == 'no?')]
        hex_ids = list(tht_df['ORF'])
        scores = self.fetch_scores(hex_ids)
        df = pd.DataFrame({'ORF': hex_ids, 'LC Score': scores})
        df = df.sort_values(by='LC Score', ascending=False)
        df = df.reset_index(drop=True)
        print(len(df))
        return df

    def no_puncta_hex(self):
        print("Proteins that dissolved with hexandiol")
        tht_df = pd.read_excel(self.tht_fpi, sheetname='Hoja3')
        tht_df = tht_df[(tht_df['180708 48h'] == 'yes') | (tht_df['180708 48h'] == 'yes?')]
        tht_df = tht_df[(tht_df['180803 48h HD 1h'] == 'no') | (tht_df['180803 48h HD 1h'] == 'no*')]
        hex_ids = list(set(tht_df['ORF']))
        scores = self.fetch_scores(hex_ids)
        df = pd.DataFrame({'ORF': hex_ids, 'LC Score': scores})
        df = df.sort_values(by='LC Score', ascending=False)
        df = df.reset_index(drop=True)
        print(len(df))
        return df

    def yes_tht_stain(self):
        print("Proteins that stained with Tht")
        tht_df = pd.read_excel(self.tht_fpi, sheetname='Hoja3')
        tht_df = tht_df[(tht_df['180809 ThT'] == 'yes')]
        tht_ids = list(set(tht_df['ORF']))
        scores = self.fetch_scores(tht_ids)
        df = pd.DataFrame({'ORF': tht_ids, 'LC Score': scores})
        df = df.sort_values(by='LC Score', ascending=False)
        df = df.reset_index(drop=True)
        print(len(df))
        return df

    def no_tht_stain(self):
        print("Proteins that did not stain with Tht")
        tht_df = pd.read_excel(self.tht_fpi, sheetname='Hoja3')
        tht_df = tht_df[(tht_df['180809 ThT'] == 'no')]
        tht_ids = list(set(tht_df['ORF']))
        scores = self.fetch_scores(tht_ids)
        df = pd.DataFrame({'ORF': tht_ids, 'LC Score': scores})
        df = df.sort_values(by='LC Score', ascending=False)
        df = df.reset_index(drop=True)
        print(len(df))
        return df

    def clusters(self):
        print("Clusters are yes hexanediol and no Tht")
        tht_df = pd.read_excel(self.tht_fpi, sheetname='Hoja3')
        tht_df = tht_df[(tht_df['180809 ThT'] == 'no')]
        tht_df = tht_df[(tht_df['180803 48h HD 1h'] == 'yes') | (
        tht_df['180803 48h HD 1h'] == 'yes?') | (
                        tht_df['180803 48h HD 1h'] == 'no?')]
        tht_ids = list(set(tht_df['ORF']))
        scores = self.fetch_scores(tht_ids)
        df = pd.DataFrame({'ORF': tht_ids, 'LC Score': scores})
        df = df.sort_values(by='LC Score', ascending=False)
        df = df.reset_index(drop=True)
        print(len(df))
        return df

    def stress_granules(self):
        drop_df = self.no_puncta_hex()
        agg_df = self.yes_tht_stain()
        clus_df = self.clusters()
        sg = pd.read_csv(self.sg_ann, sep='\t')
        sg_orfs = list(sg['Gene Systematic Name'])
        sg_drop = drop_df[drop_df['ORF'].isin(sg_orfs)]
        sg_agg = agg_df[agg_df['ORF'].isin(sg_orfs)]
        sg_clus = clus_df[clus_df['ORF'].isin(sg_orfs)]
        sequences, genes, orfs, descriptions = tools_fasta.get_yeast_desc_from_ids(self.orf_trans, list(sg_clus['ORF']))
        lengths = [len(seq) for seq in sequences]
        ns = NormScore()
        scores = ns.lc_norm_score(sequences)
        ndf = pd.DataFrame({'Descriptions': descriptions, 'Gene': genes,
                            'Length': lengths, 'ORF': orfs, 'LC Score': scores})
        ndf.to_csv(self.sg_clus_out, sep='\t')
        sequences, genes, orfs, descriptions = tools_fasta.get_yeast_desc_from_ids(self.orf_trans, list(sg_drop['ORF']))
        lengths = [len(seq) for seq in sequences]
        ns = NormScore()
        scores = ns.lc_norm_score(sequences)
        ndf = pd.DataFrame({'Descriptions': descriptions, 'Gene': genes,
                            'Length': lengths, 'ORF': orfs, 'LC Score': scores})
        ndf.to_csv(self.sg_drop_out, sep='\t')
        return sg_orfs

    def fetch_scores(self, pids):
        df = pd.read_csv(self.puncta_fpi, sep='\t')
        df = df[df['ORF'].isin(pids)]
        df_orfs = set(list(df['ORF']))
        print(set(pids) - df_orfs)
        return list(df['LC Score'])

    def plotting(self):
        puncta_hex = self.yes_tht_stain()
        npuncta_hex = self.clusters()
        plt.hist(puncta_hex['LC Score'], bins=20, range=(-60, 200), alpha=0.5,
                 normed=True)
        plt.hist(npuncta_hex['LC Score'], bins=20, range=(-60, 200), alpha=0.5, normed=True)
        plt.show()


def main():
    p = Proteins()
    p.stress_granules()


if __name__ == '__main__':
    main()