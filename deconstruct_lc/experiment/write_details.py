import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc.experiment.labels import Labels


class WriteDetails(object):
    """
    ORF, Gene, LC Score, Length, Description, Stress Granule, P body,
    'hydrolase', 'isomerase', 'ligase', 'lyase', 'oxidoreductase', 'transferase',
    RNA binding
    """
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.tht_fpi = os.path.join(data_dp, 'experiment', '180803_ThT.xls')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.desc_fpo = os.path.join(data_dp, 'experiment', 'full_descriptions.tsv')

    def read_files(self):
        df_dict = {'ORF': [], 'Gene': [], 'LC Score': [], 'Length': [],
                   'Description': [], 'cytoplasmic_stress_granule': [], 'Pbody': [],
                   'hydrolase': [], 'isomerase': [], 'ligase': [],
                   'lyase': [], 'oxidoreductase': [], 'transferase': [], 'RNA_binding': []}
        cols = ['Gene', 'ORF', 'LC Score', 'Length', 'Description', 'cytoplasmic_stress_granule',
                'Pbody', 'hydrolase', 'isomerase', 'ligase', 'lyase', 'oxidoreductase',
                'transferase', 'RNA_binding']
        l = Labels()
        enz_dict = l.get_enzyme_lists()
        lab_dict = l.get_labels()
        tht_df = pd.read_excel(self.tht_fpi, sheetname='Hoja3')
        for i, row in tht_df.iterrows():
            df_dict['ORF'].append(row['ORF'])
            df_dict['Gene'].append(row['plate 1'])
            length, desc, score = self.fetch_seq(row['ORF'])
            df_dict['Length'].append(length)
            df_dict['Description'].append(desc)
            df_dict['LC Score'].append(score)
            df_dict = self.fetch_labels(row['ORF'], enz_dict, lab_dict, df_dict)
        df_out = pd.DataFrame(df_dict, columns=cols)
        df_out.to_csv(self.desc_fpo, sep='\t')

    def fetch_seq(self, orf):
        result = tools_fasta.get_one_yeast_desc(self.orf_trans, orf)
        if result:
            ns = NormScore()
            seq = result[0]
            score = ns.lc_norm_score([seq])[0]
            length = len(seq)
            desc = result[1]
            return length, desc, score
        else:
            raise Exception("pid not present in fasta file")

    def fetch_labels(self, orf, enz_dict, lab_dict, df_dict):
        enzymes = ['hydrolase', 'isomerase', 'ligase', 'lyase',
                   'oxidoreductase', 'transferase']
        labels = ['RNA_binding', 'cytoplasmic_stress_granule', 'Pbody']
        for enz in enzymes:
            if orf in enz_dict[enz]:
                df_dict[enz].append('yes')
            else:
                df_dict[enz].append('no')
        for lab in labels:
            if orf in lab_dict[lab]:
                df_dict[lab].append('yes')
            else:
                df_dict[lab].append('no')
        return df_dict




def main():
    wd = WriteDetails()
    wd.read_files()


if __name__ == '__main__':
    main()