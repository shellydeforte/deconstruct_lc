import os
import pandas as pd

from deconstruct_lc import read_config


class Pilot(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.fpi = os.path.join(data_dp, 'experiment', '180614_glucose_starvation_2h.xlsx')
        self.fpo = os.path.join(data_dp, 'experiment', 'pilot.tsv')
        self.yeast_scores = os.path.join(data_dp, 'scores', 'all_yeast.tsv')

    def read_file(self):
        df = pd.read_excel(self.fpi, sheetname='Hoja1')
        print(df.head())
        yeast_df = pd.read_csv(self.yeast_scores, sep='\t', index_col=['ORF'])
        yeast_df = yeast_df[['Score']]
        yeast_dict = yeast_df.to_dict()
        yeast_dict['Score']['S288C'] = 'na'
        scores = []
        df_orfs = df['ORF']
        for orf in df_orfs:
            scores.append(yeast_dict['Score'][orf])
        df['Score'] = scores
        df.to_csv(self.fpo, sep='\t')


def main():
    p = Pilot()
    p.read_file()


if __name__ == '__main__':
    main()