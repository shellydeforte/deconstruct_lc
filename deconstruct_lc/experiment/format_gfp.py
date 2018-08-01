import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore


class FormatGfp(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        experiment_dp = os.path.join(data_dp, 'experiment')
        self.gfp_fp = os.path.join(experiment_dp, 'allOrfData_yeastgfp.txt')
        self.marcotte_fpi = os.path.join(experiment_dp,
                                         'marcotte_puncta_proteins.xlsx')
        self.temp_out = os.path.join(experiment_dp, 'marcotte_notpuncta.tsv')

    def read_files(self):
        marc_df = pd.read_excel(self.marcotte_fpi, 'ST1')
        marc_ids = set(marc_df['ORF'])
        df = pd.read_csv(self.gfp_fp, sep='\t', index_col=False)
        df = df[df['localization summary'] == 'cytoplasm']
        df = df[~df['yORF'].isin(marc_ids)]
        ndf = pd.DataFrame({'Gene': df['gene name'], 'ORF': df['yORF']})
        ndf.to_csv(self.temp_out, sep='\t')


def main():
    fg = FormatGfp()
    fg.read_files()


if __name__ == '__main__':
    main()