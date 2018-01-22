import configparser
import os
import pandas as pd
from deconstruct_lc.params.write_raw import WriteRaw

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class PipeRaw(object):
    def __init__(self):
        self.train_fp = config['filepaths']['train_fp']
        self.param_dp = os.path.join(config['filepaths']['data_dp'], 'params')
        self.all_ids, self.all_seqs, self.all_lens, self.y = self.get_seqs()
        #self.kr = (3, 21)
        self.kr = (3, 4)
        #self.alph = 'SGEQAPDTNKRL'
        self.alph = 'SGE'

    def get_seqs(self):
        df = pd.read_csv(self.train_fp, sep='\t', index_col=0)
        all_ids = list(df['Protein ID'])
        all_seqs = list(df['Sequence'])
        all_lens = list(df['Length'])
        y = list(df['y'])
        return all_ids, all_seqs, all_lens, y

    def all_lca(self):
        df_dict = self.init_df_dict()
        for k in range(self.kr[0], self.kr[1]):
            fno = 'raw_{}_lca.tsv'.format(k)
            fpo = os.path.join(self.param_dp, fno)
            wr = WriteRaw(k, self.all_seqs, df_dict)
            df = wr.write_lca(self.alph)
            df.to_csv(fpo, sep='\t')

    def all_lce(self):
        df_dict = self.init_df_dict()
        for k in range(self.kr[0], self.kr[1]):
            fno = 'raw_{}_lce.tsv'.format(k)
            fpo = os.path.join(self.param_dp, fno)
            wr = WriteRaw(k, self.all_seqs, df_dict)
            df = wr.write_lce()
            df.to_csv(fpo, sep='\t')

    def init_df_dict(self):
        df_dict = {'Protein ID': self.all_ids, 'Length': self.all_lens,
                   'y': self.y}
        return df_dict

def main():
    pipe = PipeRaw()
    pipe.all_lca()
    pipe.all_lce()


if __name__ == '__main__':
    main()