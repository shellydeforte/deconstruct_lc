from Bio import SeqIO
import configparser
import os
import pandas as pd
from deconstruct_lc import tools_fasta


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class WriteDataTrain(object):
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.bc_dp = self.pdb_dp = os.path.join(config['filepaths'][
                                                    'data_dp'], 'bc_prep')
        self.pdb_fpi = os.path.join(self.pdb_dp, 'pdb_train_cd90.fasta')
        self.bc_fpi = os.path.join(self.bc_dp, 'bc_train_cd90.fasta')

    def train_df(self):
        pids, seqs = tools_fasta.fasta_to_id_seq(self.train_fpi)
        lens = tools_fasta.get_lengths(seqs)
        df_dict = {'Protein ID': pids, 'Sequence': seqs, 'Length': lens}
        cols = ['Protein ID', 'Sequence', 'Length']
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.train_fpo, sep='\t')


def main():
    pass


if __name__ == '__main__':
    main()