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
        self.train_fpo = os.path.join(config['filepaths']['data_dp'],
                                      'train.tsv')

    def train_df(self):
        pdb_pids, pdb_seqs = tools_fasta.fasta_to_id_seq(self.bc_fpi)
        pdb_lens = tools_fasta.get_lengths(pdb_seqs)
        bc_pids, bc_seqs = tools_fasta.fasta_to_id_seq(self.bc_fpi)
        bc_lens = tools_fasta.get_lengths(bc_seqs)
        lens = bc_lens + pdb_lens
        pids = bc_pids + pdb_pids
        seqs = bc_seqs + pdb_seqs
        y = [0]*len(bc_pids) + [1]*len(pdb_pids)
        df_dict = {'Protein ID': pids, 'Sequence': seqs, 'Length': lens,
                   'y': y}
        cols = ['Protein ID', 'y', 'Sequence', 'Length']
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.train_fpo, sep='\t')


def main():
    wt = WriteDataTrain()
    wt.train_df()


if __name__ == '__main__':
    main()