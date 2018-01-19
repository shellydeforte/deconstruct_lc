import configparser
import os
import pandas as pd
from scipy.stats import linregress
from deconstruct_lc import tools_fasta


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class WriteTrain(object):
    """
    Concatenate PDB and BC data into a single file
    PID, y, seq, len
    """
    def __init__(self):
        self.dp = os.path.join(config['filepaths']['data_dp'])
        self.pdb_fpi = os.path.join(self.dp, 'pdb_prep',
                                      'pdb_train_cd90.tsv')
        self.bc_fpi = os.path.join(self.dp, 'bc_prep', 'bc_train_cd90.fasta')
        self.fpo = os.path.join(self.dp, 'train.tsv')

    def concat_train(self):
        bc_pids, bc_seqs = tools_fasta.fasta_to_id_seq(self.bc_fpi)
        bc_lens = tools_fasta.get_lengths(bc_seqs)
        pdb_df = pd.read_csv(self.pdb_fpi, sep='\t', index_col=0)
        pdb_pids = list(pdb_df['Protein ID'])
        pdb_seqs = list(pdb_df['Sequence'])
        pdb_lens = list(pdb_df['Length'])
        pids = bc_pids + pdb_pids
        seqs = bc_seqs + pdb_seqs
        lens = bc_lens + pdb_lens
        y = [0]*len(bc_pids) + [1]*len(pdb_pids)
        cols = ['Protein ID', 'y', 'Length', 'Sequence']
        df_dict = {'Protein ID': pids, 'Sequence': seqs, 'Length': lens,
                   'y': y}
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.fpo, sep='\t')


def main():
    wt = WriteTrain()
    wt.concat_train()


if __name__ == '__main__':
    main()