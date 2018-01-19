import configparser
import os
import pandas as pd
from deconstruct_lc import tools_fasta
from deconstruct_lc.data_pdb import histag


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class FastaTsv(object):
    """
    Fasta files are converted to tsv files with His-Tags removed
    """
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.norm_fpi = os.path.join(self.pdb_dp, 'pdb_norm_cd100.fasta')
        self.train_fpi = os.path.join(self.pdb_dp, 'pdb_train_cd90.fasta')
        self.norm_fpo = os.path.join(self.pdb_dp, 'pdb_norm_cd100.tsv')
        self.train_fpo = os.path.join(self.pdb_dp, 'pdb_train_cd90.tsv')

    def write_tsv(self):
        norm_df = self.new_df(self.norm_fpi)
        norm_df.to_csv(self.norm_fpo, sep='\t')
        train_df = self.new_df(self.train_fpi)
        train_df.to_csv(self.train_fpo, sep='\t')

    def new_df(self, fpi):
        pids, seqs = tools_fasta.fasta_to_id_seq(fpi)
        nseqs, num_his = histag.remove_histag(seqs)
        lens = tools_fasta.get_lengths(nseqs)
        df_dict = {'Protein ID': pids, 'Sequence': nseqs, 'Length': lens}
        cols = ['Protein ID', 'Sequence', 'Length']
        df = pd.DataFrame(df_dict, columns=cols)
        print("{} PDB chains had at least one His-Tag removed from {}".format(
            num_his, fpi))
        return df

    def df_stats(self):
        norm_df = pd.read_csv(self.norm_fpo, sep='\t', index_col=0)
        train_df = pd.read_csv(self.train_fpo, sep='\t', index_col=0)
        print("There are {} entries in norm_df".format(len(norm_df)))
        print("There are {} entries in train_df".format(len(train_df)))


def main():
    ft = FastaTsv()
    ft.write_tsv()
    ft.df_stats()


if __name__ == '__main__':
    main()