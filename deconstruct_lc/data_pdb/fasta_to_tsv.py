from Bio import SeqIO
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
        self.all_fpi = os.path.join(self.pdb_dp, 'pdb_all.fasta')
        self.norm_fpo = os.path.join(self.pdb_dp, 'pdb_norm_cd100.tsv')
        self.train_fpo = os.path.join(self.pdb_dp, 'pdb_train_cd90.tsv')
        self.all_fpo = os.path.join(self.pdb_dp, 'pdb_all.tsv')
        self.all_seq = os.path.join(self.pdb_dp, 'all_seqs.fasta')
        self.all_dis = os.path.join(self.pdb_dp, 'all_dis.fasta')
        self.all_ss = os.path.join(self.pdb_dp, 'all_ss.fasta')

    def write_tsv(self):
        self.train_df()
        self.write_full(self.all_fpi, self.all_fpo)
        self.write_full(self.norm_fpi, self.norm_fpo)

    def train_df(self):
        pids, seqs = tools_fasta.fasta_to_id_seq(self.train_fpi)
        lens = tools_fasta.get_lengths(seqs)
        df_dict = {'Protein ID': pids, 'Sequence': seqs, 'Length': lens}
        cols = ['Protein ID', 'Sequence', 'Length']
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.train_fpo, sep='\t')

    def df_stats(self):
        norm_df = pd.read_csv(self.norm_fpo, sep='\t', index_col=0)
        train_df = pd.read_csv(self.train_fpo, sep='\t', index_col=0)
        print("There are {} entries in norm_df".format(len(norm_df)))
        print("There are {} entries in train_df".format(len(train_df)))

    def write_full(self, fasta, fpo):
        """
        Write sequence, missing, secondary structure if in the list of pids.
        """
        all_pids = self.get_pids(fasta)
        count = 0
        with open(self.all_seq, 'r') as seq_fi, \
             open(self.all_dis, 'r') as dis_fi, \
             open(self.all_ss, 'r') as ss_fi:
            with open(fpo, 'w') as fo:
                fo.write('Protein ID\tSequence\tMissing\tSecondary '
                         'Structure\n')
                for seq_rec, dis_rec, ss_rec in zip(SeqIO.parse(seq_fi, 'fasta'),
                                                    SeqIO.parse(dis_fi, 'fasta'),
                                                    SeqIO.parse(ss_fi, 'fasta')):
                    pid = tools_fasta.id_cleanup(seq_rec.id)
                    if pid in all_pids:
                        count += 1
                        print(count)
                        seq = str(seq_rec.seq)
                        mseq = str(dis_rec.seq)
                        ss_seq = str(ss_rec.seq)
                        assert len(seq) == len(mseq) == len(ss_seq)
                        fo.write('{}\t{}\t{}\t{}\n'.format(pid, seq, mseq,
                                                           ss_seq))

    def get_pids(self, fasta):
        pids, seqs = tools_fasta.fasta_to_id_seq(fasta)
        return pids


def main():
    ft = FastaTsv()
    ft.write_tsv()


if __name__ == '__main__':
    main()