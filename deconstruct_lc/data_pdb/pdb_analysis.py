from Bio import SeqIO
import pandas as pd
import configparser
import os
from deconstruct_lc import tools_fasta


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class PdbAnalysis(object):
    def __init__(self):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.pdb_miss_fp = os.path.join(self.pdb_dp, 'pdb_norm_all.fasta')
        self.all_seq = os.path.join(self.pdb_dp, 'all_seqs.fasta')
        self.all_dis = os.path.join(self.pdb_dp, 'all_dis.fasta')
        self.pdb_all_df = os.path.join(self.pdb_dp, 'pdb_norm_all.tsv')
        self.pdb_an_df = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')

    def write_analysis(self):
        df = pd.read_csv(self.pdb_all_df, sep='\t')
        df = df.drop_duplicates(subset=['Sequence', 'Missing'], keep=False)
        df = df.reset_index()
        df.to_csv(self.pdb_an_df, sep='\t')

    def get_pids(self):
        pids, seqs = tools_fasta.fasta_to_id_seq(self.pdb_miss_fp)
        return pids

    def write_full(self):
        all_pids = self.get_pids()
        count = 0
        with open(self.all_seq, 'r') as seq_fi, \
             open(self.all_dis, 'r') as dis_fi:
            with open(self.pdb_all_df, 'w') as fo:
                fo.write('Protein ID\tSequence\tMissing\n')
                for seq_rec, dis_rec in zip(SeqIO.parse(seq_fi, 'fasta'),
                                        SeqIO.parse(dis_fi, 'fasta')):
                    pid = tools_fasta.id_cleanup(seq_rec.id)
                    if pid in all_pids:
                        count += 1
                        print(count)
                        miss_seq = str(dis_rec.seq)
                        seq = str(seq_rec.seq)
                        fo.write('{}\t{}\t{}\n'.format(pid, seq, miss_seq))


def main():
    pa = PdbAnalysis()
    # pa.write_full() # This will take a while
    pa.write_analysis()


if __name__ == '__main__':
    main()