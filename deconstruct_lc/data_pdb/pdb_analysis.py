from Bio import SeqIO
import pandas as pd
import configparser
import re
import os
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc


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
        self.k_lca = 6
        self.k_lce = 6
        self.alph_lca = 'SGEQAPDTNKR'
        self.thresh_lce = 1.6
        self.lca_label = '{}_{}'.format(self.k_lca, self.alph_lca)
        self.lce_label = '{}_{}'.format(self.k_lce, self.thresh_lce)

    def write_analysis(self):
        df = pd.read_csv(self.pdb_all_df, sep='\t')
        df = df.drop_duplicates(subset=['Sequence', 'Missing'], keep=False)
        df = df.reset_index(drop=True)
        df = self.add_raw_score(df)
        df.to_csv(self.pdb_an_df, sep='\t')

    def add_raw_score(self, df):
        seqs = list(df['Sequence'])
        miss_seqs = list(df['Missing'])
        lcas = tools_lc.calc_lca_motifs(seqs, self.k_lca, self.alph_lca)
        lces = tools_lc.calc_lce_motifs(seqs, self.k_lce, self.thresh_lce)
        lcs = tools_lc.calc_lc_motifs(seqs, self.k_lca, self.alph_lca,
                                      self.thresh_lce)
        lengths = tools_fasta.get_lengths(seqs)
        miss = self.get_missing(miss_seqs)
        df['Length'] = lengths
        df['Miss Count'] = miss
        df[self.lca_label] = lcas
        df[self.lce_label] = lces
        df['LCA+LCE'] = lcs
        return df

    def get_missing(self, miss_seqs):
        miss = []
        for seq in miss_seqs:
            miss.append(seq.count('X'))
        return miss

    def get_pids(self):
        pids, seqs = tools_fasta.fasta_to_id_seq(self.pdb_miss_fp)
        return pids

    def write_full(self):
        """Write sequence, missing. Remove histags from both."""
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
                        seq, miss_seq = self.remove_histag(str(seq_rec.seq),
                                                           str(dis_rec.seq))
                        fo.write('{}\t{}\t{}\n'.format(pid, seq, miss_seq))

    def remove_histag(self, seq, miss_seq):
        """Remove histags from the sequence and the corresponding missing
        sequence"""
        regex = r'H{6}H*'
        match = re.finditer(regex, seq)
        indexes = []
        for item in match:
            indexes.append(item.start())
            indexes.append(item.end())
        nseq = self.slice_seq(indexes, seq)
        nmseq = self.slice_seq(indexes, miss_seq)
        return nseq, nmseq

    def slice_seq(self, indexes, seq):
        if len(indexes) > 0:
            nseq = seq[:indexes[0]]
            for i in range(1, len(indexes) - 1, 2):
                nseq += seq[indexes[i]:indexes[i + 1]]
            nseq += seq[indexes[-1]:]
        else:
            nseq = seq
        return nseq


def main():
    pa = PdbAnalysis()
    # pa.write_full() # This will take a while
    pa.write_analysis()


if __name__ == '__main__':
    main()