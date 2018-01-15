from Bio import SeqIO
import configparser
import os
import pandas as pd
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class RawScores(object):
    def __init__(self, k_lca, k_lce, alph_lca, thresh_lce):
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.ss_dis_fp = os.path.join(self.pdb_dp, 'all_dis.fasta')
        self.fasta_fp = os.path.join(self.pdb_dp, 'pdb_norm.fasta')
        self.k_lca = k_lca
        self.k_lce = k_lce
        self.alph_lca = alph_lca
        self.thresh_lce = thresh_lce
        self.lca_label = '{}_{}'.format(self.k_lca, self.alph_lca)
        self.lce_label = '{}_{}'.format(self.k_lce, self.thresh_lce)

    def fasta_to_raw(self, is_pdb=True):
        pids, seqs = tools_fasta.fasta_to_id_seq(self.fasta_fp)
        assert len(pids) == len(seqs)
        lengths = self._get_lengths(seqs)
        lca_raw = tools_lc.calc_lca_motifs(seqs, self.k_lca, self.alph_lca)
        lce_raw = tools_lc.calc_lce_motifs(seqs, self.k_lce, self.thresh_lce)
        assert len(pids) == len(lca_raw) == len(lce_raw)
        cols = ['Protein ID', 'Length', self.lca_label, self.lce_label]
        df_dict = {'Protein ID': pids, 'Length': lengths,
                   self.lca_label: lca_raw, self.lce_label: lce_raw}
        df = pd.DataFrame(df_dict, columns=cols)
        if is_pdb:
            df = self._add_missing(df, pids)
        return df

    def _add_missing(self, df, pids):
        miss_dict = self._create_missing_dict(pids)
        miss_res = []
        for i, row in df.iterrows():
            pid = row['Protein ID']
            miss_res.append(miss_dict[pid].count('X'))
        df['Missing'] = miss_res
        return df

    def _create_missing_dict(self, pids):
        """Create df from pdb with pid, len, lca_raw, lce_raw"""
        ss_dis_dict = {}
        with open(self.ss_dis_fp, 'r') as ss_dis_in:
            for ss_dis_rec in SeqIO.parse(ss_dis_in, 'fasta'):
                pdb_chain = tools_fasta.id_cleanup(str(ss_dis_rec.id))
                dis_seq = str(ss_dis_rec.seq)
                if pdb_chain in pids:
                    ss_dis_dict[pdb_chain] = dis_seq
        return ss_dis_dict

    def _get_lengths(self, seqs):
        lengths = [len(seq) for seq in seqs]
        return lengths


def main():
    k_lca = 6
    k_lce = 6
    alph_lca = 'SGEQAPDTNKR'
    thresh_lce = 1.6
    rs = RawScores(k_lca, k_lce, alph_lca, thresh_lce)
    df = rs.fasta_to_raw()
    fno = 'pdbnorm_{}_{}.tsv'.format(rs.lca_label, rs.lce_label)
    fpo = os.path.join(config['filepaths']['data_fp'], 'scores', fno)
    df.to_csv(fpo, sep='\t')


if __name__ == '__main__':
    main()