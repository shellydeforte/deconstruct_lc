import pandas as pd
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc

class RawLca(object):
    def __init__(self, k, lca, fasta_fp):
        self.k = k
        self.lca = lca
        self.lca_lab = '{}_{}'.format(k, lca)
        self.fasta_fp = fasta_fp

    def get_raw(self):
        pids, seqs = tools_fasta.fasta_to_id_seq(self.fasta_fp)
        lengths = tools_fasta.get_lengths(seqs)
        scores = tools_lc.calc_lca_motifs(seqs, self.k, self.lca)
        return pids, scores, lengths

    def create_df(self):
        pids, scores, lengths = self.get_raw
        df_dict = {'Protein ID': pids, self.lca_lab: scores, 'Length': lengths}
        df_labs = ['Protein ID', self.lca_lab, 'Length']
        df = pd.DataFrame(df_dict, columns=df_labs)
        return df


class RawLce(object):
    def __init__(self, k, lce, fasta_fp):
        self.k = k
        self.lce = lce
        self.lce_lab = '{}_{}'.format(k, lce)
        self.fasta_fp = fasta_fp

    def get_raw(self):
        pids, seqs = tools_fasta.fasta_to_id_seq(self.fasta_fp)
        lengths = tools_fasta.get_lengths(seqs)
        scores = tools_lc.calc_lce_motifs(seqs, self.k, self.lce)
        return pids, scores, lengths

    def create_df(self):
        pids, scores, lengths = self.get_raw
        df_dict = {'Protein ID': pids, self.lce_lab: scores, 'Length': lengths}
        df_labs = ['Protein ID', self.lce_lab, 'Length']
        df = pd.DataFrame(df_dict, columns=df_labs)
        return df


def main():
    pass


if __name__ == '__main__':
    main()