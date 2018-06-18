import os
import pandas as pd
import matplotlib.pyplot as plt

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore


class PunctaExperiment(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        print(data_dp)
        self.gfp_xl = os.path.join(data_dp, 'experiment', 'GFP collection.xls')
        self.gfp_orfs = os.path.join(data_dp, 'experiment', 'GFP_ORFs.tsv')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.exp_out = os.path.join(data_dp, 'experiment', 'selected_proteins.tsv')

    def load_xl(self):
        df_in = pd.read_csv(self.gfp_orfs)
        all_orfs = list(df_in['ORF Name'])
        seqs, genes, orfs = tools_fasta.get_yeast_seq_gene_from_ids(self.orf_trans, all_orfs)
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        lengths = [len(seq) for seq in seqs]
        df = pd.DataFrame({'Sequence': seqs, 'Gene': genes, 'ORF': orfs, 'Score': scores, 'Length': lengths},
                          columns=['ORF', 'Gene', 'Score', 'Length', 'Sequence'])
        df_len = df[(df['Length'] >= 500) & (df['Length'] <= 1000)]
        self.sample_bins(df_len)

    def sample_bins(self, df):
        ndf = pd.DataFrame(columns=['ORF', 'Gene', 'Score', 'Length', 'Sequence'])
        for i in range(-50, 20, 3):
            dfl = df[(df['Score'] >= i) & (df['Score'] < i+3)]
            if len(dfl) > 0:
                ndf = pd.concat([ndf, dfl.sample(n=1)])
        for i in range(20, 300, 7):
            dfl = df[(df['Score'] >= i) & (df['Score'] < i+7)]
            if len(dfl) > 0:
                ndf = pd.concat([ndf, dfl.sample(n=1)])
        ndf = ndf.sort_values(by='Score', ascending=True)
        ndf.reset_index(drop=True, inplace=True)
        ndf.to_csv(self.exp_out, sep='\t')


def main():
    pe = PunctaExperiment()
    pe.load_xl()


if __name__ == '__main__':
    main()