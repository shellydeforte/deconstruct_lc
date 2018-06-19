import os
import pandas as pd
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta


class WriteNorm(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.yeast_scores = os.path.join(data_dp, 'scores', 'all_yeast.tsv')

    def write_yeast(self):
        pids, genes, seqs, descs = tools_fasta.get_pid_gene_desc_seq(self.orf_trans)
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        lengths = [len(seq) for seq in seqs]
        df_dict = {'ORF': pids, 'Gene': genes, 'Score': scores, 'Sequence': seqs, 'Length': lengths}
        df = pd.DataFrame(df_dict, columns=['ORF', 'Gene', 'Score', 'Sequence', 'Length'])
        df.to_csv(self.yeast_scores, sep='\t')


def main():
    wn = WriteNorm()
    wn.write_yeast()


if __name__ == '__main__':
    main()