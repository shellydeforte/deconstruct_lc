import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore


class FormatMarcotte(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.marcotte_fpi = os.path.join(data_dp, 'experiment', 'marcotte_puncta_proteins.xlsx')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.fpo = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')


    def read_marcotte(self):
        df = pd.read_excel(self.marcotte_fpi, 'ST1')
        yeast_ids = list(df['ORF'])
        genes = list(df['Gene'])
        seqs = tools_fasta.get_yeast_seq_from_ids(self.orf_trans, yeast_ids)
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        df_out = pd.DataFrame({'Gene': genes, 'ORF': yeast_ids,
                               'LC Score': scores, 'Sequence': seqs},
                              columns=['Gene', 'ORF', 'LC Score', 'Sequence'])
        df_out.to_csv(self.fpo, sep='\t')


def main():
    fm = FormatMarcotte()
    fm.read_marcotte()


if __name__ == '__main__':
    main()