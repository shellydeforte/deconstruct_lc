import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore


class MarcotteScores(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.marcotte_fpi = os.path.join(data_dp, 'experiment',
                                         'marcotte_puncta_proteins.xlsx')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.puncta_fpo = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')
        self.npuncta_fpo = os.path.join(data_dp, 'experiment', 'marcotte_nopuncta_scores.tsv')

    def write_puncta(self):
        df = pd.read_excel(self.marcotte_fpi, 'ST1')
        yeast_ids = list(df['ORF'])
        genes = list(df['Gene'])
        seqs = tools_fasta.get_yeast_seq_from_ids(self.orf_trans, yeast_ids)
        lengths = [len(seq) for seq in seqs]
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        df_out = pd.DataFrame({'Gene': genes, 'ORF': yeast_ids,
                               'LC Score': scores, 'Sequence': seqs,
                               'Length': lengths},
                              columns=['Gene', 'ORF', 'LC Score', 'Length', 'Sequence'])
        print(df_out.head())
        df_out.to_csv(self.puncta_fpo, sep='\t')

    def write_nopuncta(self):
        """{'YEL014C', 'YDR250C', 'YOR199W', 'YJL017W'} are not included"""
        df = pd.read_excel(self.marcotte_fpi, 'NoPuncta')
        yeast_ids = list(df['ORF'])
        seqs, ngenes, orfs = tools_fasta.get_yeast_seq_gene_from_ids(self.orf_trans, yeast_ids)
        lengths = [len(seq) for seq in seqs]
        print(set(yeast_ids) - set(orfs))
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        df_out = pd.DataFrame({'Gene': ngenes, 'ORF': orfs,
                               'LC Score': scores, 'Sequence': seqs,
                               'Length': lengths},
                              columns=['Gene', 'ORF', 'LC Score', 'Length', 'Sequence'])
        print(df_out.head())
        df_out.to_csv(self.npuncta_fpo, sep='\t')


def main():
    ms = MarcotteScores()
    ms.write_nopuncta()



if __name__ == '__main__':
    main()