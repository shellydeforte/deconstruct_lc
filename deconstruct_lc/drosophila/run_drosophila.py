import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore


class Drosophila(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.fpi = os.path.join(data_dp, '..', 'drosophila_llps', 'candidate_drosophila.fasta')
        self.cand_fpo = os.path.join(data_dp, '..', 'drosophila_llps', 'candidate_drosophila.tsv')
        self.all_fpi = os.path.join(data_dp, '..', 'drosophila_llps',
                                    'all_drosophila.fasta')
        self.all_fpo = os.path.join(data_dp, '..', 'drosophila_llps', 'all_drosophila.tsv')

    def write_scores(self):
        ids, seqs = tools_fasta.fasta_to_id_seq(self.all_fpi)
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        df_out = pd.DataFrame({'Protein ID': ids,
                               'LC Score': scores},
                              columns=['Protein ID', 'LC Score'])
        df_out = df_out.sort_values(by='LC Score', ascending=False)
        print(df_out)
        df_out.to_csv(self.all_fpo, sep='\t')


def main():
    d = Drosophila()
    d.write_scores()


if __name__ == '__main__':
    main()