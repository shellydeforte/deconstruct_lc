import os
import pandas as pd
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc.analysis_bc.write_bc_score import BcScore
from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta


class FormatBiogrid(object):
    def __init__(self):
        config = read_config.read_config()
        data = config['fps']['data_dp']
        bio_dp = os.path.join(data, 'biogrid')
        self.bio_fp = os.path.join(bio_dp, 'BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.157.mitab.txt')
        self.pbody_fp = os.path.join(bio_dp, 'Pbody_annotations.txt')
        self.interactions_fp = os.path.join(bio_dp, 'pbody.tsv')
        self.yeast_scores = os.path.join(bio_dp, 'orf_pbody_scores.tsv')
        self.yeast_fasta = os.path.join(bio_dp, 'orf_trans.fasta')

    def get_pbody(self):
        df = pd.read_csv(self.pbody_fp, sep='\t')
        pids = list(df['Gene Systematic Name'])
        return set(pids)

    def read_biogrid(self):
        df_in = pd.read_csv(self.interactions_fp, sep='\t')
        print(len(set(df_in['A'])))
        print(len(set(df_in['B'])))

    def write_biogrid(self):
        pbodies = self.get_pbody()
        df_dict = {'A': [], 'B': []}
        df = pd.read_csv(self.bio_fp, sep='\t')
        for i, row in df.iterrows():
            yida = row['Alt IDs Interactor A']
            yidb = row['Alt IDs Interactor B']
            if (len(yida.split(':')) > 3) and (len(yidb.split(':')) > 3):
                pida = yida.split(':')[3].strip()
                pidb = yidb.split(':')[3].strip()
                if pida in pbodies and pidb in pbodies:
                    if pida != pidb:
                        df_dict['A'].append(pida)
                        df_dict['B'].append(pidb)
        df_out = pd.DataFrame(df_dict)
        df_out.drop_duplicates(inplace=True)
        df_out.to_csv(self.interactions_fp, sep='\t')
        print(df_dict)
        print(len(df_dict['A']))

    def get_scores(self):
        pbodies = self.get_pbody()
        pids, seqs = tools_fasta.fasta_to_id_seq(self.yeast_fasta)
        pseqs = []
        ppids = []
        for pid, seq in zip(pids, seqs):
            if pid in pbodies:
                pseqs.append(seq)
                ppids.append(pid)
        ns = NormScore()
        scores = ns.lc_norm_score(pseqs)
        df_dict = {'Protein ID': ppids, 'LC Score': scores}
        df_out = pd.DataFrame(df_dict)
        df_out.to_csv(self.yeast_scores, sep='\t')


def main():
    bg = FormatBiogrid()
    bg.write_biogrid()


if __name__ == '__main__':
    main()