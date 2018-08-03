import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc.kelil.display_motif import Display

class MotifDisplay(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.bc_fp = os.path.join(data_dp, 'bc_analysis', 'bc_all_score.tsv')
        self.red_motif_fp = os.path.join(data_dp, 'kelil', 'rep_motifs_red.tsv')
        self.body_dp = os.path.join(data_dp, 'bc_analysis')
        self.kelil_dp = os.path.join(data_dp, 'kelil')
        self.allbc_out = os.path.join(data_dp, 'kelil', 'bc_all_motifs.tsv')

    def by_body(self):
        fns = ['Cajal_bodies_score.tsv', 'Centrosome_score.tsv',
               'Cytoplasmic_Stress_Granule_score.tsv', 'Nuclear_Speckles_score.tsv',
               'Nuclear_Stress_Granule_score.tsv', 'Nucleolus_score.tsv',
               'P_Body_score.tsv', 'Paraspeckle_score.tsv',
               'PML_Body_score.tsv']
        dm = Display(os.path.join(self.kelil_dp, 'Nucleolus_Serine.html'))
        for fn in ['Nucleolus_score.tsv']:
            print(fn)
            df = pd.read_csv(os.path.join(self.body_dp, fn), sep='\t')
            df = df[df['Organism'] == 'HUMAN']
            pids = list(df['Protein ID'])
            seqs = list(df['Sequence'])
            dm.write_body(pids, seqs)

    def by_score(self):
        df = pd.read_csv(self.bc_fp, sep='\t')
        df = df[df['Organism'] == 'HUMAN']
        low_df = df[df['LC Score'] < 0]
        hi_df = df[df['LC Score'] > 20]
        low_pids = list(low_df['Protein ID'])
        hi_pids = list(hi_df['Protein ID'])
        low_seqs = list(low_df['Sequence'])
        hi_seqs = list(hi_df['Sequence'])
        print(len(low_pids))
        print(len(hi_pids))
        dm = Display(os.path.join(self.kelil_dp, 'HighScore_Serine.html'))
        dm.write_body(hi_pids, hi_seqs)
        dm = Display(os.path.join(self.kelil_dp, 'LowScore_Serine.html'))
        dm.write_body(low_pids, low_seqs)

    def score_fun(self):
        pass


def main():
    md = MotifDisplay()
    md.by_score()


if __name__ == '__main__':
    main()