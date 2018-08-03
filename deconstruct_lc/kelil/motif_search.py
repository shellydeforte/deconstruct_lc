import os
import pandas as pd

from deconstruct_lc import read_config

class MotifSearch(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.motif_fp = os.path.join(data_dp, 'kelil', 'REPEATING_MOTIFS.dat')
        self.bc_fp = os.path.join(data_dp, 'bc_analysis', 'bc_all_score.tsv')
        self.red_motif_fp = os.path.join(data_dp, 'kelil', 'rep_motifs_red.tsv')
        self.body_dp = os.path.join(data_dp, 'bc_analysis')
        self.kelil_dp = os.path.join(data_dp, 'kelil')
        self.allbc_out = os.path.join(data_dp, 'kelil', 'bc_all_motifs.tsv')

    def motifs_human(self):
        hum_pids = self.read_bc()
        df = pd.read_csv(self.red_motif_fp, sep='\t', index_col=0)
        df = df[df['PID'].isin(hum_pids)]
        df = df['MOT'].value_counts()
        df.to_csv(self.allbc_out, sep='\t')

    def read_motifs(self):
        df = self.get_motif_ids()
        df.to_csv(self.red_motif_fp, sep='\t')

    def get_motif_ids(self):
        df = pd.read_csv(self.motif_fp, sep='\t')
        reps = list(df['REP'])
        motifs = df['MOT']
        pids = []
        for i, row in df.iterrows():
            rpid = row['PRO']
            pid = rpid.split('|')[1]
            pids.append(pid)
            if i %100 == 0:
                print(i)
        ndf = pd.DataFrame({'PID': pids, 'REP': reps, 'MOT': motifs})
        return ndf

    def read_bc(self):
        df = pd.read_csv(self.bc_fp, sep='\t')
        df = df[df['Organism'] == 'HUMAN']
        return list(set(df['Protein ID']))

    def by_body(self):
        fns = ['Cajal_bodies_score.tsv', 'Centrosome_score.tsv',
               'Cytoplasmic_Stress_Granule_score.tsv', 'Nuclear_Speckles_score.tsv',
               'Nuclear_Stress_Granule_score.tsv', 'Nucleolus_score.tsv',
               'P_Body_score.tsv', 'Paraspeckle_score.tsv',
               'PML_Body_score.tsv']
        mot_df = pd.read_csv(self.red_motif_fp, sep='\t', index_col=0)
        for fn in fns:
            print(fn)
            df = pd.read_csv(os.path.join(self.body_dp, fn), sep='\t')
            df = df[df['Organism'] == 'HUMAN']
            hum_pids = df['Protein ID']
            nmot_df = mot_df[mot_df['PID'].isin(hum_pids)]
            ndf = nmot_df['MOT'].value_counts()
            fno = 'motifs_' + fn[:-9] + '.tsv'
            fpo = os.path.join(self.kelil_dp, fno)
            ndf.to_csv(fpo, sep='\t')

    def by_score(self):
        df = pd.read_csv(self.bc_fp, sep='\t')
        df = df[df['Organism'] == 'HUMAN']
        low_df = df[df['LC Score'] < 0]
        med_df = df[(df['LC Score'] >= 0) & (df['LC Score'] <= 20)]
        hi_df = df[df['LC Score'] > 20]
        low_pids = list(low_df['Protein ID'])
        med_pids = list(med_df['Protein ID'])
        hi_pids = list(hi_df['Protein ID'])
        print(len(low_pids))
        print(len(med_pids))
        print(len(hi_pids))

        # low_fp = os.path.join(self.kelil_dp, 'low_score_motifs.tsv')
        # med_fp = os.path.join(self.kelil_dp, 'med_score_motifs.tsv')
        # hi_fp = os.path.join(self.kelil_dp, 'high_score_motifs.tsv')
        #
        # mot_df = pd.read_csv(self.red_motif_fp, sep='\t', index_col=0)
        #
        # low_mot_df = mot_df[mot_df['PID'].isin(low_pids)]
        # low_mot_df = low_mot_df['MOT'].value_counts()
        # low_mot_df.to_csv(hi_fp, sep='\t')
        #
        # med_mot_df = mot_df[mot_df['PID'].isin(med_pids)]
        # med_mot_df = med_mot_df['MOT'].value_counts()
        # med_mot_df.to_csv(med_fp, sep='\t')
        #
        # hi_mot_df = mot_df[mot_df['PID'].isin(hi_pids)]
        # hi_mot_df = hi_mot_df['MOT'].value_counts()
        # hi_mot_df.to_csv(low_fp, sep='\t')


def main():
    ms = MotifSearch()
    ms.by_score()


if __name__ == '__main__':
    main()