import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc import tools_lc


class Hexandiol(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.marcotte_fpi = os.path.join(data_dp, 'experiment',
                                         'marcotte_puncta_proteins.xlsx')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.hex_fpi = os.path.join(data_dp, 'experiment', '180803_HD.xls')
        self.puncta_fpi = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')
        self.npuncta_fpi = os.path.join(data_dp, 'experiment', 'marcotte_nopuncta_scores.tsv')
        self.nofasta = os.path.join(data_dp, 'experiment', 'hex_nop.fasta')
        self.yesfasta = os.path.join(data_dp, 'experiment', 'hex_yesp.fasta')

    def score_hist(self):
        lc_df = pd.read_csv(self.puncta_fpi, sep='\t', index_col=0)
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        hex_df = hex_df[(hex_df['180708 48h'] == 'yes') & (hex_df['180706 6h '] == 'no')]
        no_df = hex_df[(hex_df['180803 48h HD 1h'] == 'no')]
        no_orf = list(no_df['ORF'])
        yes_df = hex_df[(hex_df['180803 48h HD 1h'] == 'yes')]
        yes_orf = list(yes_df['ORF'])
        no_lc = lc_df[lc_df['ORF'].isin(no_orf)]
        yes_lc = lc_df[lc_df['ORF'].isin(yes_orf)]
        yes_scores = list(yes_lc['LC Score'])
        no_scores = list(no_lc['LC Score'])
        plt.xlabel('LC scores')
        plt.ylabel('Number of proteins')
        plt.hist(yes_scores, bins=20, range=(-60, 200), alpha=0.5, label='Does not dissolve with hexanediol')
        plt.hist(no_scores, bins=20, range=(-60, 200), alpha=0.5, label='Dissolves with hexanediol')
        plt.legend()
        plt.show()

    def write_fasta(self):
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        hex_df = hex_df[(hex_df['180708 48h'] == 'yes') & (hex_df['180706 6h '] == 'no')]
        no_df = hex_df[(hex_df['180803 48h HD 1h'] == 'no')]
        yes_df = hex_df[(hex_df['180803 48h HD 1h'] == 'yes')]
        no_orf = list(no_df['ORF'])
        yes_orf = list(yes_df['ORF'])
        tools_fasta.yeast_write_fasta_from_ids(self.orf_trans, no_orf, self.nofasta)
        tools_fasta.yeast_write_fasta_from_ids(self.orf_trans, yes_orf, self.yesfasta)

    def read_files(self):
        lc_df = pd.read_csv(self.puncta_fpi, sep='\t', index_col=0)
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        no_df = hex_df[(hex_df['180803 48h HD 1h'] == 'no') & (hex_df['180708 48h'] == 'yes') & (hex_df['180706 6h '] == 'no')]
        yes_df = hex_df[(hex_df['180803 48h HD 1h'] == 'yes') & (hex_df['180708 48h'] == 'yes') & (hex_df['180706 6h '] == 'no')]
        qyes_df = hex_df[hex_df['180803 48h HD 1h'] == 'yes?']
        yyy_df = hex_df[(hex_df['180706 6h '] == 'yes')]
        no_no_df = hex_df[hex_df['180708 48h'] == 'no']
        no_orf = list(no_df['ORF'])
        yes_orf = list(yes_df['ORF'])
        qyes_orf = list(qyes_df['ORF'])
        nono_orf = list(no_no_df['ORF'])
        all_orf = list(hex_df['ORF'])
        yyy_orf = list(yyy_df['ORF'])
        no_scores = []
        no_lens = []
        for item in no_orf:
            ndf = lc_df[lc_df['ORF'] == item]
            if len(ndf) > 0:
                lc_score = float(ndf['LC Score'])
                no_scores.append(lc_score)
                no_lens.append(len(str(list(ndf['Sequence'])[0])))
        yes_scores = []
        yes_lens = []
        for item in yes_orf:
            ndf = lc_df[lc_df['ORF'] == item]
            if len(ndf) > 0:
                lc_score = float(ndf['LC Score'])
                yes_scores.append(lc_score)
                yes_lens.append(len(str(list(ndf['Sequence'])[0])))
        print(no_lens)
        print(yes_lens)
        qyes_scores = []
        for item in qyes_orf:
            ndf = lc_df[lc_df['ORF'] == item]
            if len(ndf) > 0:
                lc_score = float(ndf['LC Score'])
                qyes_scores.append(lc_score)
        nono_scores = []
        for item in nono_orf:
            ndf = lc_df[lc_df['ORF'] == item]
            if len(ndf) > 0:
                lc_score = float(ndf['LC Score'])
                nono_scores.append(lc_score)
        all_scores = []
        for item in all_orf:
            ndf = lc_df[lc_df['ORF'] == item]
            if len(ndf) > 0:
                lc_score = float(ndf['LC Score'])
                all_scores.append(lc_score)
        yyy_scores = []
        for item in yyy_orf:
            ndf = lc_df[lc_df['ORF'] == item]
            if len(ndf) > 0:
                lc_score = float(ndf['LC Score'])
                yyy_scores.append(lc_score)
        print(len(yes_scores))
        print(len(no_scores))
        print(np.mean(yes_scores))
        print(np.mean(no_scores))
        #plt.hist(all_scores, bins=20, range=(-60, 200), normed=True)
        plt.hist(yes_lens, bins=20, normed=True, alpha=0.5)
        #plt.hist(nono_scores, bins=10, range=(-60, 200), alpha=0.5)
        #plt.hist(qyes_scores, bins=10)
        plt.hist(no_lens, bins=20, alpha=0.5, normed=True)
        print(yyy_scores)
        plt.show()

    def stubborn_puncta(self):
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        yyy_df = hex_df[(hex_df['180708 48h'] == 'yes') & (hex_df['180706 6h '] == 'yes') & (hex_df['180803 48h HD 1h'] == 'yes')]
        yyy_orf = list(yyy_df['ORF'])
        lc_df = pd.read_csv(self.puncta_fpi, sep='\t', index_col=0)
        ndf = lc_df[lc_df['ORF'].isin(yyy_orf)]
        print(list(ndf['Sequence']))

    def high_scoring_agg(self):
        lc_df = pd.read_csv(self.puncta_fpi, sep='\t', index_col=0)
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        hex_df = hex_df[(hex_df['180708 48h'] == 'yes') & (hex_df['180706 6h '] == 'no')]
        no_df = hex_df[(hex_df['180803 48h HD 1h'] == 'no')]
        no_orf = list(no_df['ORF'])
        yes_df = hex_df[(hex_df['180803 48h HD 1h'] == 'yes')]
        yes_orf = list(yes_df['ORF'])
        no_lc = lc_df[lc_df['ORF'].isin(no_orf)]
        yes_lc = lc_df[lc_df['ORF'].isin(yes_orf)]
        yes_lc = yes_lc[yes_lc['LC Score'] > 0]
        no_lc = no_lc[no_lc['LC Score'] > 0]
        no_seqs = list(no_lc['Sequence'])
        yes_seqs = list(yes_lc['Sequence'])
        for seq in no_seqs:
            analysed_seq = ProteinAnalysis(seq)
            adict = analysed_seq.get_amino_acids_percent()
            qn = adict['Q'] + adict['N']
            if qn > 0.15:
                print(seq)
        print()
        for seq in yes_seqs:
            analysed_seq = ProteinAnalysis(seq)
            adict = analysed_seq.get_amino_acids_percent()
            qn = adict['Q'] + adict['N']
            if qn > 0.15:
                print(seq)

    def ml_approach(self):
        lc_df = pd.read_csv(self.puncta_fpi, sep='\t', index_col=0)
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        hex_df = hex_df[
            (hex_df['180708 48h'] == 'yes') & (hex_df['180706 6h '] == 'no')]
        no_df = hex_df[(hex_df['180803 48h HD 1h'] == 'no')]
        no_orf = list(no_df['ORF'])
        yes_df = hex_df[(hex_df['180803 48h HD 1h'] == 'yes')]
        yes_orf = list(yes_df['ORF'])
        no_lc = lc_df[lc_df['ORF'].isin(no_orf)]
        yes_lc = lc_df[lc_df['ORF'].isin(yes_orf)]
        yes_lc = yes_lc[yes_lc['LC Score'] > 0]
        no_lc = no_lc[no_lc['LC Score'] > 0]
        no_seqs = list(no_lc['Sequence'])
        no_scores = list(no_lc['LC Score'])
        yes_seqs = list(yes_lc['Sequence'])
        yes_scores = list(yes_lc['LC Score'])
        all_vals = {'R': [], 'T': [], 'L': [], 'S': [], 'V': [], 'Y': [],
                    'M': [], 'W': [], 'E': [], 'K': [], 'G': [], 'F': [],
                    'Q': [], 'I': [], 'C': [], 'P': [], 'H': [], 'score': [],
                    'D': [], 'N': [], 'A': []}
        aclass = []
        for seq, score in zip(no_seqs, no_scores):
            analysed_seq = ProteinAnalysis(seq)
            adict = analysed_seq.get_amino_acids_percent()
            adict['score'] = score
            for item in adict:
                all_vals[item].append(adict[item])
            aclass.append(0)
        for seq, score in zip(yes_seqs, yes_scores):
            analysed_seq = ProteinAnalysis(seq)
            adict = analysed_seq.get_amino_acids_percent()
            adict['score'] = score
            for item in adict:
                all_vals[item].append(adict[item])
            aclass.append(1)
        df = pd.DataFrame(all_vals)
        df = df[['Y']]
        #print(df.head())
        #print(df.describe())
        clf = RandomForestClassifier(n_estimators=10, random_state=1)
        clf = clf.fit(df, aclass)
        #yp = clf.predict(df, aclass)
        scores = cross_val_score(clf, df, aclass)
        print(scores)
        print(scores.mean())


class TyrMotifs(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'
        self.lc_m = 0.06744064704548541
        self.lc_b = 16.5
        self.hex_fpi = os.path.join(data_dp, 'experiment', '180803_HD.xls')
        self.puncta_fpi = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')

    def count_tyr(self):
        yes_seqs, no_seqs = self.load_seqs()
        all_tyr = []
        tyr_counts = []
        asp_counts = []
        for seq in yes_seqs:
            tyr_motifs = self.detect_tyr(seq)
            all_tyr.append(tyr_motifs)
            tyr_counts.append(seq.count('Y'))
            asp_counts.append(seq.count('N'))
        print(np.mean(all_tyr))
        print(all_tyr)
        print(np.mean(tyr_counts))
        print(tyr_counts)
        print(np.mean(asp_counts))
        print(asp_counts)
        all_tyr = []
        tyr_counts = []
        asp_counts = []
        for seq in no_seqs:
            tyr_motifs = self.detect_tyr(seq)
            all_tyr.append(tyr_motifs)
            tyr_counts.append(seq.count('Y'))
            asp_counts.append(seq.count('N'))
        print(np.mean(all_tyr))
        print(all_tyr)
        print(np.mean(tyr_counts))
        print(tyr_counts)
        print(np.mean(asp_counts))
        print(asp_counts)

    def detect_tyr(self, seq):
        tyr_motifs = 0
        indexes = tools_lc.lc_to_indexes(seq, self.k, self.lca, self.lce)
        tyr_ind = [pos for pos, char in enumerate(seq) if char == 'Y']
        for i in tyr_ind:
            for j in range(i-2, i+3):
                if j in indexes:
                    tyr_motifs += 1
                    break
        return tyr_motifs

    def load_seqs(self):
        lc_df = pd.read_csv(self.puncta_fpi, sep='\t', index_col=0)
        hex_df = pd.read_excel(self.hex_fpi, sheetname='Hoja2')
        hex_df = hex_df[(hex_df['180708 48h'] == 'yes') & (hex_df['180706 6h '] == 'no')]
        no_df = hex_df[(hex_df['180803 48h HD 1h'] == 'no')]
        no_orf = list(no_df['ORF'])
        yes_df = hex_df[(hex_df['180803 48h HD 1h'] == 'yes')]
        yes_orf = list(yes_df['ORF'])
        no_lc = lc_df[lc_df['ORF'].isin(no_orf)]
        yes_lc = lc_df[lc_df['ORF'].isin(yes_orf)]
        yes_lc = yes_lc[yes_lc['LC Score'] > 100]
        no_lc = no_lc[no_lc['LC Score'] > 100]
        no_seqs = list(no_lc['Sequence'])
        no_scores = list(no_lc['LC Score'])
        yes_seqs = list(yes_lc['Sequence'])
        yes_scores = list(yes_lc['LC Score'])
        seqs = list(yes_lc['Sequence'])
        for seq in seqs:
            print(seq)
        print(yes_lc.head())
        return yes_seqs, no_seqs

    def n_clumps(self, seq):
        pass


def main():
    tm = TyrMotifs()
    tm.load_seqs()



if __name__ == '__main__':
    main()