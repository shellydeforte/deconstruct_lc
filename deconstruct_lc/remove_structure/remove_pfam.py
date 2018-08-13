import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc
from deconstruct_lc.scores.norm_score import NormScore

class RemovePfam(object):
    def __init__(self):
        config = read_config.read_config()
        self.data_dp = os.path.join(config['fps']['data_dp'])
        self.puncta = os.path.join(self.data_dp, 'experiment', 'puncta_uni.fasta')
        self.nopuncta = os.path.join(self.data_dp, 'experiment', 'nopuncta_uni.fasta')
        self.pfam_puncta = os.path.join(self.data_dp, 'experiment', 'puncta_pfam.tsv')
        self.pfam_nopuncta = os.path.join(self.data_dp, 'experiment', 'nopuncta_pfam.tsv')
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'
        self.lc_m = 0.06744064704548541
        self.lc_b = 16.5

    def run_percent_pfam(self):
        puncta_perc = os.path.join(self.data_dp, 'experiment', 'puncta_percent_pfam.tsv')
        self.percent_pfam(self.puncta, self.pfam_puncta, puncta_perc)
        nopuncta_perc = os.path.join(self.data_dp, 'experiment', 'nopuncta_percent_pfam.tsv')
        self.percent_pfam(self.nopuncta, self.pfam_nopuncta, nopuncta_perc)

    def percent_pfam(self, fasta_fp, pfam_fp, fpo):
        df = pd.read_csv(pfam_fp, sep='\t')
        pids, seqs = tools_fasta.fasta_to_id_seq(fasta_fp)
        frac_pfam = []
        for id, seq in zip(pids, seqs):
            ndf = df[df['uniprot_acc'] == id]
            ndf = ndf.sort_values(by='seq_start')
            segmented = self.segment_seq(seq, ndf)
            len_seg = 0
            for seg in segmented:
                len_seg += len(seg)
            frac_pfam.append(float(len(seq) - len_seg)/float(len(seq)))
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        df_out = pd.DataFrame({'Uniprot ID': pids, 'LC Score': scores,
                               'Pfam Fraction': frac_pfam}, columns=['Uniprot ID', 'LC Score', 'Pfam Fraction'])
        df_out = df_out.sort_values(by='LC Score', ascending=False)
        df_out.to_csv(fpo, sep='\t')
        print(np.mean(frac_pfam))

    def run_with_pfam(self):
        puncta_out = os.path.join(self.data_dp, 'experiment', 'puncta_nopfam.tsv')
        self.with_pfam(self.puncta, self.pfam_puncta, puncta_out)
        nopuncta_out = os.path.join(self.data_dp, 'experiment', 'nopuncta_nopfam.tsv')
        self.with_pfam(self.nopuncta, self.pfam_nopuncta, nopuncta_out)

    def with_pfam(self, fasta_fp, pfam_fp, fpo):
        """
        How many proteins in the set have pfam domains?
        What is the fraction occupied by pfam domains?"""
        df = pd.read_csv(pfam_fp, sep='\t')
        pfam_ids = list(set(df['uniprot_acc']))
        pids, seqs = tools_fasta.fasta_to_id_seq(fasta_fp)
        print(len(pids))
        nopfam_ids = list(set(pids) - set(pfam_ids))
        nopfam_seqs = []
        for pid, seq in zip(pids, seqs):
            if pid in nopfam_ids:
                nopfam_seqs.append(seq)
        ns = NormScore()
        scores = ns.lc_norm_score(nopfam_seqs)
        df_out = pd.DataFrame({'UniProt ID': nopfam_ids, 'LC Score': scores}, columns=['UniProt ID', 'LC Score'])
        df_out = df_out.sort_values(by='LC Score', ascending=False)
        df_out.to_csv(fpo, sep='\t')

    def fetch_score(self, df, pids):
        scores = []
        for pid in pids:
            df = df[df['Protein ID'] == pid]
            scores.append(list(df['LC Score'])[0])
        return scores

    def score_in_pfam(self):
        ids, seqs = tools_fasta.fasta_to_id_seq(self.nopuncta)
        df = pd.read_csv(self.pfam_nopuncta, sep='\t', index_col=0)
        below = 0
        above = 0
        norm_scores = []
        fl_norm_scores = []
        for id, seq in zip(ids, seqs):
            ndf = df[df['uniprot_acc'] == id]
            ndf = ndf.sort_values(by='seq_start')
            segmented = self.pfam_segments(seq, ndf)
            total = 0
            for item in segmented:
                total += len(item)
            if total >= 100:
                above += 1
                fl_score, fl_length = self.get_segment_scores([seq])
                fl_norm = self.norm_function([fl_score], [fl_length])
                raw_score, length = self.get_segment_scores(segmented)
                norm_score = self.norm_function([raw_score], [length])
                norm_scores.append(norm_score[0])
                fl_norm_scores.append(fl_norm[0])
            else:
                below += 1
        print(above)
        print(below)
        print(np.mean(norm_scores))
        print(np.mean(fl_norm_scores))
        print(np.median(norm_scores))
        print(np.median(fl_norm_scores))
        plt.hist(fl_norm_scores, alpha=0.5, bins=20, range=(-100, 200), label='Full length scores')
        plt.hist(norm_scores, alpha=0.5, bins=20, range=(-100, 200), label='Inside Pfam scores')
        plt.legend()
        plt.show()

    def run(self):
        ids, seqs = tools_fasta.fasta_to_id_seq(self.puncta)
        df = pd.read_csv(self.pfam_puncta, sep='\t', index_col=0)
        new_seqs = []
        below = 0
        above = 0
        norm_scores = []
        fl_norm_scores = []
        for id, seq in zip(ids, seqs):
            ndf = df[df['uniprot_acc'] == id]
            ndf = ndf.sort_values(by='seq_start')
            segmented = self.segment_seq(seq, ndf)
            total = 0
            for item in segmented:
                total += len(item)
            if total >= 100:
                above += 1
                fl_score, fl_length = self.get_segment_scores([seq])
                fl_norm = self.norm_function([fl_score], [fl_length])
                raw_score, length = self.get_segment_scores(segmented)
                norm_score = self.norm_function([raw_score], [length])
                norm_scores.append(norm_score[0])
                fl_norm_scores.append(fl_norm[0])
            else:
                below += 1
        print(above)
        print(below)
        print(np.mean(norm_scores))
        print(np.mean(fl_norm_scores))
        print(np.median(norm_scores))
        print(np.median(fl_norm_scores))
        plt.hist(fl_norm_scores, alpha=0.5, bins=20, range=(-100, 200), label='Full length scores')
        plt.hist(norm_scores, alpha=0.5, bins=20, range=(-100, 200), label='Outside Pfam scores')
        plt.legend()
        plt.show()

    def pfam_segments(self, seq, df):
        new_seq = []
        for i, row in df.iterrows():
            new_seq.append(seq[row['seq_start']: row['seq_end']+1])
        return new_seq

    def segment_seq(self, seq, df):
        """Given intervals, pull out the domain, and segment around it"""
        start = 0
        new_seq = []
        for i, row in df.iterrows():
            new_seq.append(seq[start:row['seq_start']])
            start = row['seq_end'] + 1
        new_seq.append(seq[start:])
        return new_seq

    def pfam_in_common(self):
        df = pd.read_csv(self.pfam_puncta, sep='\t', index_col=0)
        print(df['pfamA_acc'].value_counts())

    def get_segment_scores(self, segment_seq):
        total_motifs = 0
        total_length = 0
        for seq in segment_seq:
            motifs = tools_lc.count_lc_motifs(seq, self.k, self.lca, self.lce)
            total_motifs += motifs
            total_length += len(seq)
        return total_motifs, total_length

    def norm_function(self, raw_scores, lengths):
        norm_scores = []
        for raw_score, length in zip(raw_scores, lengths):
            norm_score = raw_score - ((self.lc_m * length) + self.lc_b)
            norm_scores.append(norm_score)
        return norm_scores


def main():
    rp = RemovePfam()
    rp.pfam_in_common()


if __name__ == '__main__':
    main()