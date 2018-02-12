import os
import pandas as pd
from Bio import SeqIO
from deconstruct_lc import read_config
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc import tools_fasta


class BcScore(object):
    def __init__(self):
        self.config = read_config.read_config()
        self.data_dp = self.config['fps']['data_dp']
        self.bc_dp = os.path.join(self.data_dp, 'bc_prep')
        self.bc_an_dp = os.path.join(self.data_dp, 'bc_analysis')
        # Use fasta file with all bc sequences
        self.fasta = os.path.join(self.bc_dp, 'quickgo_bc.fasta')
        self.bc_ss = os.path.join(self.bc_dp, 'quickgo_bc.xlsx')
        self.bc_score_fp = os.path.join(self.bc_an_dp, 'bc_all_score.tsv')

    def compile_bcs(self):
        bc_pids = self.create_bc_dict()
        df_in = pd.read_csv(self.bc_score_fp, sep='\t', index_col=0)
        for bc in bc_pids:
            fno = '{}_score.tsv'.format(bc)
            fpo = os.path.join(self.bc_an_dp, fno)
            pids = bc_pids[bc]
            ndf = df_in[df_in['Protein ID'].isin(pids)]
            ndf.to_csv(fpo, sep='\t')

    def create_bc_dict(self):
        fns = self.get_sheets()
        bc_pids = {}
        for sheet in fns:
            df_in = pd.read_excel(self.bc_ss, sheetname=sheet)
            bc_pids[sheet] = list(df_in['Protein ID'])
        return bc_pids

    def get_sheets(self):
        ex = pd.ExcelFile(self.bc_ss)
        sheet_names = ex.sheet_names
        return sorted(sheet_names)

    def write_scores(self):
        """Write ID, length, score from fasta file"""
        df_dict = {'Protein ID': [], 'Sequence': [], 'Organism': []}
        with open(self.fasta, 'r') as fi:
            for record in SeqIO.parse(fi, 'fasta'):
                rec_id = record.id.split('|')
                pid = rec_id[1]
                gene_org = rec_id[2]
                org = gene_org.split('_')[1]
                seq = str(record.seq)
                df_dict['Protein ID'].append(pid)
                df_dict['Sequence'].append(seq)
                df_dict['Organism'].append(org)
        seqs = df_dict['Sequence']
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        lengths = tools_fasta.get_lengths(seqs)
        df_dict['LC Score'] = scores
        df_dict['Length'] = lengths
        cols = ['Protein ID', 'Organism', 'Length', 'LC Score', 'Sequence']
        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(self.bc_score_fp, sep='\t')

def main():
    bc = BcScore()
    bc.compile_bcs()


if __name__ == '__main__':
    main()