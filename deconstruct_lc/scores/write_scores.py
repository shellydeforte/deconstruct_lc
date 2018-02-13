import os
import pandas as pd
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc.analysis_bc.write_bc_score import BcScore
from deconstruct_lc import read_config


class WriteNorm(object):
    def __init__(self):
        config = read_config.read_config()
        data = config['fps']['data_dp']
        prot_dp = os.path.join(data, 'proteomes')
        self.bc_dp = os.path.join(data, 'bc_analysis')

        self.train_fpi = os.path.join(data, 'train.tsv')
        self.yeast_fp = os.path.join(prot_dp, 'UP000002311_559292_Yeast.fasta')
        self.human_fp = os.path.join(prot_dp, 'UP000005640_9606_Human.fasta')
        self.fpo = os.path.join(data, 'scores', 'pdb_bc_scores.tsv')

    def write_scores(self):
        """Write tsv that is pid, proteome, org, lc score"""
        pdb_pids, pdb_proteome, pdb_org, pdb_scores = self.get_pdb()
        bc_pids, bc_proteome, bc_org, bc_scores = self.get_bcs()
        df_dict = {'Protein ID': pdb_pids + bc_pids,
                   'Proteome': pdb_proteome + bc_proteome,
                   'Organism': pdb_org + bc_org,
                   'LC Score': pdb_scores + bc_scores}
        cols = ['Protein ID', 'Proteome', 'Organism', 'LC Score']
        df_out = pd.DataFrame(df_dict, columns=cols)
        df_out.to_csv(self.fpo, sep='\t')

    def get_pdb(self):
        df = pd.read_csv(self.train_fpi, sep='\t', index_col=0)
        pdb_df = df[df['y'] == 1]
        pids = list(pdb_df['Protein ID'])
        proteome = ['PDB']*len(pids)
        org = ['PDB']*len(pids)
        seqs = list(pdb_df['Sequence'])
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        return pids, proteome, org, scores

    def get_bcs(self):
        bcs = BcScore()
        bc_names = bcs.get_sheets()
        pids = []
        proteome = []
        org = []
        scores = []
        for bc_name in bc_names:
            bc_fp = os.path.join(self.bc_dp, '{}_score.tsv'.format(bc_name))
            df = pd.read_csv(bc_fp, sep='\t', index_col=0)
            pids += list(df['Protein ID'])
            proteome += [bc_name]*len(list(df['Protein ID']))
            org += list(df['Organism'])
            scores += list(df['LC Score'])
        return pids, proteome, org, scores


def main():
    wn = WriteNorm()
    wn.write_scores()


if __name__ == '__main__':
    main()