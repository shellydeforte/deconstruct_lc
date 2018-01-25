from collections import defaultdict
import configparser
import os
import pandas as pd
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc
from deconstruct_lc.scores.norm_score import NormScore

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class UniPdbScore(object):
    def __init__(self):
        self.exp_dp = os.path.join(config['filepaths']['data_dp'],
                                   'experiment')
        self.exp_fpi = os.path.join(self.exp_dp, 'bc_pdb.tsv')
        self.exp_fpo = os.path.join(self.exp_dp, 'bc_pdb_scores.tsv')
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'

    def write_scores(self):
        df = pd.read_csv(self.exp_fpi, sep='\t', index_col=0)
        df_dict = {'PDB ID': [], 'Uni ID': [], 'Organism': [],
                   'Compartment': [], 'Uni Score': [], 'PDB Score': [],
                   'Miss Score': []}
        for i, row in df.iterrows():
            df_dict['PDB ID'].append(row['PDB ID'])
            df_dict['Uni ID'].append(row['Uni ID'])
            df_dict['Organism'].append(row['Organism'])
            df_dict['Compartment'].append(row['Compartment'])
            uni_ns = NormScore([row['Uni Seq']])
            uni_score = uni_ns.lc_norm_score()
            pdb_ns = NormScore([row['PDB Seq']])
            pdb_score = pdb_ns.lc_norm_score()
            miss_ns = NormScore([row['PDB Seq']])
            miss_score = miss_ns.lc_miss_norm([row['Miss Seq']])
            df_dict['Uni Score'].append(uni_score)
            df_dict['PDB Score'].append(pdb_score)
            df_dict['Miss Score'].append(miss_score)
        cols = ['PDB ID', 'Uni ID', 'Organism', 'Compartment', 'Uni Score',
                'PDB Score', 'Miss Score']
        df_out = pd.DataFrame(df_dict, columns=cols)
        df_out.to_csv(self.exp_fpo, sep='\t')

    def miss_score(self, pdb_seq, miss_seq):
        pdb_kmers = tools_lc.seq_to_kmers(pdb_seq, self.k)
        miss_kmers = tools_lc.seq_to_kmers(miss_seq, self.k)
        motif_count = 0
        for pdb, miss in zip(pdb_kmers, miss_kmers):
            if miss.count('X') == 0:
                if tools_lc.lca_motif(pdb, self.lca):
                    motif_count += 1
                elif tools_lc.lce_motif(pdb, self.lce):
                    motif_count += 1
                else:
                    pass


def main():
    pipe = UniPdbScore()
    pipe.write_scores()


if __name__ == '__main__':
    main()