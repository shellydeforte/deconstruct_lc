from Bio import pairwise2
from Bio.pairwise2 import format_alignment
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
        self.exp_fpo = os.path.join(self.exp_dp, 'bc_pdb_align.txt')
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'

    def read_pdb(self):
        df = pd.read_csv(self.exp_fpi, sep='\t', index_col=0)
        df = df[df['Organism'] == 'YEAST']
        df = df[(df['Compartment'] == '[\'Cytoplasmic_Stress_Granule\']') |
                (df['Compartment'] == '[\'P_Body\']') |
                (df['Compartment'] == '[\'Cytoplasmic_Stress_Granule\', \'P_Body\']')]
        with open(self.exp_fpo, 'w') as fo:
            for i, row in df.iterrows():
                pdb_seq = row['PDB Seq']
                uni_seq = row['Uni Seq']
                uni_id = row['Uni ID']
                pdb_id = row['PDB ID']
                comp = row['Compartment']
                uni_ns = NormScore([row['Uni Seq']])
                uni_score = uni_ns.lc_norm_score()[0]
                pdb_ns = NormScore([row['PDB Seq']])
                pdb_score = pdb_ns.lc_norm_score()[0]
                mot_seq = tools_lc.display_lc(uni_seq, self.k, self.lca,
                                              self.lce)
                if uni_score > 5 and pdb_score < -5:
                    alignments = pairwise2.align.localms(uni_seq, pdb_seq, 1,-1,-2,-2)
                    als = format_alignment(*alignments[0])
                    fo.write('{}\t{}\t{}\t{}\t{}\n'.format(uni_id, pdb_id,
                                                           comp, uni_score,
                                                           pdb_score))
                    fo.write(mot_seq)
                    fo.write('\n')
                    fo.write(als)
                    fo.write('\n')



def main():
    up = UniPdbScore()
    up.read_pdb()


if __name__ == '__main__':
    main()