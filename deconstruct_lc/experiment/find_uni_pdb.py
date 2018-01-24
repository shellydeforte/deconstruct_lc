import configparser
import os
import pandas as pd
from deconstruct_lc.data_pdb import pdb
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores import calc_score

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class UniPdb(object):
    def __init__(self):
        self.bc_dp = os.path.join(config['filepaths']['data_dp'], 'bc_prep')
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.exp_dp = os.path.join(config['filepaths']['data_dp'],
                                   'experiment')
        self.exp_fpo = os.path.join(self.exp_dp, 'pdb_yeast.tsv')
        self.bc_fp = os.path.join(self.bc_dp, 'quickgo_bc.xlsx')
        self.bc_fasta = os.path.join(self.bc_dp, 'quickgo_bc_len.fasta')
        self.pdb_uni_fp = os.path.join(self.pdb_dp, 'pdb_chain_uniprot.tsv')
        self.pdb_an_fp = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')

    def get_pdb_uni(self):
        df = pd.read_csv(self.pdb_uni_fp, sep='\t', header=1)
        bc_yeast, all_pids = self.get_pids_from_cb()
        pdb_uni = {}
        for i, row in df.iterrows():
            pid = row['SP_PRIMARY']
            if pid in all_pids:
                print(pid)
                try:
                    pdb_id = row['PDB'].upper()
                    chain = row['CHAIN'].upper()
                    pdb_chain = '{}_{}'.format(pdb_id, chain)
                    pdb_uni[pdb_chain] = pid
                except:
                    pass
        return pdb_uni

    def get_seq(self):
        df_dict = {'PDB ID': [], 'Uni ID': [], 'PDB Seq': [], 'Uni Seq': [],
                   'Miss Seq': [], 'PDB Norm': [], 'Uni Norm': []}
        pdb_uni = self.get_pdb_uni()
        pdb_an_df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        uni_seq = self.get_uni_seq()
        for i, row in pdb_an_df.iterrows():
            pdb_id = row['Protein ID']
            if pdb_id in pdb_uni:
                uni = pdb_uni[pdb_id]
                print(uni)
                if pdb_uni[pdb_id] in uni_seq:
                    df_dict['PDB ID'].append(pdb_id)
                    df_dict['Uni ID'].append(uni)
                    df_dict['Uni Seq'].append(uni_seq[uni])
                    df_dict['PDB Seq'].append(row['Sequence'])
                    df_dict['Miss Seq'].append(row['Missing'])
                    uni_cs = calc_score.CalcScore([uni_seq[uni]])
                    uni_norm = uni_cs.lc_norm_score()
                    df_dict['Uni Norm'].append(uni_norm)
                    pdb_cs = calc_score.CalcScore([row['Sequence']])
                    pdb_norm = pdb_cs.lc_norm_score()
                    df_dict['PDB Norm'].append(pdb_norm)
        cols = ['PDB ID', 'Uni ID', 'PDB Norm', 'Uni Norm', 'Uni Seq',
                'PDB Seq', 'Miss Seq']
        df_out = pd.DataFrame(df_dict, columns=cols)
        df_out.to_csv(self.exp_fpo, sep='\t')

    def get_uni_seq(self):
        uni_seq = {}
        pids, seqs = tools_fasta.fasta_to_id_seq(self.bc_fasta)
        for pid, seq in zip(pids, seqs):
            uni_seq[pid] = seq
        return uni_seq

    def get_sheets(self, fp):
        ex = pd.ExcelFile(fp)
        sheet_names = ex.sheet_names
        return sorted(sheet_names)

    def get_pids_from_cb(self):
        fns = self.get_sheets(self.bc_fp)
        bc_yeast = {}
        all_pids = []
        for sheet in fns:
            df_in = pd.read_excel(self.bc_fp, sheetname=sheet)
            df_yeast = df_in[df_in['Organism'] == 'YEAST']
            pids = list(df_yeast['Protein ID'])
            bc_yeast[sheet] = pids
            all_pids += pids
        return bc_yeast, set(all_pids)

    def check_yeast_pdb(self):
        df = pd.read_csv(self.exp_fpo, sep='\t', index_col=0)
        df = df[['PDB ID', 'Uni ID', 'PDB Norm', 'Uni Norm']]
        pids = set(list(df['Uni ID']))
        print(len(pids))
        for pid in pids:
            ndf = df[df['Uni ID'] == pid]
            print(ndf)



def main():
    pipe = UniPdb()
    pipe.check_yeast_pdb()


if __name__ == '__main__':
    main()


