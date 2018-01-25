from collections import defaultdict
import configparser
import os
import pandas as pd
from deconstruct_lc import tools_fasta

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
        self.exp_fpo = os.path.join(self.exp_dp, 'bc_pdb.tsv')
        self.bc_fp = os.path.join(self.bc_dp, 'quickgo_bc.xlsx')
        self.bc_fasta = os.path.join(self.bc_dp, 'quickgo_bc_len.fasta')
        self.pdb_uni_fp = os.path.join(self.pdb_dp, 'pdb_chain_uniprot.tsv')
        self.pdb_an_fp = os.path.join(self.pdb_dp, 'pdb_analysis.tsv')

    def write_uni_pdb_org_comp(self):
        """Write a tsv with pdb_chains of BC proteins"""
        df_dict = {'PDB ID': [], 'Uni ID': [], 'PDB Seq': [], 'Uni Seq': [],
                   'Miss Seq': [], 'Organism': [], 'Compartment': []}
        pdb_uni = self.get_pdb_uni()
        uni_seq = self.get_uni_seq()
        uni_comp, uni_org = self.get_uni_comp_org()
        pdb_an_df = pd.read_csv(self.pdb_an_fp, sep='\t', index_col=0)
        for i, row in pdb_an_df.iterrows():
            pdb_id = row['Protein ID']
            if pdb_id in pdb_uni:
                uni = pdb_uni[pdb_id]
                if pdb_uni[pdb_id] in uni_seq:
                    org = uni_org[uni]
                    comp = uni_comp[uni]
                    df_dict['PDB ID'].append(pdb_id)
                    df_dict['Uni ID'].append(uni)
                    df_dict['Uni Seq'].append(uni_seq[uni])
                    df_dict['PDB Seq'].append(row['Sequence'])
                    df_dict['Miss Seq'].append(row['Missing'])
                    df_dict['Organism'].append(org)
                    df_dict['Compartment'].append(comp)
        cols = ['PDB ID', 'Uni ID', 'Organism', 'Compartment', 'Uni Seq',
                'PDB Seq', 'Miss Seq']
        df_out = pd.DataFrame(df_dict, columns=cols)
        df_out.to_csv(self.exp_fpo, sep='\t')

    def get_uni_seq(self):
        uni_seq = {}
        pids, seqs = tools_fasta.fasta_to_id_seq(self.bc_fasta)
        for pid, seq in zip(pids, seqs):
            uni_seq[pid] = seq
        return uni_seq

    def get_pdb_uni(self):
        """If Uni ID is in BC, then create pdb_chain: uni dict."""
        df = pd.read_csv(self.pdb_uni_fp, sep='\t', header=1)
        uni_comp, uni_org = self.get_uni_comp_org()
        pdb_uni = {}
        for i, row in df.iterrows():
            pid = row['SP_PRIMARY']
            if pid in uni_comp:
                try:
                    pdb_id = row['PDB'].upper()
                    chain = row['CHAIN'].upper()
                    pdb_chain = '{}_{}'.format(pdb_id, chain)
                    pdb_uni[pdb_chain] = pid
                except:
                    print(row)
        return pdb_uni

    def get_uni_comp_org(self):
        bcs = self._get_sheets(self.bc_fp)
        uni_comp = defaultdict(list)
        uni_org = {}
        for sheet in bcs:
            df_in = pd.read_excel(self.bc_fp, sheetname=sheet)
            for i, row in df_in.iterrows():
                pid = row['Protein ID']
                uni_comp[pid].append(sheet)
                org = row['Organism']
                uni_org[pid] = org
        return uni_comp, uni_org

    def _get_sheets(self, fp):
        ex = pd.ExcelFile(fp)
        sheet_names = ex.sheet_names
        return sorted(sheet_names)


def main():
    pipe = UniPdb()
    pipe.write_uni_pdb_org_comp()


if __name__ == '__main__':
    main()


