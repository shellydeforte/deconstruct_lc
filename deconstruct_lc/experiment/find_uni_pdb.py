import configparser
import os
import pandas as pd

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class UniPdb(object):
    def __init__(self):
        self.bc_dp = os.path.join(config['filepaths']['data_dp'], 'bc_prep')
        self.bc_fp = os.path.join(self.bc_dp, 'quickgo_bc.xlsx')
        self.pdb_tax_fp = os.path.join(config['filepaths']['data_dp'],
                                       'pdb_prep', 'pdb_chain_uniprot.tsv')

    def get_yeast_uni(self):
        bc_yeast = self.get_pids_from_cb()
        pdb_uni = self.get_pdb_uni()
        for bc in bc_yeast:
            print(bc)
            bc_pids = set(bc_yeast[bc])
            print(pdb_uni & bc_pids)

    def get_pdb_uni(self):
        df = pd.read_csv(self.pdb_tax_fp, sep='\t', header=1)
        pids = set(list(df['SP_PRIMARY']))
        return pids

    def get_sheets(self, fp):
        ex = pd.ExcelFile(fp)
        sheet_names = ex.sheet_names
        return sorted(sheet_names)

    def get_pids_from_cb(self):
        fns = self.get_sheets(self.bc_fp)
        bc_yeast = {}
        for sheet in fns:
            df_in = pd.read_excel(self.bc_fp, sheetname=sheet)
            df_yeast = df_in[df_in['Organism'] == 'YEAST']
            pids = list(df_yeast['Protein ID'])
            bc_yeast[sheet] = pids
        return bc_yeast


def main():
    pipe = UniPdb()
    pipe.get_yeast_uni()


if __name__ == '__main__':
    main()


