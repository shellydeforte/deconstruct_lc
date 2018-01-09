import configparser
import os
import pandas as pd
from Bio import SeqIO

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class WriteGO(object):

    def __init__(self):
        self.fd = os.path.join(config['filepaths']['data_fp'], 'bc_prep')
        self.cb_fp = os.path.join(self.fd, 'quickgo_bc.xlsx')
        self.fasta_in = os.path.join(self.fd, 'quickgo_bc.fasta')
        self.pids_fp = os.path.join(self.fd, 'pids.txt')

    def go_to_ss(self):
        pid_gene_org = self.create_org_dict()
        writer = pd.ExcelWriter(self.cb_fp, engine='xlsxwriter')
        fns = ['Cajal_bodies', 'Centrosome', 'Cytoplasmic_Stress_Granule',
               'Nuclear_Speckles', 'Nuclear_Stress_Granule', 'Nucleolus',
               'P_Body', 'P_granule', 'Paraspeckle', 'PML_Body']
        for sheet in fns:
            df_dict = {'Protein ID': [], 'Reference': [], 'Source': [], 'Gene ID': [], 'Organism': []}
            go_fp = os.path.join(self.fd, '{}.tsv'.format(sheet))
            go_df = pd.read_csv(go_fp, sep='\t', comment='!', header=None)
            go_ids = set(list(go_df[1]))
            for go_id in go_ids: # Take just the first entry info
                if go_id in pid_gene_org:
                    gene = pid_gene_org[go_id][0]
                    org = pid_gene_org[go_id][1]
                    fgo_df = go_df[go_df[1] == go_id]
                    pmid = list(fgo_df[4])[0]
                    source = list(fgo_df[9])[0]
                    df_dict['Protein ID'].append(go_id)
                    df_dict['Reference'].append(pmid)
                    df_dict['Source'].append(source)
                    df_dict['Gene ID'].append(gene)
                    df_dict['Organism'].append(org)
                else:
                    print(go_id)
            df_out = pd.DataFrame(df_dict, columns=['Protein ID', 'Gene ID', 'Organism', 'Reference', 'Source'])
            df_out.to_excel(writer, sheet_name=sheet, index=False)

    def create_org_dict(self):
        pid_gene_org = {}
        with open(self.fasta_in, 'r') as fasta_in:
            for record in SeqIO.parse(fasta_in, 'fasta'):
                rec_id = record.id.split('|')
                pid = rec_id[1]
                gene_org = rec_id[2]
                gene = gene_org.split('_')[0]
                org = gene_org.split('_')[1]
                pid_gene_org[pid] = (gene, org)
        return pid_gene_org

    def get_pids_from_cb(self):
        fns = self.get_sheets()
        all_pids = []
        for sheet in fns:
            df_in = pd.read_excel(self.cb_fp, sheetname=sheet)
            all_pids += list(df_in['Protein ID'])
        return set(all_pids)

    def get_pids_from_qg(self):
        all_pids = []
        fns = ['Cajal_bodies', 'Centrosome', 'Cytoplasmic_Stress_Granule',
               'Nuclear_Speckles', 'Nuclear_Stress_Granule', 'Nucleolus',
               'P_Body', 'P_granule', 'Paraspeckle', 'PML_Body']
        for fn in fns:
            go_fp = os.path.join(self.fd, '{}.tsv'.format(fn))
            go_df = pd.read_csv(go_fp, sep='\t', comment='!', header=None)
            go_ids = set(list(go_df[1]))
            all_pids += go_ids
        return set(all_pids)

    def write_pids(self, pids):
        with open(self.pids_fp, 'w') as fo:
            for pid in pids:
                fo.write('{}\n'.format(pid))

    def get_sheets(self):
        ex = pd.ExcelFile(self.cb_fp)
        sheet_names = ex.sheet_names
        return sorted(sheet_names)


def main():
    lg = WriteGO()
    lg.go_to_ss()


if __name__ == '__main__':
    main()