"""
1. Create excel with Protein ID, PubMed ID, source (Manual if manual).
2. Convert yeast IDs
"""
import os
import pandas as pd
from Bio import SeqIO
from predict_llps import tools_fasta

class LoadGO(object):

    def __init__(self):
        self.fd = os.path.join(os.path.dirname(__file__), '..', 'data', 'quickgo')
        self.cb_fp = os.path.join(self.fd, 'quickgo_cb.xlsx')
        self.yeast_pid_fp = os.path.join(self.fd, 'yeast_pids.txt')
        self.fasta_in = os.path.join(self.fd, 'quickgo_cb.fasta')
        self.yeast_map_fp = os.path.join(self.fd, 'yeast_map.xlsx')
        self.pids_fp = os.path.join(self.fd, 'pids.txt')

    def write_go(self):
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

    def write_yeast_pids(self):
        fns = ['Cajal_bodies', 'Centrosome', 'P_granule', 'Nuclear_Speckles', 'Nucleolus',
                    'P_Body', 'PML_Body', 'Paraspeckle']
        yeast_pids = []
        for sheet in fns:
            df_in = pd.read_excel(self.new_cb_fp, sheetname=sheet)
            yeast_df_in = df_in[df_in['Organism'] == 'YEAST']
            yeast_pids += list(yeast_df_in['Protein ID'])
        with open(self.yeast_pid_fp, 'w') as fo:
            for pid in set(yeast_pids):
                fo.write('{}\n'.format(pid))

    def create_yeast_dict(self):
        """
        Yeast map comes from here http://www.uniprot.org/docs/yeast.txt, imported into excel, preserve 2 cols
        checked one entry by hand in the returned dictionary
        """
        yeast_map = {}
        df = pd.read_excel(self.yeast_map_fp, sheetname='Sheet1', header=None)
        for i, row in df.iterrows():
            yeast_map[row[1]] = row[0]
        return yeast_map

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

    def compare_cb(self):
        """Not using"""
        sheet_fn = {'Cajal Bodies': 'Cajal_bodies', 'Centrosomes': 'Centrosome',
                    'Germinal Granules': 'P_granule', 'Nuclear Speckle': 'Nuclear_Speckles',
                    'Nucleolus': 'Nucleolus', 'P bodies': 'P_Body', 'PML body': 'PML_Body',
                    'Paraspeckle': 'Paraspeckle'}
        for sheet in sheet_fn:
            yeast_map = self.create_yeast_dict()
            cb_df = pd.read_excel(self.cb_fp, sheetname=sheet)
            cb_ids = set(list(cb_df['Protein ID']))
            go_fp = os.path.join(self.fd, '{}.tsv'.format(sheet_fn[sheet]))
            go_df = pd.read_csv(go_fp, sep='\t', comment='!', header=None)
            go_ids = set(list(go_df[1]))
            new_go_ids = []
            for pid in go_ids:
                if pid in yeast_map:
                    new_go_ids.append(yeast_map[pid])
                else:
                    new_go_ids.append(pid)
            new_go_ids = set(new_go_ids)
            print(sheet)
            print(cb_ids - new_go_ids)
            print(len(cb_ids - new_go_ids))


def main():
    lg = LoadGO()
    pids = lg.write_go()


if __name__ == '__main__':
    main()