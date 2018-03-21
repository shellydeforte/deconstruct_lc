from datetime import datetime
import os
import pandas as pd
from Bio import SeqIO

from deconstruct_lc import read_config
from deconstruct_lc.data_bc import pull_uni


class WriteGO(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.minlen = config['dataprep'].getint('minlen')
        self.maxlen = config['dataprep'].getint('maxlen')
        self.fd = os.path.join(data_dp, 'bc_prep')
        self.now = datetime.now().strftime("%y%m%d")
        self.cb_fp = os.path.join(self.fd, '{}quickgo_bc.xlsx'.format(
            self.now))
        self.fasta_in = os.path.join(self.fd, 'quickgo_bc.fasta')
        self.fasta_len = os.path.join(self.fd, 'quickgo_bc_len.fasta')
        self.pids_fp = os.path.join(self.fd, '{}pids.txt'.format(self.now))
        # Alternatively spliced proteins must be dealt with separately
        self.pids_alt_fp = os.path.join(self.fd, '{}pids_alt.txt'.format(
            self.now))
        self.alt_fasta = os.path.join(self.fd,
                                      '{}quickgo_bc_alt.fasta'.format(self.now))

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
                fgo_df = go_df[go_df[1] == go_id]
                if go_id not in pid_gene_org:
                    go_id = go_id.split('-')[0]
                gene = pid_gene_org[go_id][0]
                org = pid_gene_org[go_id][1]
                pmid = list(fgo_df[4])[0]
                source = list(fgo_df[9])[0]
                df_dict['Protein ID'].append(go_id)
                df_dict['Reference'].append(pmid)
                df_dict['Source'].append(source)
                df_dict['Gene ID'].append(gene)
                df_dict['Organism'].append(org)

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
        with open(self.pids_fp, 'w') as fo, open(self.pids_alt_fp, 'w') as fao:
            for pid in pids:
                if '-' in pid:
                    fao.write('{}\n'.format(pid))
                else:
                    fo.write('{}\n'.format(pid))

    def read_pids(self, fp):
        pids = []
        with open(fp, 'r') as fpi:
            for line in fpi:
                pids.append(line.strip())
        return pids

    def get_sheets(self):
        ex = pd.ExcelFile(self.cb_fp)
        sheet_names = ex.sheet_names
        return sorted(sheet_names)

    def filter_fasta(self):
        """Filter CB fasta for length"""
        new_records = []
        with open(self.fasta_in, 'r') as cb_in:
            for seq_rec in SeqIO.parse(cb_in, 'fasta'):
                sequence = str(seq_rec.seq)
                prot_len = len(sequence)
                if self.minlen <= prot_len <= self.maxlen:
                    if self.standard_aa(sequence):
                        new_records.append(seq_rec)
        with open(self.fasta_len, 'w') as seq_fo:
            SeqIO.write(new_records, seq_fo, 'fasta')
        count = 0
        with open(self.fasta_len, 'r') as handle:
            for _ in SeqIO.parse(handle, 'fasta'):
                count += 1
        print('There are {} records'.format(count))

    def standard_aa(self, sequence):
        aas = 'ADKERNTSQYFLIVMCWHGP'
        for c in sequence:
            if c not in aas:
                return False
        return True


def main():
    lg = WriteGO()
    pids = lg.get_pids_from_qg()
    lg.write_pids(pids)
    alt_pids = lg.read_pids(lg.pids_alt_fp)
    pull_uni.write_fasta(alt_pids, lg.alt_fasta)
    ###########################################################################
    # Here the PID list must be manually uploaded to uniprot to get the       #
    # fasta file and then concatenated with the alt pids                      #
    # Do this before creating the spreadsheet and filtering the fasta file    #
    ###########################################################################
    # lg.go_to_ss()
    # lg.filter_fasta()


if __name__ == '__main__':
    main()