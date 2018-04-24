import os
import pandas as pd
from Bio import SeqIO

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.display.display_lc import Display


class Atg(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.atg_fp = os.path.join(data_dp, 'atg', 'atg.xlsx')
        self.atg_out = os.path.join(data_dp, 'atg', 'atg_gene_orf_seq.tsv')
        self.atg_fasta = os.path.join(data_dp, 'atg', 'atg.fasta')
        self.atg_display = os.path.join(data_dp, 'atg', 'atg_display.html')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')

    def readxl(self):
        sns = ['Cvt', 'Starvation-Induced', 'Core', 'Pexophagy', 'Subgroups']
        all_genes = []
        for sn in sns:
            df = pd.read_excel(self.atg_fp, sheetname=sn)
            all_genes += list(df['Gene'])
        all_genes = list(set(all_genes))
        all_genes = [gene.upper() for gene in all_genes]
        seqs, genes, orfs = self.get_yeast_seq_gene_from_ids(self.orf_trans, all_genes)
        print(set(all_genes) - set(genes))
        tools_fasta.yeast_write_fasta_from_ids(self.orf_trans, orfs, self.atg_fasta)
        self.display()


    def get_yeast_seq_gene_from_ids(self, orf_trans_fp, gene_ids):
        sequences = []
        genes = []
        orfs = []
        with open(orf_trans_fp, 'r') as fasta_in:
            for record in SeqIO.parse(fasta_in, 'fasta'):
                pid = str(record.id)
                full_description = str(record.description)
                fd_sp = full_description.split(',')
                gene = fd_sp[0].split(' ')[1]
                if gene in gene_ids:
                    sequences.append(str(record.seq))
                    genes.append(gene)
                    orfs.append(pid)
        return sequences, genes, orfs

    def display(self):
        ds = Display(self.atg_fasta, self.atg_display, color=False)
        ds.write_body()


def main():
    a = Atg()
    a.readxl()


if __name__ == '__main__':
    main()