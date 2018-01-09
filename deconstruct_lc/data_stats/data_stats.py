import os
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from predict_llps import tools_fasta

class DataStats(object):
    def __init__(self):
        self.fd = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.cb90 = os.path.join(self.fd, 'quickgo',
                                           'quickgo_cb_cd90.fasta')
        self.pdb90 = os.path.join(self.fd, 'pdb_nomiss_cd90.fasta')

    def length_stats(self):
        cb_seqs = tools_fasta.fasta_to_seq(self.cb90, minlen=100,
                                             maxlen=2000, unique=True)
        pdb_seqs = tools_fasta.fasta_to_seq(self.pdb90, minlen=100,
                                             maxlen=2000, unique=True)
        cb_lens = tools_fasta.get_lengths(cb_seqs)
        pdb_lens = tools_fasta.get_lengths(pdb_seqs)
        cb_heights, cb_bins = np.histogram(cb_lens, bins=20)
        cbn_heights = cb_heights / sum(cb_heights)


        pdb_heights, pdb_bins = np.histogram(pdb_lens, bins=20)
        pdbn_heights = pdb_heights / sum(pdb_heights)
        plt.bar(pdb_bins[:-1], pdbn_heights, width=(max(pdb_bins) - min(
            pdb_bins)) / len(pdb_bins), color="darkblue", alpha=0.7,
                label='PDB')

        plt.bar(cb_bins[:-1], cbn_heights, width=(max(cb_bins) - min(
            cb_bins)) / len(cb_bins), color="orangered", alpha=0.7,
                label='Condensates')

        #plt.hist(pdb_lens, 10, normed=1, edgecolor='darkblue', alpha=0.5,
        #         cumulative=False, histtype='step', lw=3, label='PDB')
        #plt.hist(cb_lens, 10, normed=1, edgecolor='orangered', alpha=0.5,
        #         cumulative=False, histtype='step', lw=3, label='Condensates')
        #plt.yticks([], [])
        plt.xlabel('Protein Length')
        plt.ylabel('Relative Fraction')
        plt.legend()
        plt.show()

    def composition(self):
        cb_seqs = tools_fasta.fasta_to_seq(self.cb90, minlen=100,
                                             maxlen=2000, unique=True)
        pdb_seqs = tools_fasta.fasta_to_seq(self.pdb90, minlen=100,
                                             maxlen=2000, unique=True)
        aas = 'SGEQAPDTNKRLHVYFIMCW'
        aas_list = [aa for aa in aas]
        ind = range(len(aas))
        pdb_seq = ''
        for seq in pdb_seqs:
            pdb_seq += seq
        cb_seq = ''
        for seq in cb_seqs:
            cb_seq += seq
        an_pdb_seq = ProteinAnalysis(pdb_seq)
        pdb_dict = an_pdb_seq.get_amino_acids_percent()
        an_cb_seq = ProteinAnalysis(cb_seq)
        cb_dict = an_cb_seq.get_amino_acids_percent()
        pdb_bins = []
        cb_bins = []
        for aa in aas:
            pdb_bins.append(pdb_dict[aa])
            cb_bins.append(cb_dict[aa])
        plt.bar(ind, pdb_bins, color='darkblue', alpha=0.7, label='PDB')
        plt.bar(ind, cb_bins, color='orangered', alpha=0.7,
                label='Condensates')
        plt.xticks(ind, aas_list)
        plt.legend()
        plt.xlabel('Amino Acids')
        plt.ylabel('Relative Fraction in full dataset')
        plt.show()

    def overlap(self):
        pdb_uni = self.read_pdb_uni()
        cb_ids, cb_seqs = tools_fasta.fasta_to_id_seq(self.cb90)
        cb_pdb_unis = {}
        cb_pdbs = []
        for id in cb_ids:
            if id in pdb_uni:
                cb_pdb_unis[pdb_uni[id]] = id
                cb_pdbs.append(pdb_uni[id])
        print(cb_pdb_unis)
        pdb_ids, pdb_seqs = tools_fasta.fasta_to_id_seq(self.pdb90)
        for pdb_id in pdb_ids:
            if pdb_id in cb_pdbs:
                print(pdb_id)
                print(cb_pdb_unis[pdb_id])



    def read_pdb_uni(self):
        pdb_uni = {}
        pdb_uni_fp = os.path.join(self.fd, 'pdb_chain_uniprot.tsv')
        with open(pdb_uni_fp, 'r') as fi:
            for i in range(0, 2):
                next(fi)
            for line in fi:
                ls = line.split('\t')
                pdb = ls[0].upper()
                chain = ls[1].upper()
                uni = ls[2]
                pdb_chain = '{}_{}'.format(pdb, chain)
                pdb_uni[uni] = pdb_chain
        return pdb_uni


def main():
    ds = DataStats()
    ds.overlap()


if __name__ == '__main__':
    main()


