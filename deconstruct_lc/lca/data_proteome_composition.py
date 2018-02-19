from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import defaultdict
import os
import pandas as pd

from deconstruct_lc import read_config


class LcProteome(object):
    """
    For each proteome, record the amino acid composition for that proteome as
    a continuous sequence string.
    """
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.seg_dpi = os.path.join(data_dp, 'proteomes', 'euk_seg')
        self.fns = os.listdir(self.seg_dpi)
        self.fpo = os.path.join(data_dp, 'proteomes_analysis', 'lc_composition.tsv')

    def write_all_comps(self):
        aa_order = 'SGEQAPDTNKRLHVYFIMCW'
        all_perc = defaultdict(list)
        for fasta_in in self.fns:
            aa_dict = self._one_organism(fasta_in)
            for aa in aa_order:
                all_perc[aa].append(aa_dict[aa])
        cols = ['Filename'] + [aa for aa in aa_order]
        all_perc['Filename'] = self.fns
        df = pd.DataFrame(all_perc, columns=cols)
        df.to_csv(self.fpo, sep='\t')

    def _one_organism(self, fasta_in):
        all_aa = ''
        fasta_in = os.path.join(self.seg_dpi, fasta_in)
        with open(fasta_in, 'r') as fasta_in:
            for record in SeqIO.parse(fasta_in, 'fasta'):
                sequence = str(record.seq)
                for aa in sequence:
                    if aa.islower():
                        all_aa += aa
        analyzed_sequence = ProteinAnalysis(all_aa)
        return analyzed_sequence.get_amino_acids_percent()


def main():
    lcp = LcProteome()
    lcp.write_all_comps()


if __name__ == '__main__':
    main()