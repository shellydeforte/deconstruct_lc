from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import defaultdict
import numpy as np
import os
import pandas as pd

from deconstruct_lc import read_config


class LcProteome(object):
    """
    Read all seg proteome files and output a tsv file with the mean
    composition and standard deviation. The composition is calculated for
    the entire organism and the mean and standard deviation are taken from
    the list of these values.
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
        means = []
        stds = []
        aas = []
        for aa in aa_order:
            means.append(np.mean(all_perc[aa]))
            stds.append(np.std(all_perc[aa]))
            aas.append(aa)
        df_dict = {'Amino Acid': aas, 'Mean': means, 'Standard Deviation':
            stds}
        df = pd.DataFrame(df_dict, columns=['Amino Acid',
                                            'Mean', 'Standard Deviation'])
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

    def write_n(self):
        all_perc = []
        for fasta_in in self.fns:
            aa_dict = self._one_organism(fasta_in)
            all_perc.append(aa_dict['N'])
        df_dict = {'Filename': self.fns, 'Fraction of Asparagine': all_perc}
        df = pd.DataFrame(df_dict, columns=['Filename', 'Fraction of '
                                                        'Asparagine'])
        fpo = os.path.join(self.fpo, 'asparagine.tsv')
        df.to_csv(fpo, sep='\t')


def main():
    lcp = LcProteome()
    lcp.write_all_comps()


if __name__ == '__main__':
    main()