import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.display.display_lc import Display


class EisosomeScores(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.eis_ids = os.path.join(data_dp, 'puncta', 'eisosome_annotations.txt')
        self.fasta_fp = os.path.join(data_dp, 'puncta', 'eisosome_fasta.fsa')
        self.display_fp = os.path.join(data_dp, 'puncta', 'eisosome.html')

    def display(self):
        ds = Display(self.fasta_fp, self.display_fp, color=True)
        ds.write_body()

    def write_fasta(self):
        df = pd.read_csv(self.eis_ids, sep='\t', header=7)
        pids = set(list(df['Gene Systematic Name']))
        tools_fasta.yeast_write_fasta_from_ids(self.orf_trans, pids, self.fasta_fp)


def main():
    es = EisosomeScores()
    es.display()


if __name__ == '__main__':
    main()