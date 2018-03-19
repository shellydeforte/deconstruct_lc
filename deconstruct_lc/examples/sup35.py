import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.display import display_lc


class Sup35(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.sup_fp = os.path.join(data_dp, 'examples', 'S288C_YDR172W_SUP35_protein.fsa')

    def format_seq(self):
        seq = tools_fasta.fasta_to_seq(self.sup_fp)
        disp = display_lc.Display(seq, 'sup35.html')
        disp.write_body()


def main():
    sup = Sup35()
    sup.format_seq()


if __name__ == '__main__':
    main()