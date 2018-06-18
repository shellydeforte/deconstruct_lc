import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.display import display_lc
from deconstruct_lc.scores.norm_score import NormScore


class Sup35(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.sup_fp = os.path.join(data_dp, 'examples', 'S288C_YDR172W_SUP35_protein.fsa')

    def format_seq(self):
        seq = tools_fasta.fasta_to_seq(self.sup_fp)
        disp = display_lc.Display(seq, 'sup35.html')
        disp.write_body()

    def score_by_section(self):
        seq = tools_fasta.fasta_to_seq(self.sup_fp)[0][0:-1]
        ns = NormScore()
        nm_domain = seq[0:253]
        c_domain = seq[253:]
        mc_domain = seq[123:]
        print(nm_domain)
        print(c_domain)
        print(mc_domain)
        print(ns.lc_norm_score([nm_domain]))
        print(ns.lc_norm_score([c_domain]))
        print(ns.lc_norm_score([mc_domain]))


def main():
    sup = Sup35()
    sup.score_by_section()


if __name__ == '__main__':
    main()