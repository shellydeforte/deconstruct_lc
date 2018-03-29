import os

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores import norm_score

class Sandbox(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        exp_dp = os.path.join(data_dp, 'experiment')
        self.lsm4_fp = os.path.join(exp_dp, 'S288C_YER112W_LSM4_protein.fsa')
        self.dhh1_fp = os.path.join(exp_dp, 'S288C_YDL160C_DHH1_protein.fsa')
        self.edc3_fp = os.path.join(data_dp, 'S288C_YEL015W_EDC3_protein.fsa')
        self.pat1_fp = os.path.join(data_dp, 'S288C_YCR077C_PAT1_protein.fsa')

    def scores(self):
        """
        edc3 and pat1 are currently tagged with fluoresence
        """
        #lsm4_seq = tools_fasta.fasta_to_seq(self.lsm4_fp)
        #dhh1_seq = tools_fasta.fasta_to_seq(self.dhh1_fp)
        edc3_seq = tools_fasta.fasta_to_seq(self.edc3_fp)
        pat1_seq = tools_fasta.fasta_to_seq(self.pat1_fp)
        print(edc3_seq)
        print(pat1_seq)
        ns = norm_score.NormScore()
        #print(ns.lc_norm_score(lsm4_seq))
        #print(ns.lc_norm_score(dhh1_seq))
        print(ns.lc_norm_score(edc3_seq))
        print(ns.lc_norm_score(pat1_seq))
        #lsm4_del = lsm4_seq[0][0:91] + lsm4_seq[0][187:]
        #dhh1_del = dhh1_seq[0][0:434] + dhh1_seq[0][506:]
        edc3_del = edc3_seq[0][0:65] + edc3_seq[0][186:]
        pat1_del = pat1_seq[0][0:133] + pat1_seq[0][200:]
        print()
        #print(ns.lc_norm_score([lsm4_del]))
        #print(ns.lc_norm_score([dhh1_del]))
        print(ns.lc_norm_score([edc3_del]))
        print(ns.lc_norm_score([pat1_del]))


def main():
    sb = Sandbox()
    sb.scores()


if __name__ == '__main__':
    main()