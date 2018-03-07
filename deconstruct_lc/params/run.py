from deconstruct_lc import read_config
from deconstruct_lc.params import raw_scores
from deconstruct_lc.params import raw_svm
from deconstruct_lc.params import raw_top
from deconstruct_lc.params import write_mb
from deconstruct_lc.params import raw_norm


class RunRaw(object):
    def __init__(self):
        self.config = read_config.read_config()

    def run_raw_scores(self):
        pr = raw_scores.RawScores(self.config)
        pr.write_lca()
        pr.write_lce()

    def run_svm(self):
        rs = raw_svm.RawSvm(self.config)
        rs.svm_lca_lce()

    def run_rawtop(self):
        rm = raw_top.RawTop(self.config)
        rm.write_top()

    def run_writemb(self):
        mb = write_mb.WriteMb(self.config)
        #mb.write_mb_solo()
        mb.write_mb_combos()

    def run_rawnorm(self):
        rn = raw_norm.RawNorm(self.config)


def main():
    rr = RunRaw()
    rr.run_writemb()


if __name__ == '__main__':
    main()