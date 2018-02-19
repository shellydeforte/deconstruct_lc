from deconstruct_lc import read_config
from deconstruct_lc.params import raw_scores
from deconstruct_lc.params import raw_svm
from deconstruct_lc.params import raw_mb


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

    def run_rawmb(self):
        rm = raw_mb.RawMb(self.config)
        rm.write_top()


def main():
    rr = RunRaw()
    rr.run_rawmb()


if __name__ == '__main__':
    main()