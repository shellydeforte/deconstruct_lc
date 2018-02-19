from deconstruct_lc import read_config
from deconstruct_lc.params import raw_scores
from deconstruct_lc.params import raw_svm


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


def main():
    rr = RunRaw()


if __name__ == '__main__':
    main()