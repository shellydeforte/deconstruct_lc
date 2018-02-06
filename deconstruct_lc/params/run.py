from deconstruct_lc import read_config
from deconstruct_lc.params import raw_scores


class RunRaw(object):
    def __init__(self):
        self.config = read_config.read_test_config()

    def run(self):
        pr = raw_scores.PipeRaw(self.config)
        pr.write_lca()
        pr.write_lce()


def main():
    rr = RunRaw()
    rr.run()


if __name__ == '__main__':
    main()