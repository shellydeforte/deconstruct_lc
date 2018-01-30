import configparser
import os
from deconstruct_lc.params import lc_labels

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

class PrintStats(object):
    def __init__(self):
        self.k = 6
        self.alph = 'SGEQAPDTNKRL'

    def iterations(self):
        lc_labs = lc_labels.GetLabels(self.k)
        k_lcas = lc_labs.create_lcas(self.alph)
        print(len(k_lcas))


def main():
    ps = PrintStats()
    ps.iterations()


if __name__ == '__main__':
    main()