import os
import configparser
from deconstruct_lc.params import write_raw

config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))

bc_fp = config['filepaths']['cb_fp']
pdb_fp = config['filepaths']['pdb_train_fp']

def main():
    k = 4
    alph = 'SGE'
    wr = write_raw.WriteRaw(k, pdb_fp, bc_fp)
    df = wr.write_lce()
    print(df)


if __name__ == '__main__':
    main()