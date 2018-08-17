import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta


class Labels(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.label_dp = os.path.join(data_dp, 'annotations')
        self.enz_labels = ['hydrolase', 'isomerase', 'ligase', 'lyase',
                    'oxidoreductase', 'transferase']
        self.labels = ['Pbody', 'RNA_binding', 'cytoplasmic_stress_granule']

    def get_enzyme_lists(self):
        enz_dict = {}
        apps = ['co', 'mc', 'ht']
        for fn in self.enz_labels:
            orfs = []
            for app in apps:
                ffn = "{}_activity_annotations_{}.txt".format(fn, app)
                ffp = os.path.join(self.label_dp, ffn)
                if os.path.exists(ffp):
                    df = pd.read_csv(ffp, sep='\t', comment='!')
                    orfs += list(df['Gene Systematic Name'])
            enz_dict[fn] = set(orfs)
        return enz_dict

    def get_labels(self):
        lab_dict = {}
        apps = ['co', 'mc', 'ht']
        for fn in self.labels:
            orfs = []
            for app in apps:
                ffn = "{}_annotations_{}.txt".format(fn, app)
                ffp = os.path.join(self.label_dp, ffn)
                if os.path.exists(ffp):
                    df = pd.read_csv(ffp, sep='\t', comment='!')
                    orfs += list(df['Gene Systematic Name'])
            lab_dict[fn] = set(orfs)
        return lab_dict


def main():
    en = Labels()
    print(en.get_labels())


if __name__ == '__main__':
    main()