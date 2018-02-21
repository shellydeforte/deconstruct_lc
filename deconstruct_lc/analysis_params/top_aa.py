import os
import pandas as pd

from deconstruct_lc import read_config

class TopAa(object):
    """Find the most important amino acids in the LCA"""
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        params_dp = os.path.join(data_dp, 'params')
        self.fpi = os.path.join(params_dp, 'top_svm_lca.tsv')

    def top_aa(self):
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        labels = df['Label']
        aas = {}
        for label in labels:
            lab = label.split('_')[1]
            for aa in lab:
                if aa in aas:
                    aas[aa] += 1
                else:
                    aas[aa] = 1
        naas = {}
        for aa in aas:
            naas[aa] = aas[aa]/300
        self.dict_to_df(naas)

    def dict_to_df(self, adict):
        vals = []
        aas = []
        for item in adict:
            vals.append(adict[item])
            aas.append(item)
        df = pd.DataFrame({'AA': aas, 'Fraction': vals})
        result = df.sort(['Fraction'], ascending=False)
        print(result)


def main():
    taa = TopAa()
    taa.top_aa()


if __name__ == '__main__':
    main()