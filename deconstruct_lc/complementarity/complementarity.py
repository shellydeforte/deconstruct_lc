"""
The goal here is to see if the blobs that I'm measuring might have some kind
of complementarity.

Ideally, it would be good to check this within a droplet, however, maybe a
good place to start is within the protein.

It wouldn't be too hard to check within a droplet. I could check within
stress granules for instance. I could see which proteins seem to interact
with each other and see if they have complementary blobs.

But maybe first see if I find *anything* within the proteins. I should find
something.

So what would complementarity look like? Well I have a list of complementary amino acids.

I am going to start by looking in the motifs.

Start with charge complementarity within a stretch
"""
import os
import pandas as pd
from deconstruct_lc.complementarity import motif_seq

class Complementarity(object):

    def __init__(self, nmo_fpi, lca_label, lce_label):
        self.nmo_fpi = nmo_fpi
        self.lca_label = lca_label
        self.lce_label = lce_label

    def check_comp(self):
        """What motifs occur together?"""
        nmo_df = pd.read_csv(self.nmo_fpi, sep='\t', index_col=0)
        nmo_df = nmo_df[(nmo_df[self.lca_label] > 100)]
        nmo_df = nmo_df.sort_values(by=[self.lca_label])
        nmo_df = nmo_df.reset_index(drop=True)
        for i, row in nmo_df.iterrows():
            sequence = row['Sequence']
            print(sequence)
            print(row['Protein ID'])
            ms = motif_seq.LcSeq(sequence, 6, 'SGEQAPDTNKR', 'lca')
            motifs = ms.list_motifs()
            print(len(motifs))
            alphs = self.get_motif_comp(motifs)
            print(alphs)

    def get_motif_comp(self, motifs):
        alphs = {'ST': 0, 'ED': 0, 'RK': 0, 'QN': 0, 'GA': 0, 'P': 0, 'ev': 0}
        for motif in motifs:
            flag = True
            ch = self.get_net_charge(motif)
            if ch <= -1:
                alphs['ED'] += 1
                flag = False
            if ch >= 1:
                alphs['RK'] += 1
                flag = False
            if motif.count('P') >= 3:
                alphs['P'] += 1
                flag = False
            if (motif.count('Q') + motif.count('N')) >= 3:
                alphs['QN'] += 1
                flag = False
            if (motif.count('S') + motif.count('T')) >= 3:
                alphs['ST'] += 1
                flag = False
            if (motif.count('G') + motif.count('A')) >= 4:
                alphs['GA'] += 1
                flag = False
            if flag:
                alphs['ev'] += 1
        return alphs


    def get_net_charge(self, motif):
        charges = {'E': -1, 'D': -1, 'R': 1, 'K': 1}
        charge_total = 0
        for aa in motif:
            if aa in charges:
                charge_total += charges[aa]
        return charge_total



class Pipeline(object):

    def __init__(self):
        self.base_fp = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.nmo_fpi = os.path.join(self.base_fp, 'scores',
                                    'nmo_6_SGEQAPDTNKR_6_1.6_seq_scores.tsv')
        self.lca_label = '6_SGEQAPDTNKR'
        self.lce_label = '6_1.6'

    def run_charge(self):
        comp = Complementarity(self.nmo_fpi, self.lca_label, self.lce_label)
        comp.check_comp()


def main():
    pipe = Pipeline()
    pipe.run_charge()


if __name__ == '__main__':
    main()