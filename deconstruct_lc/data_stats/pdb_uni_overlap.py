import os

from deconstruct_lc import tools_fasta
from deconstruct_lc import read_config


class Overlap(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.bc90 = os.path.join(data_dp, 'bc_prep', 'bc_train_cd90.fasta')
        self.pdb90 = os.path.join(data_dp, 'pdb_prep', 'pdb_train_cd90.fasta')
        self.pdb_chain = os.path.join(data_dp, 'pdb_prep', 'pdb_chain_uniprot.tsv')

    def overlap(self):
        pdb_uni = self.read_pdb_uni()
        cb_ids, cb_seqs = tools_fasta.fasta_to_id_seq(self.bc90)
        cb_pdb_unis = {}
        cb_pdbs = []
        for id in cb_ids:
            if id in pdb_uni:
                cb_pdb_unis[pdb_uni[id]] = id
                cb_pdbs.append(pdb_uni[id])
        pdb_ids, pdb_seqs = tools_fasta.fasta_to_id_seq(self.pdb90)
        print("Proteins overlapping between the PDB and BC datasets")
        for pdb_id in pdb_ids:
            if pdb_id in cb_pdbs:
                print(pdb_id)
                print(cb_pdb_unis[pdb_id])

    def read_pdb_uni(self):
        pdb_uni = {}
        with open(self.pdb_chain, 'r') as fi:
            for i in range(0, 2):
                next(fi)
            for line in fi:
                ls = line.split('\t')
                pdb = ls[0].upper()
                chain = ls[1].upper()
                uni = ls[2]
                pdb_chain = '{}_{}'.format(pdb, chain)
                pdb_uni[uni] = pdb_chain
        return pdb_uni


def main():
    ol = Overlap()
    ol.overlap()


if __name__ == '__main__':
    main()