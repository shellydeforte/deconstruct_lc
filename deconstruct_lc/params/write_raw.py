import pandas as pd
from deconstruct_lc.params import lc_labels
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc

class WriteRaw(object):
    def __init__(self, k, pdb_fp, bc_fp):
        self.pdb_fp = pdb_fp
        self.bc_fp = bc_fp
        self.all_ids, self.all_seqs, self.y = self.get_seqs()
        self.all_lens = tools_fasta.get_lengths(self.all_seqs)
        self.k = k

    def get_seqs(self):
        pdb_ids, pdb_seqs = tools_fasta.fasta_to_id_seq(self.pdb_fp)
        cb_ids, cb_seqs = tools_fasta.fasta_to_id_seq(self.bc_fp)
        all_seqs = cb_seqs + pdb_seqs
        all_ids = cb_ids + pdb_ids
        y = [0]*len(cb_ids) + [1]*len(pdb_ids)
        return all_ids, all_seqs, y


    def write_lca(self, alph):
        lc_labs = lc_labels.GetLabels(self.k)
        lcas = lc_labs.create_lcas(alph)
        df_dict = {}
        df_dict['Protein ID'] = self.all_ids
        df_dict['Length'] = self.all_lens
        df_dict['y'] = self.y
        for k_lca in lcas:
            lca = k_lca.split('_')[1]
            scores = tools_lc.calc_lca_motifs(self.all_seqs, self.k, lca)
            df_dict[k_lca] = scores
        df = pd.DataFrame(df_dict, columns=['Protein ID', 'Length', 'y']+lcas)
        return df

    def write_lce(self):
        lc_labs = lc_labels.GetLabels(self.k)
        lces = lc_labs.create_lces(self.all_seqs)
        df_dict = {}
        df_dict['Protein ID'] = self.all_ids
        df_dict['Length'] = self.all_lens
        df_dict['y'] = self.y
        for k_lce in lces:
            lce = float(k_lce.split('_')[1])
            scores = tools_lc.calc_lce_motifs(self.all_seqs, self.k, lce)
            df_dict[k_lce] = scores
        df = pd.DataFrame(df_dict, columns=['Protein ID', 'Length', 'y']+lces)
        return df


def main():
    pass


if __name__ == '__main__':
    main()