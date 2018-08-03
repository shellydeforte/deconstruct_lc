import os
import pandas as pd

from deconstruct_lc import read_config


class BcProteome(object):
    """Create contingency tables for BC data vs. proteome, with BC prtoteins
    removed from the proteome"""
    def __init__(self, bc, organism):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        scores_dp = os.path.join(data_dp, 'scores')
        self.bc = bc
        self.organism = organism
        self.fpi = os.path.join(scores_dp, 'pdb_bc_scores.tsv')

    def get_cont_table(self):
        """Organism can be 'Human', 'Yeast'"""
        bc_df = self.get_bc()
        prot_df = self.get_proteome()
        bc_counts = get_bins(bc_df)
        prot_counts = get_bins(prot_df)
        return bc_counts, prot_counts

    def get_bc_ids(self):
        """
        Given the proteome label and the BC label, return the BC Uniprot IDs
        """
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        bc_ids = list(set(df[(df['Proteome'] == self.bc) & (df['Organism'] == self.organism)]['Protein ID']))
        return bc_ids

    def get_proteome(self):
        """
        Given the proteome label and the BC label, return the proteome
        scores minus the values that were in the BC.
        """
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        bc_ids = self.get_bc_ids()
        df = df[(df['Organism'] == self.organism) & (df['Proteome'] == self.organism)]
        df = df[~df['Protein ID'].isin(bc_ids)]
        return df

    def get_bc(self):
        df = pd.read_csv(self.fpi, sep='\t', index_col=0)
        df = df[(df['Organism'] == self.organism.upper()) & (df['Proteome'] == self.bc)]
        return df


def get_bins(df):
    ndf = df[df['LC Score'] < 0]
    lt = len(ndf)
    ndf = df[(df['LC Score'] >= 0) & (df['LC Score'] <= 20)]
    m = len(ndf)
    ndf = df[df['LC Score'] > 20]
    gt = len(ndf)
    counts = (lt, m, gt)
    return counts


def main():
    bc = BcProteome('P_Body', 'Human')
    bc_counts, prot_counts = bc.get_cont_table()
    print(bc_counts)
    print(prot_counts)



if __name__ == '__main__':
    main()