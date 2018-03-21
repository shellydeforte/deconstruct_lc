from deconstruct_lc import read_config
from deconstruct_lc.data_pdb.ssdis_to_fasta import SsDis
from deconstruct_lc.data_pdb.filter_pdb import PdbFasta
from deconstruct_lc.data_pdb.norm_all_to_tsv import FastaTsv
from deconstruct_lc.data_pdb.write_pdb_analysis import PdbAnalysis


class RunPdbPrep(object):
    def __init__(self):
        self.config = read_config.read_config()

    def run_ssdis(self):
        """
        Convert ss_dis.txt to three fasta files: disorder, secondary structure,
        and sequence
        """
        ssd = SsDis(self.config)
        ssd.seq_dis_to_fasta()
        ssd.ss_to_fasta()
        ssd.verify_ss_dis_to_fasta()

    def run_filterpdb(self):
        """
        Filter pdb by x-ray only, eukaryote only, standard amino acid alphabet,
        and then create two files, one including sequences that have missing
        residues, and one that does not.
        """
        pdb = PdbFasta(self.config)
        pdb.create_pdb_miss()
        pdb.create_pdb_nomiss()

    def run_allnormtsv(self):
        """
        Write tsv files that include secondary structure for the norm and all
        datasets. The all dataset will be used to create the PDB analysis dataset.
        """
        ft = FastaTsv(self.config)
        ft.write_tsv()

    def run_pdbanalysis(self):
        pa = PdbAnalysis(self.config)
        pa.write_analysis()


def main():
    rr = RunPdbPrep()
    rr.run_ssdis()
    rr.run_filterpdb()
    rr.run_allnormtsv()
    rr.run_pdbanalysis()


if __name__ == '__main__':
    main()