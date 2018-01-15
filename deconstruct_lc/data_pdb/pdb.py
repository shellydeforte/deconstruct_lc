from Bio import SeqIO
import configparser
import os
from deconstruct_lc import tools_fasta


config = configparser.ConfigParser()
cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                   'config.cfg'))
config.read_file(open(cfg_fp, 'r'))


class PdbFasta(object):
    def __init__(self):
        self.minlen = 100
        self.maxlen = 2000
        self.pdb_dp = os.path.join(config['filepaths']['data_dp'], 'pdb_prep')
        self.pdb_miss_fp = os.path.join(self.pdb_dp, 'pdb_norm.fasta')
        self.pdb_nomiss_fp = os.path.join(self.pdb_dp, 'pdb_train.fasta')
        self.ss_dis_fp = os.path.join(self.pdb_dp, 'ss_dis.txt')
        self.all_dis_fp = os.path.join(self.pdb_dp, 'all_dis.fasta')
        self.all_seq_fp = os.path.join(self.pdb_dp, 'all_seqs.fasta')
        self.entry_type_fp = os.path.join(self.pdb_dp, 'pdb_entry_type.txt')
        self.taxonomy_fp = os.path.join(self.pdb_dp, 'pdb_chain_taxonomy.tsv')
        self.speclist_fp = os.path.join(self.pdb_dp, 'speclist.txt')

    def create_pdb_miss(self):
        """
        Apply the following filtering:
        x-ray
        eukaryote
        standard amino acid alphabet
        do not apply any missing region or length filtering
        """
        diffraction = self.get_diffraction()
        eukaryote = self.get_euk_pdb()
        new_records = []
        with open(self.all_seq_fp, 'r') as seq_fi:
            for seq_rec in SeqIO.parse(seq_fi, 'fasta'):
                pdb_chain = tools_fasta.id_cleanup(str(seq_rec.id))
                pdb = pdb_chain.split('_')[0]
                sequence = str(seq_rec.seq)
                if pdb_chain in eukaryote:
                    if pdb in diffraction:
                        if self.standard_aa(sequence):
                            new_records.append(seq_rec)
        with open(self.pdb_miss_fp, 'w') as seq_fo:
            SeqIO.write(new_records, seq_fo, 'fasta')
        count = 0
        with open(self.pdb_miss_fp, 'r') as handle:
            for _ in SeqIO.parse(handle, 'fasta'):
                count += 1
        print('There are {} records with missing regions'.format(count))

    def create_pdb_nomiss(self):
        """
        load pdb_miss and apply no missing regions filter, and length filter
        """
        pdb_miss = self.get_missing()
        new_records = []
        with open(self.pdb_miss_fp, 'r') as miss_fi:
            for seq_rec in SeqIO.parse(miss_fi, 'fasta'):
                pdb_chain = tools_fasta.id_cleanup(str(seq_rec.id))
                sequence = str(seq_rec.seq)
                prot_len = len(sequence)
                if pdb_chain in pdb_miss:
                    if pdb_miss[pdb_chain] == 0:
                        if self.minlen <= prot_len <= self.maxlen:
                            new_records.append(seq_rec)
        with open(self.pdb_nomiss_fp, 'w') as seq_fo:
            SeqIO.write(new_records, seq_fo, 'fasta')
        count = 0
        with open(self.pdb_nomiss_fp, 'r') as handle:
            for _ in SeqIO.parse(handle, 'fasta'):
                count += 1
        print('There are {} records without missing regions'.format(count))

    def get_euk_pdb(self):
        """
        Create a list of PDB IDs that are eukaryotes, ie.
        {'3V6M_A', '3WGU_E', '4ITZ_B',...
        """
        euks = self.get_euk_tax()
        euk_pdbs = []
        with open(self.taxonomy_fp, 'r') as fi:
            next(fi)
            next(fi)
            for line in fi:
                ls = line.split()
                if ls[2] in euks:
                    pdb_chain = '{}_{}'.format(ls[0].upper(), ls[1].upper())
                    euk_pdbs.append(pdb_chain)
        return set(euk_pdbs)

    def get_euk_tax(self):
        """
        Create a set of all taxonomic identifiers that are 'E' for eukaryote
        ie. {'348046', '160085', '143180',...
        """
        tax_org = self.read_speclist()
        euks = []
        for item in tax_org:
            if tax_org[item] == 'E':
                euks.append(item)
        return set(euks)

    def read_speclist(self):
        """
        Read spec list and create a tax: org dictionary, ie. '3320': 'E'
        """
        tax_org = {}
        with open(self.speclist_fp, 'r') as fi:
            for i in range(59):
                next(fi)
            for line in fi:
                ls = line.split()
                if len(ls) >= 3:
                    if ls[2][-1] == ':':
                        tax = ls[2][0:-1]
                        org = ls[1]
                        tax_org[tax] = org
        return tax_org

    def get_diffraction(self):
        """
        Read diffraction file and return a set of PDB IDs that are
        from crystal structures. Note this does not include chain info. ie.
        {'3S5T', '1F8W', '1PCX',...
        """
        xray_ids = []
        with open(self.entry_type_fp, 'r') as fi:
            for line in fi:
                ls = line.split()
                pdb_id = ls[0]
                method = ls[2]
                if method == 'diffraction':
                    xray_ids.append(pdb_id.upper())
        return set(xray_ids)

    def get_missing(self):
        """
        Create a dictionary with pdb_chain: # missing regions
        """
        pdb_miss = {}
        with open(self.all_dis_fp, 'r') as dis_fo:
            for dis_rec in SeqIO.parse(dis_fo, 'fasta'):
                pdb_miss[tools_fasta.id_cleanup(str(dis_rec.id))] = \
                    self.missing_count(str(dis_rec.seq))
        return pdb_miss

    def missing_count(self, dis_seq):
        return dis_seq.count('X')

    def standard_aa(self, sequence):
        aas = 'ADKERNTSQYFLIVMCWHGP'
        for c in sequence:
            if c not in aas:
                return False
        return True

    def ss_dis_to_fasta(self):
        """
        Read ss_dis.txt and create fasta files for sequence and disorder.
        """
        sequence = []
        disorder = []
        with open(self.ss_dis_fp, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                rid = str(record.id)
                if 'disorder' in rid:
                    disorder.append(record)
                elif 'sequence' in rid:
                    sequence.append(record)
                else:
                    pass
        with open(self.all_seq_fp, 'w') as output_sequence:
            SeqIO.write(sequence, output_sequence, 'fasta')
        with open(self.all_dis_fp, 'w') as output_disorder:
            SeqIO.write(disorder, output_disorder, 'fasta')

    def verify_ss_dis_to_fasta(self):
        """
        Confirm that protein IDs and sequence lengths are the same
        """
        with open(self.all_seq_fp, 'r') as seq_fasta:
            with open(self.all_dis_fp, 'r') as dis_fasta:
                for seq_rec, dis_rec in zip(SeqIO.parse(seq_fasta, 'fasta'),
                                            SeqIO.parse(dis_fasta, 'fasta')):
                    seq_id = tools_fasta.id_cleanup(seq_rec.id)
                    dis_id = tools_fasta.id_cleanup(dis_rec.id)
                    assert seq_id == dis_id
                    assert len(seq_rec.seq) == len(dis_rec.seq)


def main():
    pdb = PdbFasta()
    pdb.ss_dis_to_fasta()
    pdb.create_pdb_miss()
    pdb.create_pdb_nomiss()


if __name__ == '__main__':
    main()