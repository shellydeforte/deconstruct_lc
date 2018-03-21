from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

from deconstruct_lc import tools_fasta


class SsDis(object):
    def __init__(self, config):
        data_dp = config['fps']['data_dp']
        pdb_dp = os.path.join(data_dp, 'data_pdb')
        self.ss_dis_fp = os.path.join(pdb_dp, 'outside_data', 'ss_dis.txt')
        self.all_dis_fp = os.path.join(pdb_dp, 'all_dis.fasta')
        self.all_seq_fp = os.path.join(pdb_dp, 'all_seqs.fasta')
        self.all_ss_fp = os.path.join(pdb_dp, 'all_ss.fasta')

    def seq_dis_to_fasta(self):
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
        print("Done writing disorder and sequence files")

    def ss_to_fasta(self):
        """
        Read ss_dis.text. For the secondary structure file, add 'P' where
        there is a blank
        """
        ss_fp = self.ss_dis_fp
        ss_fpo = self.all_ss_fp
        new_fasta = []
        with open(ss_fp, 'r') as ss_fi:
            for line in ss_fi:
                if 'secstr' in line:
                    nid = line[1:].strip()
                    line = next(ss_fi)
                    nseq = ''
                    while line[0] != '>':
                        nseq += line[:-1]
                        line = next(ss_fi)
                    new_seq = self._add_p(nseq)
                    new_record = SeqRecord(Seq(new_seq, IUPAC.protein), id=nid,
                                           description='')
                    new_fasta.append(new_record)
        with open(ss_fpo, 'w') as output_handle:
            SeqIO.write(new_fasta, output_handle, 'fasta')
        print("Done writing secondary structure")

    def _add_p(self, sequence):
        new_seq = ''
        for aa in sequence:
            if aa == ' ':
                new_seq += 'P'
            else:
                new_seq += aa
        return new_seq

    def verify_ss_dis_to_fasta(self):
        """
        Confirm that protein IDs and sequence lengths are the same
        """
        total_entries = 0
        with open(self.all_seq_fp, 'r') as seq_fasta:
            with open(self.all_dis_fp, 'r') as dis_fasta:
                with open(self.all_ss_fp, 'r') as ss_fasta:
                    for seq_rec, dis_rec, ss_rec in \
                            zip(SeqIO.parse(seq_fasta, 'fasta'),
                                SeqIO.parse(dis_fasta, 'fasta'),
                                SeqIO.parse(ss_fasta, 'fasta')):
                        seq_id = tools_fasta.id_cleanup(seq_rec.id)
                        dis_id = tools_fasta.id_cleanup(dis_rec.id)
                        ss_id = tools_fasta.id_cleanup(ss_rec.id)
                        assert seq_id == dis_id == ss_id
                        assert len(seq_rec.seq) == len(dis_rec.seq) == len(
                            ss_rec.seq)
                        total_entries += 1
        print("ss_dis fasta files verified.")
        print("There are {} total entries from ss_dis.txt".format(total_entries))