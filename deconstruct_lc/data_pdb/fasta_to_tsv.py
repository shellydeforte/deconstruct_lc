from Bio import SeqIO
import os

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta


class FastaTsv(object):
    """
    Fasta files are converted to tsv files with His-Tags removed
    """
    def __init__(self):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.pdb_dp = os.path.join(data_dp, 'pdb_prep')
        self.norm_fpi = os.path.join(self.pdb_dp, 'pdb_norm_cd100.fasta')
        self.all_fpi = os.path.join(self.pdb_dp, 'pdb_all.fasta')
        self.norm_fpo = os.path.join(self.pdb_dp, 'pdb_norm_cd100.tsv')
        self.all_fpo = os.path.join(self.pdb_dp, 'pdb_all.tsv')
        self.all_seq = os.path.join(self.pdb_dp, 'all_seqs.fasta')
        self.all_dis = os.path.join(self.pdb_dp, 'all_dis.fasta')
        self.all_ss = os.path.join(self.pdb_dp, 'all_ss.fasta')

    def write_tsv(self):
        self.write_full(self.all_fpi, self.all_fpo)
        self.write_full(self.norm_fpi, self.norm_fpo)

    def write_full(self, fasta, fpo):
        """
        Write sequence, missing, secondary structure if in the list of pids.
        """
        all_pids = self.get_pids(fasta)
        count = 0
        with open(self.all_seq, 'r') as seq_fi, \
             open(self.all_dis, 'r') as dis_fi, \
             open(self.all_ss, 'r') as ss_fi:
            with open(fpo, 'w') as fo:
                fo.write('Protein ID\tSequence\tMissing\tSecondary '
                         'Structure\n')
                for seq_rec, dis_rec, ss_rec in zip(SeqIO.parse(seq_fi, 'fasta'),
                                                    SeqIO.parse(dis_fi, 'fasta'),
                                                    SeqIO.parse(ss_fi, 'fasta')):
                    pid = tools_fasta.id_cleanup(seq_rec.id)
                    if pid in all_pids:
                        count += 1
                        print(count)
                        seq = str(seq_rec.seq)
                        mseq = str(dis_rec.seq)
                        ss_seq = str(ss_rec.seq)
                        assert len(seq) == len(mseq) == len(ss_seq)
                        fo.write('{}\t{}\t{}\t{}\n'.format(pid, seq, mseq,
                                                           ss_seq))

    def get_pids(self, fasta):
        pids, seqs = tools_fasta.fasta_to_id_seq(fasta)
        return pids


def main():
    ft = FastaTsv()
    ft.write_tsv()


if __name__ == '__main__':
    main()