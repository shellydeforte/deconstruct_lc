import os
from Bio import SeqIO


def fasta_to_seq(fasta_fp, minlen=0, maxlen=float('inf'), unique=False):
    """
    Return a list of sequences from the given fasta file.

    if unique = True, this function will only include the first unique
    sequence in the file.
    """
    sequences = []
    with open(fasta_fp, 'r') as file_in:
        for record in SeqIO.parse(file_in, 'fasta'):
            sequence = str(record.seq)
            if minlen <= len(sequence) <= maxlen:
                if unique:
                    if sequence not in sequences:
                        sequences.append(sequence)
                else:
                    sequences.append(sequence)
    return sequences


def fasta_to_id_seq(fasta_fp, minlen=0, maxlen=float('inf'), unique=False):
    """
    Return a list of ids and a list of the corresponding sequences from a
    fasta file.

    if unique = True, this function will only include the first unique
    sequence and id in the file.
    """
    sequences = []
    pids = []
    with open(fasta_fp, 'r') as file_in:
        for record in SeqIO.parse(file_in, 'fasta'):
            sequence = str(record.seq)
            if minlen <= len(sequence) <= maxlen:
                pid = id_cleanup(str(record.id))
                if unique:
                    if sequence not in sequences:
                        pids.append(pid)
                        sequences.append(sequence)
                else:
                    pids.append(pid)
                    sequences.append(sequence)
    return pids, sequences


def id_cleanup(protein_id):
    if '|' in protein_id:
        nid = protein_id.split('|')[1]
    elif ':' in protein_id:
        ps = protein_id.split(':')
        nid = '{}_{}'.format(ps[0], ps[1])
    else:
        nid = protein_id
    return nid


def get_pid_gene_desc_seq(fasta_fp):
    """This is specific to yeast fasta files"""
    pids = []
    genes = []
    seqs = []
    descs = []
    with open(fasta_fp, 'r') as fasta_in:
        for record in SeqIO.parse(fasta_in, 'fasta'):
            full_description = str(record.description)
            fd_sp = full_description.split(',')
            pid = str(record.id)
            gene = fd_sp[0].split(' ')[1]
            seq = str(record.seq)
            fd_sp_q = full_description.split('"')
            desc = fd_sp_q[1]
            pids.append(pid)
            genes.append(gene)
            seqs.append(seq)
            descs.append(desc)
    return pids, genes, seqs, descs


def get_lengths(seqs):
    lengths = [len(seq) for seq in seqs]
    return lengths


def main():
    pass


if __name__ == '__main__':
    main()


