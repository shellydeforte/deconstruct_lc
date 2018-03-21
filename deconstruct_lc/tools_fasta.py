import os
import re
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


def fasta_to_head_seq(fasta_fp, minlen=0, maxlen=float('inf'), unique=False):
    sequences = []
    headers = []
    with open(fasta_fp, 'r') as file_in:
        for record in SeqIO.parse(file_in, 'fasta'):
            sequence = str(record.seq)
            if minlen <= len(sequence) <= maxlen:
                desc = str(record.description)
                if unique:
                    if sequence not in sequences:
                        headers.append(desc)
                        sequences.append(sequence)
                else:
                    headers.append(desc)
                    sequences.append(sequence)
    return headers, sequences


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


def get_yeast_seq_from_ids(orf_trans_fp, orf_ids):
    sequences = []
    with open(orf_trans_fp, 'r') as fasta_in:
        for record in SeqIO.parse(fasta_in, 'fasta'):
            pid = str(record.id)
            if pid in orf_ids:
                sequences.append(str(record.seq))
    return sequences


def get_yeast_seq_gene_from_ids(orf_trans_fp, orf_ids):
    sequences = []
    genes = []
    with open(orf_trans_fp, 'r') as fasta_in:
        for record in SeqIO.parse(fasta_in, 'fasta'):
            pid = str(record.id)
            if pid in orf_ids:
                full_description = str(record.description)
                fd_sp = full_description.split(',')
                gene = fd_sp[0].split(' ')[1]
                sequences.append(str(record.seq))
                genes.append(gene)
    return sequences, genes


def yeast_write_fasta_from_ids(orf_trans_fp, orf_ids, fasta_out):
    records = []
    with open(orf_trans_fp, 'r') as fasta_in:
        for record in SeqIO.parse(fasta_in, 'fasta'):
            pid = str(record.id)
            if pid in orf_ids:
                records.append(record)
    SeqIO.write(records, fasta_out, "fasta")


def get_lengths(seqs):
    lengths = [len(seq) for seq in seqs]
    return lengths


def remove_all_histags(seqs):
    nseqs = []
    for seq in seqs:
        nseqs.append(remove_histag(seq))
    return nseqs


def remove_histag(seq):
    """If H*6 or greater, remove from sequence"""
    regex = r'H{6}H*'
    #nseq = re.sub(regex, '', seq)
    match = re.finditer(regex, seq)
    indexes = []
    for item in match:
        indexes.append(item.start())
        indexes.append(item.end())
    nseq = seq[:indexes[0]]
    for i in range(1,len(indexes)-1, 2):
        nseq += seq[indexes[i]:indexes[i+1]]
    nseq += seq[indexes[-1]:]
    return nseq


def main():
    seq = 'ABCHHHHHHHREBHHHHHHQ'
    nseq = remove_histag(seq)
    print(nseq)


if __name__ == '__main__':
    main()


