import re


def remove_histags(seqs):
    """Remove histags from a list of sequences and return a list of the new
    histags as well as the number of sequences with histags removed"""
    nseqs = []
    num_his = 0
    for seq in seqs:
        indexes = find_histag(seq)
        if len(indexes) > 0:
            num_his += 1
        nseq = slice_seq(indexes, seq)
        nseqs.append(nseq)
    return nseqs, num_his


def remove_histag_miss(seq, miss_seq):
    """Remove histags from the sequence and the corresponding missing
    sequence"""
    indexes = find_histag(seq)
    nseq = slice_seq(indexes, seq)
    nmseq = slice_seq(indexes, miss_seq)
    return nseq, nmseq


def find_histag(seq):
    regex = r'H{6}H*'
    match = re.finditer(regex, seq)
    indexes = []
    for item in match:
        indexes.append(item.start())
        indexes.append(item.end())
    return indexes


def slice_seq(indexes, seq):
    if len(indexes) > 0:
        nseq = seq[:indexes[0]]
        for i in range(1, len(indexes) - 1, 2):
            nseq += seq[indexes[i]:indexes[i + 1]]
        nseq += seq[indexes[-1]:]
    else:
        nseq = seq
    return nseq


def main():
    pass


if __name__ == '__main__':
    main()