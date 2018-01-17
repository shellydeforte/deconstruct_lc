"""
Created by Shelly DeForte, Michnick Lab, University of Montreal 2018
"""
import math
from deconstruct_lc import tools


def seq_to_kmers(sequence, k):
    """Given a sequence, return a list of all overlapping k-mers"""
    i = 0
    len_sequence = len(sequence)
    kmers = []
    while i+k <= len_sequence:
        kmers.append(sequence[i:i+k])
        i += 1
    return kmers


def calc_lce_motifs(sequences, k, lce):
    """Given a list of sequences, return a list of motif counts for each
    sequence"""
    motif_counts = []
    for sequence in sequences:
        motif_count = count_lce_motifs(sequence, k, lce)
        motif_counts.append(motif_count)
    return motif_counts


def calc_lca_motifs(sequences, k, lca):
    """Given a list of sequences, return a list of motif counts for each
    sequence"""
    motif_counts = []
    for sequence in sequences:
        motif_count = count_lca_motifs(sequence, k, lca)
        motif_counts.append(motif_count)
    return motif_counts


def calc_lc_motifs(sequences, k, lca, lce):
    """Calculate the total number of unique lca or lce motifs
    k must be the same"""
    motif_counts = []
    for sequence in sequences:
        motif_count = count_lc_motifs(sequence, k, lca, lce)
        motif_counts.append(motif_count)
    return motif_counts


def count_lc_motifs(sequence, k, lca, lce):
    kmers = seq_to_kmers(sequence, k)
    motif_count = 0
    for kmer in kmers:
        if lca_motif(kmer, lca):
            motif_count += 1
        elif lce_motif(kmer, lce):
            motif_count += 1
        else:
            pass
    return motif_count


def count_lca_motifs(sequence, k, lca):
    """Count the number of lca motifs of length k in a sequence"""
    kmers = seq_to_kmers(sequence, k)
    motif_count = 0
    for kmer in kmers:
        if lca_motif(kmer, lca):
            motif_count += 1
    return motif_count


def count_lce_motifs(sequence, k, threshold):
    """Count the number of lce motifs of length k in a sequence"""
    kmers = seq_to_kmers(sequence, k)
    motif_count = 0
    for kmer in kmers:
        if lce_motif(kmer, threshold):
            motif_count += 1
    return motif_count


def lca_motif(kmer, lca):
    """Checks to see if the sequence contains only those amino acids as
    defined in the lca (a string)'"""
    in_motif = set(lca)
    if set(kmer) <= in_motif:
        return True
    else:
        return False


def lce_motif(kmer, threshold):
    """Checks if the sequence is les than or equal to the threshold in its
    shannon entropy score"""
    h = shannon(kmer)
    if h <= threshold:
        return True
    else:
        return False


def shannon(astring):
    """Calculates shannon entropy for any string with log base 2"""
    entropy = 0
    len_str = float(len(astring))
    unique = set(astring)
    for c in unique:
        p_x = float(astring.count(c))/len_str
        entropy += p_x*math.log(p_x, 2)
    if not entropy == 0.0:
        entropy = -(entropy)
    return entropy


def lca_to_interval(sequence, k, lca):
    """
    Returns inclusive interval, where all numbers are in the motif, ie (0, 6)
    and not (0, 7)
    """
    kmers = seq_to_kmers(sequence, k)
    indexes = set()
    for i, kmer in enumerate(kmers):
        if lca_motif(kmer, lca):
            for j in range(i, i+k):
                indexes.add(j)
    intervals = tools.ints_to_ranges(sorted(list(indexes)))
    return intervals


def lce_to_interval(sequence, k, lce):
    kmers = seq_to_kmers(sequence, k)
    indexes = set()
    for i, kmer in enumerate(kmers):
        if lce_motif(kmer, lce):
            for j in range(i, i+k):
                indexes.add(j)
    intervals = tools.ints_to_ranges(sorted(list(indexes)))
    return intervals


def display_lce(sequence, thresh_lce, k_lce):
    """Given a sequence, mark the motifs with a 'O'"""
    kmers_lce = seq_to_kmers(sequence, k_lce)
    indexes = set()
    for i, kmer in enumerate(kmers_lce):
        if lce_motif(kmer, thresh_lce):
            for item in range(i, i+k_lce):
                indexes.add(item)
    new_sequence = ''
    for i, let in enumerate(sequence):
        if i in indexes:
            new_sequence += 'O'
        else:
            new_sequence += let
    return new_sequence


def display_lca(sequence, alph_lca, k_lca):
    """Given a sequence, mark the motifs with a 'O'"""
    kmers_lca = seq_to_kmers(sequence, k_lca)
    indexes = set()
    for i, kmer in enumerate(kmers_lca):
        if lca_motif(kmer, alph_lca):
            for item in range(i, i+k_lca):
                indexes.add(item)
    new_sequence = ''
    for i, let in enumerate(sequence):
        if i in indexes:
            new_sequence += 'O'
        else:
            new_sequence += let
    return new_sequence


def main():
    pass


if __name__ == '__main__':
    main()


