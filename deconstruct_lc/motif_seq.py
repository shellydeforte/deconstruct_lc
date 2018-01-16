from deconstruct_lc import tools_lc

class LcSeq(object):

    def __init__(self, sequence, k, lc, lctype):
        self.sequence = sequence
        self.k = k
        self.lc = lc
        self.lctype = lctype

    def overlapping_seq_in_motif(self):
        """Counts overlapping composition"""
        inside_seq = ''
        outside_seq = ''
        kmers = tools_lc.seq_to_kmers(self.sequence, self.k)
        for kmer in kmers:
            if self.lctype == 'lca':
                if tools_lc.lca_motif(kmer, self.lc):
                    inside_seq += kmer
                else:
                    outside_seq += kmer
            elif self.lctype == 'lce':
                if tools_lc.lce_motif(kmer, self.lc):
                    inside_seq += kmer
                else:
                    outside_seq += kmer
            else:
                raise Exception("lctype must be lca or lce")
        return inside_seq, outside_seq

    def list_motifs(self):
        motifs = []
        kmers = tools_lc.seq_to_kmers(self.sequence, self.k)
        for kmer in kmers:
            if self.lctype == 'lca':
                if tools_lc.lca_motif(kmer, self.lc):
                    motifs.append(kmer)
            elif self.lctype == 'lce':
                if tools_lc.lce_motif(kmer, self.lc):
                    motifs.append(kmer)
        return motifs

    def seq_in_motif(self):
        ind_in, ind_out = self._get_motif_indexes()
        seq_in = ''.join([self.sequence[i] for i in ind_in])
        seq_out = ''.join([self.sequence[i] for i in ind_out])
        return seq_in, seq_out

    def _get_motif_indexes(self):
        kmers = tools_lc.seq_to_kmers(self.sequence, self.k)
        ind_in = set()
        for i, kmer in enumerate(kmers):
            if self.lctype == 'lca':
                if tools_lc.lca_motif(kmer, self.lc):
                    for j in range(i, i + self.k):
                        ind_in.add(j)
            elif self.lctype == 'lce':
                if tools_lc.lce_motif(kmer, self.lc):
                    for j in range(i, i + self.k):
                        ind_in.add(j)
            else:
                raise Exception("lctype must be lca or lce")
        ind_out = set(list(range(len(self.sequence)))) - ind_in
        return ind_in, ind_out

def main():
    k = 6
    lc = 'SEQAPDTNKR'
    lctype = 'lca'
    #k = 6
    #lc = 1.6
    #lctype = 'lce'
    seq = 'RSQLTSLEKDCSLRAIEKNDDNSCRNPEHTDVIDELEEEEDIDTK'
    print(seq[2])
    ls = LcSeq(seq, k, lc, lctype)
    ind_in, ind_out = ls._get_motif_indexes()
    print(ind_in)
    print(ind_out)
    seq_in, seq_out = ls.seq_in_motif()
    print(seq_in)
    print(seq_out)



if __name__ == '__main__':
    main()