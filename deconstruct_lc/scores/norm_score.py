from deconstruct_lc import tools_lc

class NormScore(object):
    def __init__(self, seqs):
        self.seqs = seqs
        self.lens = [len(seq) for seq in seqs]
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'
        self.lc_m = 0.066213297264721263
        self.lc_b = 16.5

    def lc_norm_score(self):
        scores = tools_lc.calc_lc_motifs(self.seqs, self.k, self.lca, self.lce)
        lc_norm = self.norm_function(self.lc_m, self.lc_b, scores, self.lens)
        return lc_norm

    def lc_miss_norm(self, miss_seqs):
        """For PDB chains, do not count the kmer if it has a missing residue"""
        scores = tools_lc.calc_lc_motifs_nomiss(self.seqs, miss_seqs, self.k,
                                                self.lca, self.lce)
        lc_norm = self.norm_function(self.lc_m, self.lc_b, scores, self.lens)
        return lc_norm

    def norm_function(self, m, b, raw_scores, lengths):
        norm_scores = []
        for raw_score, length in zip(raw_scores, lengths):
            norm_score = raw_score - ((m * length) + b)
            norm_scores.append(norm_score)
        return norm_scores


def main():
    pass


if __name__ == '__main__':
    main()