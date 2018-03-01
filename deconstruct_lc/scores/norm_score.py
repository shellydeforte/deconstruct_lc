from deconstruct_lc import tools_lc


class NormScore(object):
    def __init__(self):
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'
        self.lc_m = 0.06744064704548541
        self.lc_b = 16.5

    def lc_norm_score(self, seqs):
        lens = [len(seq) for seq in seqs]
        scores = tools_lc.calc_lc_motifs(seqs, self.k, self.lca, self.lce)
        lc_norm = self.norm_function(self.lc_m, self.lc_b, scores, lens)
        return lc_norm

    def lc_miss_norm(self, seqs, miss_seqs):
        """For PDB chains, do not count the kmer if it has a missing residue"""
        lens = [len(seq) for seq in seqs]
        scores = tools_lc.calc_lc_motifs_nomiss(seqs, miss_seqs, self.k,
                                                self.lca, self.lce)
        lc_norm = self.norm_function(self.lc_m, self.lc_b, scores, lens)
        return lc_norm

    @staticmethod
    def norm_function(m, b, raw_scores, lengths):
        norm_scores = []
        for raw_score, length in zip(raw_scores, lengths):
            norm_score = raw_score - ((m * length) + b)
            norm_scores.append(norm_score)
        return norm_scores


def main():
    pass


if __name__ == '__main__':
    main()