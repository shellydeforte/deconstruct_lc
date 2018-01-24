from deconstruct_lc import tools_lc

class CalcScore(object):
    def __init__(self, seqs):
        self.seqs = seqs
        self.lens = [len(seq) for seq in seqs]
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'
        self.m = 0.066213297264721263
        self.b = 1.7520712972708843

    def lc_norm_score(self):
        scores = tools_lc.calc_lc_motifs(self.seqs, self.k, self.lca, self.lce)
        lc_norm = self.norm_function(self.m, self.b, scores, self.lens)
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