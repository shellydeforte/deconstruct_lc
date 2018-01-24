from deconstruct_lc import tools_lc

class NormScore(object):
    def __init__(self, seqs):
        self.seqs = seqs
        self.lens = [len(seq) for seq in seqs]
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'
        self.lc_m = 0.066213297264721263
        #self.lc_b = 1.7520712972708843 + 15
        self.lc_b = 16.5
        self.lca_m = 0.049569081348185169
        self.lca_b = 1.5446026962158523
        self.lce_m = 0.022592872067308627
        self.lce_b = 0.37785221922942025

    def lc_norm_score(self):
        scores = tools_lc.calc_lc_motifs(self.seqs, self.k, self.lca, self.lce)
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