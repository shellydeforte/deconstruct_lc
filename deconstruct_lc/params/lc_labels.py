from itertools import combinations
from deconstruct_lc import tools_lc


class GetLabels(object):
    def __init__(self, k):
        self.k = k

    def format_labels(self, lcs):
        labels = ['{}_{}'.format(self.k, lc) for lc in lcs]
        return labels

    def create_lcas(self, alph):
        """
        alph should usually be 'SGEQAPDTNKRL'
        Return all combinations of the base LCA from 2 to the full length as
        a list of strings.

        Example: create_lcas('SGE') returns ['k_SG', 'k_SE', 'k_GE', 'Sk_GE']
        """
        lcas = []
        for i in range(2, len(alph) + 1):
            lca_combos = combinations(alph, i)
            for lca in lca_combos:
                lca_str = ''.join(lca)
                lcas.append(lca_str)
        lca_labels = self.format_labels(lcas)
        return lca_labels

    def find_entropy_values(self, all_seqs):
        """
        Return all the possible shannon entropies in my data set for the given
        k-mer length, rounded up to the nearest 0.1.
        """
        all_shannon = set()
        for seq in all_seqs:
            kmers = tools_lc.seq_to_kmers(seq, self.k)
            for kmer in kmers:
                s = tools_lc.shannon(kmer)
                all_shannon.add(s)
        new_scores = []
        all_shannon = sorted(list(all_shannon), reverse=True)
        for score in all_shannon[1:]:
            new_score = self._round_up(score)
            new_scores.append(new_score)
            new_scores = sorted(list(set(new_scores)), reverse=True)
        lce_labels = self.format_labels(new_scores)
        return lce_labels

    def _round_up(self, score):
        """
        Given a float, round to the nearest 10th place, and then format into
        string.
        """
        new_score = round(score, 1)
        if new_score < score:
            new_score = new_score + 0.1
            new_score = round(new_score, 1)  # trailing 0.999..
        return new_score


def main():
    gl = GetLabels(4)
    lcas = gl.create_lcas('SGE')
    print(lcas)


if __name__ == '__main__':
    main()