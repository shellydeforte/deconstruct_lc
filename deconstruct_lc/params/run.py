from deconstruct_lc import read_config
from deconstruct_lc.params import raw_scores
from deconstruct_lc.params import raw_svm
from deconstruct_lc.params import raw_top
from deconstruct_lc.params import write_mb
from deconstruct_lc.params import raw_norm
from deconstruct_lc.params import norm_svm
from deconstruct_lc.params import ran_forest


class RunRaw(object):
    def __init__(self):
        self.config = read_config.read_config()

    def run_raw_scores(self):
        """
        Write files with the un-normalized LCA/LCE sums for each possible
        LCA and LCE for each value of k
        """
        pr = raw_scores.RawScores(self.config)
        pr.write_lca()
        pr.write_lce()

    def run_svm(self):
        """
        Write the accuracy as obtained by an SVC on the un-normalized LC sums
        """
        rs = raw_svm.RawSvm(self.config)
        rs.svm_lca_lce()

    def run_rawtop(self):
        """
        Write only those k, LCA, LCE combos with an un-normalized accuracy > 0.82
        """
        rm = raw_top.RawTop(self.config)
        rm.write_top()

    def run_writemb(self):
        """
        After selecting a representative set of LCA proteins by hand, calculate
        the normalization parameters both for the individual LCA/LCE values,
        but also in combinations
        """
        mb = write_mb.WriteMb(self.config)
        mb.write_mb_solo()
        mb.write_mb_combos()

    def run_rawnorm(self):
        """
        Write the normalized scores based on the calculated normalization
        parameters, when pearson's correlation coefficient was > 0.7
        """
        rn = raw_norm.RawNorm(self.config)
        rn.solo_norm()
        rn.combo_norm()

    def run_normsvm(self):
        """
        Run an SVM classifier on the normalized scores for LCA/LCE and combinations
        """
        ns = norm_svm.NormSvm(self.config)
        ns.oned_svm()

    def run_ran_forest(self):
        """
        Run random forest with all parameters, and in various combinations
        to check for best parameters and upper cap
        """
        rf = ran_forest.BestFeatures(self.config)
        rf.ran_forest()


def main():
    rr = RunRaw()
    rr.run_ran_forest()


if __name__ == '__main__':
    main()