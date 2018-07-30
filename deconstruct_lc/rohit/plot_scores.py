import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc.display.display_lc import Display

class PlotRohit(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.dp = os.path.join(data_dp, 'FUS')

    def read_fasta(self):
        ns = NormScore()
        files = os.listdir(self.dp)
        files = [afile for afile in files if '.fasta' in afile]
        seqs = []
        sizes = []
        for afile in files:
            seq = tools_fasta.fasta_to_seq(os.path.join(self.dp, afile))
            seqs.append(seq[0])
        lc_scores = ns.lc_norm_score(seqs)
        for score in lc_scores:
            sizes.append(score*.8)
        labels = []
        for file, score in zip(files, lc_scores):
            labels.append(file.split('_')[0] + ' ' + str(score).split('.')[0])
        arg, tyr = self.arg_tyr(seqs)
        df = pd.DataFrame({'Num_Arg': arg, 'Num_Tyr': tyr, 'group': labels})
        p1 = sns.regplot(data=df, x='Num_Tyr', y='Num_Arg', fit_reg=False,
                    color="skyblue", scatter_kws={"s": sizes})
        for line in range(0, df.shape[0]):
            p1.text(df.Num_Tyr[line] + 0.2, df.Num_Arg[line], df.group[line],
                    horizontalalignment='left', size='medium', color='black',
                    weight='semibold')
        #sns.plt.show()
        fasta_in = os.path.join(self.dp, 'HNRNPA1_P09651.fasta')
        fn_out = os.path.join(self.dp, 'HNRNPA1.html')
        dis = Display(fasta_in, fn_out, color=True)
        dis.write_body()

        #fig, ax = plt.subplots()
        #for i, txt in enumerate(labels):
        #    ax.annotate(txt, (tyr[i], arg[i]))
        #ax.scatter(tyr, arg, alpha=0.5, sizes=sizes)
        #plt.show()

    def arg_tyr(self, seqs):
        arg = []
        tyr = []
        for seq in seqs:
            arg.append(seq.count('R'))
            tyr.append(seq.count('Y'))
        return arg, tyr

def main():
    pr = PlotRohit()
    pr.read_fasta()


if __name__ == '__main__':
    main()