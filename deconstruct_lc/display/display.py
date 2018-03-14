import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc.scores import norm_score
from deconstruct_lc import tools_lc
from deconstruct_lc import tools

class Display(object):
    """
    input a sequence, a list of sequences, or a fasta file.
    Create an html file that will put the LC motifs in bold
    Create another html file that will highlight charge
    Create another html file that will highlight Q/N
    Create another html file that will highlight aromatics
    Output score
    1. Just write html page with sequences and score
    """
    def __init__(self):
        config = read_config.read_config()
        self.data_dp = config['fps']['data_dp']
        self.bc_dp = os.path.join(self.data_dp, 'bc_analysis', 'Cytoplasmic_Stress_Granule_score.tsv')
        self.fp_out = os.path.join(self.data_dp, 'display', 'stressgranule_yeast.html')
        self.k = config['score'].getint('k')
        self.lca = config['score'].get('lca')
        self.lce = config['score'].getfloat('lce')

    def write_body(self):
        contents = '''
        <!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
        <html>

        <head>
          <meta content="text/html; charset=ISO-8859-1"
         http-equiv="content-type">
          <title>Hello</title>
        </head>
        <body>
        '''
        form_seqs = self.read_seq()
        for seq in form_seqs:
            contents += seq
            contents += '<br>'
        contents += '''
        </body>
        </html>
        '''
        with open(self.fp_out, 'w') as fo:
            fo.write(contents)

    def read_seq(self):
        df = pd.read_csv(self.bc_dp, sep='\t', index_col=0)
        df = df[df['Organism'] == 'YEAST']
        seqs = df['Sequence']
        form_seqs = []
        for seq in seqs:
            inds = tools_lc.lc_to_indexes(seq, self.k, self.lca, self.lce)
            ranges = list(tools.ints_to_ranges(sorted(list(inds))))
            es = self.format_string(seq, ranges)
            form_seqs.append(es)
        return form_seqs

    def format_string(self, seq, c_ind):
        """
        Given a sequence, return the excel formatted color string of the form:
         'red, 'sequence region', 'sequence region2', red, sequence region3...
        """
        es = ''  # excel string
        if len(c_ind) > 0:
            if c_ind[0][0] != 0:
                es += seq[0:c_ind[0][0]]
            for i, index in enumerate(c_ind):
                es += '<b>' + seq[index[0]:index[1] + 1] + '</b>'
                if i != len(c_ind) - 1:
                    es += seq[index[1] + 1:c_ind[i + 1][0]]
            es += seq[c_ind[-1][1] + 1:]
        else:
            es = seq
        return es


def main():
    d = Display()
    d.write_body()


if __name__ == '__main__':
    main()
