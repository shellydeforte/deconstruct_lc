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
        self.bc_dp = os.path.join(self.data_dp, 'bc_analysis', 'P_Body_score.tsv')
        self.fp_out = os.path.join(self.data_dp, 'display', 'Pbody_human.html')
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
        form_seqs, scores = self.read_seq()
        sort_scores, sort_form_seqs = tools.sort_list_by_second_list(scores, form_seqs)
        for seq, score in zip(sort_form_seqs, sort_scores):
            contents += seq
            contents += '<br>'
            contents += str(score)
            contents += '<br>'
        contents += '''
        </body>
        </html>
        '''
        with open(self.fp_out, 'w') as fo:
            fo.write(contents)

    def read_seq(self):
        df = pd.read_csv(self.bc_dp, sep='\t', index_col=0)
        df = df[df['Organism'] == 'HUMAN']
        seqs = df['Sequence']
        scores = list(df['LC Score'])
        form_seqs = []
        for seq in seqs:
            inds = tools_lc.lc_to_indexes(seq, self.k, self.lca, self.lce)
            ranges = list(tools.ints_to_ranges(sorted(list(inds))))
            es = self.format_string(seq, ranges)
            ns = self.add_colors(es)
            #ns = self.color_aromatics(es)
            form_seqs.append(ns)
        return form_seqs, scores

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

    def add_colors(self, form_seq):
        """
        Given a string that has been formatted with bold, add color tags
        ED: blue, RK: red, ST: green QN: orange, AG: just leave black, P: brown
        """
        ns = ''
        for aa in form_seq:
            if aa == 'E' or aa == 'D':
                ns += '<font color=\'blue\'>' + aa + '</font>'
            elif aa == 'R' or aa == 'K':
                ns += '<font color=\'red\'>' + aa + '</font>'
            elif aa == 'Q' or aa == 'N':
                ns += '<font color=\'orange\'>' + aa + '</font>'
            elif aa == 'S' or aa == 'T':
                ns += '<font color=\'green\'>' + aa + '</font>'
            elif aa == 'P':
                ns += '<font color=\'brown\'>' + aa + '</font>'
            else:
                ns += aa
        return ns

    def color_aromatics(self, form_seq):
        ns = ''
        for aa in form_seq:
            if aa == 'Y' or aa == 'F' or aa == 'W':
                ns += '<font color=\'blue\'>' + aa + '</font>'
            elif aa == 'R' or aa == 'K':
                ns += '<font color=\'red\'>' + aa + '</font>'
            else:
                ns += aa
        return ns


def main():
    seq = 'MHQQHSKSENKPQQQRKKFEGPKREAILDLAKYKDSKIRVKLMGGKLVIGVLKGYDQLMNLVLDDTVEYMSNPDDENNTELISKNARKLGLTVIRGTILVSLSSAEGSDVLYMQK'
    d = Display()
    d.write_body()
    #inds = tools_lc.lc_to_indexes(seq, d.k, d.lca, d.lce)
    #ranges = list(tools.ints_to_ranges(sorted(list(inds))))
    #es = d.format_string(seq, ranges)
    #ns = d.add_colors(es)
    #print(ns)


if __name__ == '__main__':
    main()
