import os

from deconstruct_lc import read_config
from deconstruct_lc.scores import norm_score
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc
from deconstruct_lc import tools

class Display(object):
    def __init__(self, fn_out, color=True):
        config = read_config.read_config()
        data_dp = config['fps']['data_dp']
        self.k = config['score'].getint('k')
        self.lca = config['score'].get('lca')
        self.lce = config['score'].getfloat('lce')
        self.fp_out = os.path.join(data_dp, 'display', fn_out)
        self.color = color

    def write_body(self, headers, seqs):
        contents = '''
        <!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
        <html>

        <head>
          <meta content="text/html; charset=ISO-8859-1"
         http-equiv="content-type">
          <title>LC motifs</title>
        </head>
        <body>
        '''
        form_seqs = self.form_seq(seqs)
        ns = norm_score.NormScore()
        scores = ns.lc_norm_score(seqs)
        zip_list = zip(scores, headers, form_seqs)
        key_fun = lambda pair: pair[0]
        sorted_tup = sorted(zip_list, reverse=True, key=key_fun)
        for score, head, seq in sorted_tup:
            contents += head
            contents += '<br>'
            contents += seq
            contents += '<br>'
            contents += str(score)
            contents += '<br><br>'
        contents += '''
        </body>
        </html>
        '''
        with open(self.fp_out, 'w') as fo:
            fo.write(contents)

    def form_seq(self, seqs):
        form_seqs = []
        for seq in seqs:
            inds = tools_lc.lc_to_indexes(seq, self.k, self.lca, self.lce)
            ranges = list(tools.ints_to_ranges(sorted(list(inds))))
            es = self.format_string(seq, ranges)
            if self.color:
                ns = self.add_colors(es)
                form_seqs.append(ns)
            else:
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

    def add_colors(self, form_seq):
        """
        Given a string that has been formatted with bold, add color tags
        ED: blue, RK: red, ST: green QN: orange, AG: just leave black, P: brown
        """
        ns = ''
        for i, aa in enumerate(form_seq):
            if aa == 'S':
                if i+10 > len(form_seq):
                    if 'S' in form_seq[i+2:]:
                        ns += '<font color=\'blue\'>'+ '<b>' + aa + '</b>' + '</font>'
                    elif 'S' in form_seq[i-10:i-2]:
                        ns += '<font color=\'blue\'>' + '<b>' + aa + '</b>' + '</font>'
                    else:
                        ns += '<font color=\'red\'>' + '<b>' + aa + '</b>' + '</font>'
                elif i-10 < 0:
                    if 'S' in form_seq[i+2:i+10]:
                        ns += '<font color=\'blue\'>' + '<b>' + aa + '</b>' + '</font>'
                    elif 'S' in form_seq[0:i-2]:
                        ns += '<font color=\'blue\'>' + '<b>' + aa + '</b>' + '</font>'
                    else:
                        ns += '<font color=\'red\'>' + '<b>' + aa + '</b>' + '</font>'
                else:
                    if 'S' in form_seq[i+2:i+10]:
                        ns += '<font color=\'blue\'>' + '<b>' + aa + '</b>' + '</font>'
                    elif 'S' in form_seq[i-10:i-2]:
                        ns += '<font color=\'blue\'>' + '<b>' + aa + '</b>' + '</font>'
                    else:
                        ns += '<font color=\'red\'>' + '<b>' + aa + '</b>' + '</font>'
            else:
                ns += aa
        return ns


def main():
    d = Display('temp')
    seq = 'MDEPPLAQPLELNQHSRFIIGSVSEDNSEDEISNLVKLDLLEEKEGSLSPASVGSDTLSDLGISSLQDGLALHIRSSMS'
    ns = d.add_colors(seq)
    print(ns)


if __name__ == '__main__':
    main()