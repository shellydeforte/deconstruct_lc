"""
Look at some specific examples of proteins
"""
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc


def main():
    """420-552"""
    k = 6
    lce = 1.6
    lca = 'SGEQAPDTNKR'
    seq = tools_fasta.fasta_to_seq('../laf1.fasta')[0]
    print(seq)
    disp = tools_lc.display_lc(seq, k, lca, lce)
    print(disp)
    print('')
    #seq = seq[0:420] + 'X'*132 + seq[552:]
    #print(seq)
    #disp = tools_lc.display_lc(seq, k, lca, lce)
    #print(disp)
    lc_count = tools_lc.count_lc_motifs(seq, k, lca, lce)
    print(lc_count)


if __name__ == '__main__':
    main()