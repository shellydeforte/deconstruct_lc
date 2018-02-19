import pandas as pd
from deconstruct_lc import tools_lc


def lce_example():
    examples = ['QQQQQQ', 'DDDYDD', 'NNNNRR', 'RERERE', 'PGAPPP', 'LLSSTS',
                'AADDFF', 'RQNGGG', 'SPESLL', 'LDELTI', 'GFKAPT']
    lces = []
    for example in examples:
        s_entropy = tools_lc.shannon(example)
        lces.append(s_entropy)
    df_dict = {'Shannon information entropy': lces, 'Example region': examples}
    df = pd.DataFrame(df_dict, columns=['Shannon information entropy', 'Example region'])
    return df


def main():
    df = lce_example()
    print(df)


if __name__ == '__main__':
    main()