import os
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from itertools import combinations

from deconstruct_lc import read_config
from deconstruct_lc.chi2 import create_contingency


def format_cont(row1, row2):
    return np.array([list(row1), list(row2)])


def write_bc_vs_bg(bcs, fpo, organism):
    adict = {'BC': [], 'pval': []}
    for bc in bcs:
        bcc = create_contingency.BcProteome(bc, organism)
        bc_cont, org_cont = bcc.get_cont_table()
        if sum(list(bc_cont)) > 10:
            ct = format_cont(bc_cont, org_cont)
            pval = chi2_contingency(ct)[1]
            adict['BC'].append(bc)
            adict['pval'].append(pval)
    df_out = pd.DataFrame(adict, columns=['BC', 'pval'])
    df_out.to_csv(fpo, sep='\t')


def run_bc_vs_bg():
    bcs = ['Nuclear_Stress_Granule', 'Nucleolus', 'P_granule', 'PDB',
           'PML_Body', 'P_Body', 'Nuclear_Speckles', 'Cytoplasmic_Stress_Granule',
           'Cajal_bodies', 'Paraspeckle', 'Centrosome']
    # Compare BC to background
    config = read_config.read_config()
    data_dp = config['fps']['data_dp']
    yeast_fpo = os.path.join(data_dp, 'chi2', 'bc_vs_yeast.tsv')
    human_fpo = os.path.join(data_dp, 'chi2', 'bc_vs_human.tsv')
    write_bc_vs_bg(bcs, yeast_fpo, 'Yeast')
    write_bc_vs_bg(bcs, human_fpo, 'Human')


def write_bc_vs_bc(bcs, fpo, organism):
    combs = list(combinations(bcs, 2))
    adict = {'BC': [], 'pval': []}
    for comb in combs:
        bc1 = comb[0]
        bc2 = comb[1]
        bcc1 = create_contingency.BcProteome(bc1, organism)
        bc1_cont, org_cont = bcc1.get_cont_table()
        bcc2 = create_contingency.BcProteome(bc2, organism)
        bc2_cont, org_cont = bcc2.get_cont_table()
        if sum(list(bc1_cont)) > 10 and sum(list(bc2_cont)) > 10:
            ct = format_cont(bc1_cont, bc2_cont)
            pval = chi2_contingency(ct)[1]
            adict['BC'].append(comb)
            adict['pval'].append(pval)
    df_out = pd.DataFrame(adict, columns=['BC', 'pval'])
    df_out.to_csv(fpo, sep='\t')


def run_bc_vs_bc():
    bcs = ['Nuclear_Stress_Granule', 'Nucleolus', 'P_granule', 'PDB',
           'PML_Body', 'P_Body', 'Nuclear_Speckles', 'Cytoplasmic_Stress_Granule',
           'Cajal_bodies', 'Paraspeckle', 'Centrosome']
    # Compare BC to each other
    config = read_config.read_config()
    data_dp = config['fps']['data_dp']
    yeast_fpo = os.path.join(data_dp, 'chi2', 'bc_vs_bc_yeast.tsv')
    human_fpo = os.path.join(data_dp, 'chi2', 'bc_vs_bc_human.tsv')
    write_bc_vs_bc(bcs, yeast_fpo, 'Yeast')
    write_bc_vs_bc(bcs, human_fpo, 'Human')


def main():
    bcs = ['Nuclear_Stress_Granule', 'Nucleolus', 'P_granule', 'PDB',
           'PML_Body', 'P_Body', 'Nuclear_Speckles', 'Cytoplasmic_Stress_Granule',
           'Cajal_bodies', 'Paraspeckle', 'Centrosome']
    # Compare BC to each other
    run_bc_vs_bc()






if __name__ == '__main__':
    main()