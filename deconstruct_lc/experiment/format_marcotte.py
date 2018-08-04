import os
import pandas as pd

from deconstruct_lc import read_config
from deconstruct_lc import tools_fasta
from deconstruct_lc.scores.norm_score import NormScore


class FormatMarcotte(object):
    def __init__(self):
        config = read_config.read_config()
        data_dp = os.path.join(config['fps']['data_dp'])
        self.marcotte_fpi = os.path.join(data_dp, 'experiment', 'marcotte_puncta_proteins.xlsx')
        self.huh_fpi = os.path.join(data_dp, 'experiment', 'huh_cytoplasmic.xlsx')
        self.huh_full = os.path.join(data_dp, 'experiment', 'Huh_WK_et_al_2003_go_terms.txt')
        self.huh_fpout = os.path.join(data_dp, 'experiment', 'huh_scores.tsv')
        self.orf_trans = os.path.join(data_dp, 'proteomes', 'orf_trans.fasta')
        self.fpo = os.path.join(data_dp, 'experiment', 'marcotte_puncta_scores.tsv')
        self.manual = os.path.join(data_dp, 'experiment', 'cytoplasm_annotations_manual.txt')
        self.high_throughput = os.path.join(data_dp, 'experiment', 'cytoplasm_annotations_high_throughput.txt')
        self.computational = os.path.join(data_dp, 'experiment', 'cytoplasm_annotations_computational.txt')
        self.man_cytos = os.path.join(data_dp, 'experiment', 'cytosol_annotations_manual.txt')
        self.ht_cytos = os.path.join(data_dp, 'experiment', 'cytosol_annotations_high_throughput.txt')
        self.comp_cytos = os.path.join(data_dp, 'experiment', 'cytosol_annotations_computational.txt')
        self.marcott_not_huh = os.path.join(data_dp, 'experiment', 'marcotte_notin_huh.tsv')
        self.marcott_not_any = os.path.join(data_dp, 'experiment', 'marcotte_notcytoplasmic.tsv')

    def read_marcotte(self):
        df = pd.read_excel(self.marcotte_fpi, 'ST1')
        yeast_ids = list(df['ORF'])
        genes = list(df['Gene'])
        seqs = tools_fasta.get_yeast_seq_from_ids(self.orf_trans, yeast_ids)
        lengths = [len(seq) for seq in seqs]
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        df_out = pd.DataFrame({'Gene': genes, 'ORF': yeast_ids,
                               'LC Score': scores, 'Sequence': seqs,
                               'Length': lengths},
                              columns=['Gene', 'ORF', 'LC Score', 'Length', 'Sequence'])
        print(df_out.head())
        df_out.to_csv(self.fpo, sep='\t')

    def huh_out(self):
        huh_df = pd.read_excel(self.huh_fpi, sheetname='Sheet1')
        huh_ids = set(list(huh_df['Gene Systematic Name']))
        seqs, genes, orfs = tools_fasta.get_yeast_seq_gene_from_ids(self.orf_trans, huh_ids)
        lengths = [len(seq) for seq in seqs]
        ns = NormScore()
        scores = ns.lc_norm_score(seqs)
        df_out = pd.DataFrame({'Gene': genes, 'ORF': orfs,
                               'LC Score': scores, 'Sequence': seqs,
                               'Length': lengths},
                              columns=['Gene', 'ORF', 'LC Score', 'Length', 'Sequence'])
        print(df_out.head())
        df_out.to_csv(self.huh_fpout, sep='\t')

    def check_cytoplasmic(self):
        marcotte_df = pd.read_excel(self.marcotte_fpi, 'ST1')
        marcotte_ids = list(marcotte_df['ORF'])
        #huh_df = pd.read_excel(self.huh_fpi, 'Sheet1')
        #huh_ids = list(huh_df['Gene Systematic Name'])
        huh_df = pd.read_csv(self.huh_full, sep='\t')
        huh_ids = set(list(huh_df['Gene Systematic Name']))
        huh_genes = set(list(huh_df['Gene']))
        clean_huh = []
        for item in huh_ids:
            clean_huh.append(item[0:7])
        print(clean_huh[0:10])
        print(len(set(clean_huh)))
        not_in = set(marcotte_ids) - set(clean_huh)
        not_in_df = marcotte_df[marcotte_df['ORF'].isin(not_in)]
        not_genes = not_in_df['Gene']
        print(len(set(not_genes) - set(huh_genes)))
        print(len(not_in))
        #not_in_df.to_csv(self.marcott_not_huh, sep='\t')

    def check_annotations(self):
        marcotte_df = pd.read_excel(self.marcotte_fpi, 'ST1')
        marcotte_ids = list(marcotte_df['ORF'])
        manual_df = pd.read_csv(self.manual, sep='\t', header=7)
        manual_ids = list(manual_df['Gene Systematic Name'])
        ht_df = pd.read_csv(self.high_throughput, sep='\t', header=7)
        ht_ids = list(ht_df['Gene Systematic Name'])
        comp_df = pd.read_csv(self.computational, sep='\t', header=7)
        comp_ids = list(comp_df['Gene Systematic Name'])
        man_cyt_df = pd.read_csv(self.man_cytos, sep='\t', header=7)
        man_cyt_ids = list(man_cyt_df['Gene Systematic Name'])
        ht_cyt_df = pd.read_csv(self.ht_cytos, sep='\t', header=7)
        ht_cyt_ids = list(ht_cyt_df['Gene Systematic Name'])
        comp_cyt_df = pd.read_csv(self.comp_cytos, sep='\t', header=7)
        comp_cyt_ids = list(comp_cyt_df['Gene Systematic Name'])
        all_ids = set(manual_ids + ht_ids + comp_ids + man_cyt_ids + ht_cyt_ids)
        print(len(all_ids))
        not_in = set(marcotte_ids) - all_ids
        print(len(set(marcotte_ids) - all_ids))
        marcotte_out = marcotte_df[marcotte_df['ORF'].isin(not_in)]
        marcotte_out.to_csv(self.marcott_not_any, sep='\t')


def main():
    fm = FormatMarcotte()
    fm.check_cytoplasmic()


if __name__ == '__main__':
    main()