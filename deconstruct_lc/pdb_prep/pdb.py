"""
Create processed PDB files, ready for CD-hit
All have the following filtering applied:
x-ray
eukaryote
standard amino acid alphabet
>= 100, <= 2000

nomiss: no missing regions
miss: include missing regions

This is another file that I probably do not need to include in the final
scripts. I'm not sure, I need to think about it. If I do, I definitely need
to do some cleanup
"""
import os
import pickle
from Bio import SeqIO
#from predict_llps.lib.config import Config

config = Config()
data_fp = config.DATA


def create_pdb_miss():
    """
    Apply the following filtering:
    x-ray
    eukaryote
    standard amino acid alphabet
    >= 100, <= 2000
    do not apply any missing regions filtering
    Run on July 10, 2017
    """
    out_fp = os.path.join(data_fp, 'pdb', 'pdb_miss.fasta')
    pdb_seq = os.path.join(data_fp, 'pdb', 'ss_dis_seqs.fasta')
    diffraction = get_diffraction()
    eukaryote = euk_pdb()
    new_records = []
    with open(pdb_seq, 'r') as seq_fi:
        for seq_rec in SeqIO.parse(seq_fi, 'fasta'):
            pdb_chain = id_cleanup(str(seq_rec.id))
            pdb = pdb_chain.split('_')[0]
            sequence = str(seq_rec.seq)
            prot_len = len(sequence)
            if pdb_chain in eukaryote:
                if pdb in diffraction:
                    if standard_aa(sequence):
                        if prot_len >= 100 and prot_len <= 2000:
                            new_records.append(seq_rec)
    with open(out_fp, 'w') as seq_fo:
        SeqIO.write(new_records, seq_fo, 'fasta')
    count = 0
    with open(out_fp, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                count += 1
    print('There are {} records'.format(count))


def create_pdb_nomiss():
    """
    Apply the following filtering:
    x-ray
    eukaryote
    standard amino acid alphabet
    >= 100, <= 2000
    no missing regions
    Run on July 10, 2017
    """
    out_fp = os.path.join(data_fp, 'pdb', 'pdb_nomiss.fasta')
    pdb_seq = os.path.join(data_fp, 'pdb', 'ss_dis_seqs.fasta')
    diffraction = get_diffraction()
    eukaryote = euk_pdb()
    pdb_miss = get_missing()
    new_records = []
    with open(pdb_seq, 'r') as seq_fi:
        for seq_rec in SeqIO.parse(seq_fi, 'fasta'):
            pdb_chain = id_cleanup(str(seq_rec.id))
            pdb = pdb_chain.split('_')[0]
            sequence = str(seq_rec.seq)
            prot_len = len(sequence)
            if pdb_chain in eukaryote:
                if pdb in diffraction:
                    if standard_aa(sequence):
                        if prot_len >= 100 and prot_len <= 2000:
                            if pdb_chain in pdb_miss: # Check for missing regions
                                if pdb_miss[pdb_chain] == 0:
                                    new_records.append(seq_rec)
    with open(out_fp, 'w') as seq_fo:
        SeqIO.write(new_records, seq_fo, 'fasta')
    count = 0
    with open(out_fp, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                count += 1
    print('There are {} records'.format(count))


def create_pdb_miss_dict():
    """Create a dictionary that is pdb_chain: # missing regions, and pickle"""
    pdb_miss = get_missing()
    out_fp = os.path.join(data_fp, 'pdb', 'pdb_miss_dict.pkl')
    pickle.dump(pdb_miss, open(out_fp, 'wb'))


def get_missing():
    """Create a dictionary with pdb_chain: # missing regions"""
    pdb_miss = {}
    pdb_dis = os.path.join(data_fp, 'pdb', 'ss_dis_dis.fasta')
    with open(pdb_dis, 'r') as dis_fo:
        for dis_rec in SeqIO.parse(dis_fo, 'fasta'):
            pdb_miss[id_cleanup(str(dis_rec.id))] = missing_count(str(dis_rec.seq))
    return pdb_miss


def missing_count(dis_seq):
    return dis_seq.count('X')


def get_diffraction():
    """
    Read diffraction file and return a set of PDB IDs that are
    from crystal structures. Note this does not include chain info. ie.
    {'3S5T', '1F8W', '1PCX',...
    """
    type_fp = os.path.join(data_fp, 'pdb', 'pdb_entry_type.txt')
    xray_ids = []
    with open(type_fp) as fi:
        for line in fi:
            ls = line.split()
            pdb_id = ls[0]
            method = ls[2]
            if method == 'diffraction':
                xray_ids.append(pdb_id.upper())
    return set(xray_ids)


def euk_pdb():
    """
    Create a list of PDB IDs that are eukaryotes, ie.
    {'3V6M_A', '3WGU_E', '4ITZ_B',...
    """
    euks = euk_list()
    euk_pdbs = []
    with open(os.path.join(data_fp, 'pdb', 'pdb_chain_taxonomy.tsv')) as fi:
        next(fi)
        next(fi)
        for line in fi:
            ls = line.split()
            if ls[2] in euks:
                pdb_chain = '{}_{}'.format(ls[0].upper(), ls[1].upper())
                euk_pdbs.append(pdb_chain)
    return set(euk_pdbs)


def euk_list():
    """
    Create a set of all taxonomic identifiers that are 'E' for eukaryote
    ie. {'348046', '160085', '143180',...
    """
    tax_org = read_speclist()
    euks = []
    for item in tax_org:
        if tax_org[item] == 'E':
            euks.append(item)
    return set(euks)


def read_speclist():
    """
    Read spec list and create a tax: org dictionary, ie. '3320': 'E'
    """
    tax_org = {}
    with open(os.path.join(data_fp, 'speclist.txt')) as fi:
        for i in range(59):
            next(fi)
        for line in fi:
            ls = line.split()
            if len(ls) >= 3:
                if ls[2][-1] == ':':
                    tax = ls[2][0:-1]
                    org = ls[1]
                    tax_org[tax] = org
    return tax_org


def standard_aa(sequence):
    aas = 'ADKERNTSQYFLIVMCWHGP'
    for c in sequence:
        if c not in aas:
            return False
    return True


def id_cleanup(protein_id):
    if '|' in protein_id:
        nid = protein_id.split('|')[1]
    elif ':' in protein_id:
        nid = '{}_{}'.format(protein_id.split(':')[0], protein_id.split(':')[1])
    else:
        nid = protein_id
    return nid


def main():
    create_pdb_miss_dict()


if __name__ == '__main__':
    main()