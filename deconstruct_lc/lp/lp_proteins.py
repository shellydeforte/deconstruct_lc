from deconstruct_lc.scores.norm_score import NormScore
from deconstruct_lc import tools_fasta
from deconstruct_lc import tools_lc


class CheckPrDs(object):
    def __init__(self):
        self.fasta_fpi = 'C:\LP\lps_proteins.fasta'
        self.k = 6
        self.lce = 1.6
        self.lca = 'SGEQAPDTNKR'

    def read_fasta(self):
        pids, seqs = tools_fasta.fasta_to_id_seq(self.fasta_fpi)
        norm = NormScore()
        # ent1[211:457]
        ent1 = seqs[0]
        #print(ent1[211:457])
        ent1wo = ent1[:211] + ent1[457:]
        #print(norm.lc_norm_score([ent1wo]))
        #print(norm.lc_norm_score([ent1]))
        # ent2[224:616]
        ent2 = seqs[1]
        #print(ent2[224:616])
        ent2wo = ent2[:224] + ent2[616:]
        #print(norm.lc_norm_score([ent2wo]))
        # yap1801[351:638]
        yap1801 = seqs[2]
        #print(yap1801[351:638])
        yap1801wo = yap1801[:351] + yap1801[638:]
        #print()
        #print(norm.lc_norm_score([yap1801]))
        #print(norm.lc_norm_score([yap1801wo]))
        # yap1802[319:569]
        yap1802 = seqs[3]
        #print(yap1802[319:569])
        yap1802wo = yap1802[:319] + yap1802[569:]
        #print(norm.lc_norm_score([yap1802wo]))
        # sla1[954:1244]
        sla1 = seqs[4]
        print(len(sla1))
        print(sla1[954:1244])
        print()
        ns = tools_lc.display_lc(sla1, self.k, self.lca, self.lce)
        print(sla1)
        print(ns)
        sla1wo = sla1[:954] + sla1[1244:]
        print(norm.lc_norm_score([sla1wo]))
        print(norm.lc_norm_score([sla1]))
        #sla2[348:442]
        sla2 = seqs[5]
        #print(sla2[348:442])
        sla2wo = sla2[:348] + sla2[442:]
        #print(norm.lc_norm_score([sla2wo]))
        #print(norm.lc_norm_score([sla2]))
        # sup35[0:123]
        sup35 = seqs[6]
        #print(sup35[0:123])

        #print(norm.lc_norm_score([seq]))



def main():
    cp = CheckPrDs()
    cp.read_fasta()


if __name__ == '__main__':
    main()