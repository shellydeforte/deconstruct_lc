import numpy as np

class KappaKmers(object):
    #...................................................................................#
    def __init__(self, kmers, seq):
        """
        seq = amino acid sequence as a string
        """
        self.seq = seq
        self.len = len(seq)
        self.kmers = kmers
        self.kmer_charges = self.kmerChargePattern()
        self.chargePattern = self.seqChargePattern()

    def seqChargePattern(self):
        charges = {'K': 1, 'R': 1, 'D': -1, 'E': -1}
        nseq = []
        for aa in self.seq:
            if aa in charges:
                nseq.append(charges[aa])
            else:
                nseq.append(0)
        return np.array(nseq)

    def kmerChargePattern(self):
        charges = {'K': 1, 'R': 1, 'D': -1, 'E': -1}
        kmer_charges = []
        for kmer in self.kmers:
            nkmer = []
            for aa in kmer:
                if aa in charges:
                    nkmer.append(charges[aa])
                else:
                    nkmer.append(0)
            kmer_charges.append(nkmer)
        return np.array(kmer_charges)

    #...................................................................................#
    def deltaForm(self):
        """ Calculate the delta value as defined in REF 1
        """
        bloblen = 6

        sigma = self.sigma()
        nblobs = len(self.kmer_charges)
        ans = 0

        for kmer in self.kmer_charges:

            # get the blob charge pattern list
            blob = kmer

            # calculate a bunch of parameters for the blob
            # with the blob sigma value being the ultimate
            # goal

            bpos = np.where(blob > 0)[0].size
            bneg = np.where(blob < 0)[0].size

            bncpr = (bpos - bneg) / (bloblen + 0.0)
            bfcr = (bpos + bneg) / (bloblen + 0.0)

            if(bfcr == 0):
                bsig = 0
            else:
                bsig = bncpr**2 / bfcr

            # calculate the square deviation of the
            # blob sigma from the sequence sigma and
            # weight by the number of blobs in the sequence
            ans += (sigma - bsig)**2 / nblobs

        return ans

    #...................................................................................#
    def sigma(self):
        """ Returns the sigma value for a sequence

           \sigma = \dfrac{NCPR^2}{FCR}

           When the sequence has one or more charged residues
           sigma = (NCPR^2)/FCR

           When the sequence has no charged residues
           sigma = 0
        """
        if(self.countNeut() == self.len):
            return 0
        else:
            return self.NCPR()**2 / self.FCR()

    # ...................................................................................#
    def FCR(self):
        return (self.countPos() + self.countNeg()) / (self.len + 0.0)

    # ...................................................................................#
    def NCPR(self):
        """ Get the net charge per residue of the sequence """
        return (self.countPos() - self.countNeg()) / (self.len + 0.0)

    #...................................................................................#
    def countPos(self):
        """ Get the number of positive residues in the sequence """
        return len(np.where(self.chargePattern > 0)[0])

    #...................................................................................#
    def countNeg(self):
        """ Get the number of negative residues in the sequence """
        return len(np.where(self.chargePattern < 0)[0])

    #...................................................................................#
    def countNeut(self):
        """ Get the number of neutral residues in the sequence """
        return len(np.where(self.chargePattern == 0)[0])


def main():
    pass


if __name__ == '__main__':
    main()