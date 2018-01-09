import os
import requests, sys


class PullGo(object):

    def __init__(self, goid, fn):
        self.goid = goid
        self.fn = fn
        self.fno = '{}.tsv'.format(self.fn)
        self.fpo = os.path.join(os.path.dirname(__file__), '..', 'data', 'quickgo', self.fno)

    def query_quickgo(self):
        requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?" \
                     "goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in&" \
                     "goId={}&evidenceCode=ECO:0000269&evidenceCodeUsage=descendants&" \
                     "qualifier=part_of,colocalizes_with&geneProductType=protein".format(self.goid)
        r = requests.get(requestURL, headers={"Accept": "text/gpad"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        response_body = r.text
        self.write_response(response_body)
        return response_body

    def write_response(self, response_body):
        with open(self.fpo, 'w') as fo:
            fo.write(response_body)


def main():
    go_terms = {'Cajal_bodies': 'GO:0015030',
                'Centrosome': 'GO:0005813',
                'Nuclear_Speckles': 'GO:0016607',
                'Nucleolus': 'GO:0005730',
                'P_Body': 'GO:0000932',
                'PML_Body': 'GO:0016605',
                'Paraspeckle': 'GO:0042382',
                'Nuclear_Stress_Granule': 'GO:0097165',
                'Cytoplasmic_Stress_Granule': 'GO:0010494',
                'P_granule': 'GO:0043186'}
    for fn in go_terms:
        fpo = os.path.join(os.path.dirname(__file__), '..', 'data', 'quickgo', '{}.tsv'.format(fn))
        if not os.path.exists(fpo):
            goid = go_terms[fn]
            pg = PullGo(goid, fn)
            pg.query_quickgo()


if __name__ == '__main__':
    main()