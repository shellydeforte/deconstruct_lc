import requests


def fetch_uniprot(uniID):
    """
    status codes:
    200 -- The request was processed successfully
    400 -- Bad request. There was a problem with your input
    404 -- Not found. The resource you requested doesn't exist
    410 -- Gone. The resource you requested was removed
    500 -- Internal server error. Most likely a temporary problem
    503 -- Service not available. The server is being updated
    """
    headers = {'From': 'shelly.deforte@gmail.com'}
    try:
        r = requests.get(
            'http://www.uniprot.org/uniprot/{}.fasta'.format(uniID),
            timeout=5, headers=headers)
        if r.status_code == 200:
            return r
    except:
        return None


def write_fasta(uniIDs, fp):
    with open(fp, 'w') as fpo:
        for uniID in uniIDs:
            fasta = fetch_uniprot(uniID)
            if fasta:
                fpo.write(fasta.text)
            else:
                raise Exception("Error pulling UniProt entry")