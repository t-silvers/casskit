import sys
import json
import time
import requests
from urllib.parse import urlencode


class EnsemblRestClient(object):
    '''
    Adapted from https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client
    
    # Example usage
    species = "homo_sapiens"
    symbol = "BRCA2"
    run(species, symbol)
    '''
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            response = requests.get(self.server + endpoint, headers=hdrs)
            response.raise_for_status()
            content = response.content
            if content:
                data = json.loads(content)
            self.req_count += 1

        except requests.exceptions.HTTPError as e:
            # check if we are being rate limited by the server
            if e.response.status_code == 429:
                if 'Retry-After' in e.response.headers:
                    retry = e.response.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.response.status_code} Reason: {1.response.reason}\n'.format(endpoint, e))

        return data

    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            endpoint='/xrefs/symbol/{0}/{1}'.format(species, symbol),
            params={'object_type': 'gene'}
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants
        return None

def run(species, symbol):
    client = EnsemblRestClient()
    variants = client.get_variants(species, symbol)
    if variants:
        for v in variants:
            print('{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})'.format(**v))