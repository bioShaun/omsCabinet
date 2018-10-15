#!/usr/bin/env python

import sys
import urllib
import requests
import time


class EnsemblRestClient(object):
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
            endpoint += '?' + urllib.parse.urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        requests.adapters.DEFAULT_RETRIES = 5
        try:
            request = requests.get(self.server + endpoint, headers=hdrs)
            if request.ok:
                if hdrs['Content-Type'] == 'application/json':
                    data = request.json()
                elif 'text' in hdrs['Content-Type']:
                    data = request.text
            else:
                request.raise_for_status()
            self.req_count += 1

        except requests.exceptions.HTTPError as e:
            # check if we are being rate limited by the server
            if request.status_code == 429:
                if 'Retry-After' in request.headers:
                    retry = request.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            '/xrefs/symbol/{0}/{1}'.format(species, symbol),
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

    def get_orthologues(self, ens_id, taxon=4565):
        othos = self.perform_rest_action(
            '/homology/id/{ei}'.format(ei=ens_id),
            params={'type': 'orthologues',
                    'target_taxon': taxon}
        )
        return othos

    def get_seq(self, ens_id):
        seq_obj = self.perform_rest_action(
            '/sequence/id/{ei}'.format(ei=ens_id),
            hdrs={'Content-Type': "text/x-fasta"},
            params={'type': 'genomic'}
        )
        return seq_obj

    def get_uniprot_id(self, ens_id):
        uniprot_obj = self.perform_rest_action(
            '/xrefs/id/{ei}'.format(ei=ens_id),
            params={'external_db': 'Uniprot_gn'}
        )
        return uniprot_obj


class UniprotClient(object):
    def __init__(self, server='https://www.ebi.ac.uk', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Accept' not in hdrs:
            hdrs['Accept'] = 'application/json'

        if params:
            endpoint += '?' + urllib.parse.urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        requests.adapters.DEFAULT_RETRIES = 5
        try:
            request = requests.get(self.server + endpoint, headers=hdrs)
            if request.ok:
                if hdrs['Accept'] == 'application/json':
                    data = request.json()
                elif 'text' in hdrs['Accept']:
                    data = request.text
            else:
                request.raise_for_status()
            self.req_count += 1

        except requests.exceptions.HTTPError:
            # check if we are being rate limited by the server
            if request.status_code == 429:
                if 'Retry-After' in request.headers:
                    retry = request.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_protein_inf(self, uniprot_id):
        prot_obj = self.perform_rest_action(
            '/proteins/api/proteins/{ui}'.format(ui=uniprot_id),
        )
        return prot_obj

    def get_go_inf(self, uniprot_id, page=1):
        go_obj = self.perform_rest_action(
            '/QuickGO/services/annotation/search',
            params={'geneProductType': 'protein',
                    'geneProductId': uniprot_id,
                    'page': page}
        )
        return go_obj


def run(species, symbol):
    client = EnsemblRestClient()
    variants = client.get_variants(species, symbol)
    if variants:
        for v in variants:
            print(
                '{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})'.format(**v))


if __name__ == '__main__':
    if len(sys.argv) == 3:
        species, symbol = sys.argv[1:]
    else:
        species, symbol = 'human', 'BRAF'

    run(species, symbol)
