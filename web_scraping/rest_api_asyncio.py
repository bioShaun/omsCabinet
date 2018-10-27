#!/usr/bin/env python

import sys
import urllib
import urllib3
import requests
import time
import asyncio
import aiohttp


GENE_SERVERS_DICT = {
    'http://rest.ensembl.org': 'Ensembl',
    'https://rest.ensemblgenomes.org': 'EnsemblPlants'
}


def get_db(species):
    for each_server in GENE_SERVERS_DICT:
        ext = f"/info/assembly/{species}?"
        r = requests.get(
            each_server + ext, headers={"Content-Type": "application/json"})
        if r.ok:
            return GENE_SERVERS_DICT[each_server]
    sys.exit(f'{species} not found in Ensembl.')


class EnsemblRestClient(object):
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.server = server

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
    '''
    class to use EMBL-EBI API to download protein & go annotation

    >>> client = UniprotClient()

    test download go annotation
    >>> go_anno_obj = client.get_go_anno('GO:0008150')
    >>> go_anno_obj['results'][0]['id'] == 'GO:0008150'
    True
    >>> go_anno_obj['results'][0]['name'] == 'biological_process'
    True
    >>> go_anno_obj['results'][0]['aspect'] == 'biological_process'
    True
    '''

    def __init__(self, server='https://www.ebi.ac.uk', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    async def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Accept' not in hdrs:
            hdrs['Accept'] = 'application/json'

        if params:
            url = endpoint + '?' + urllib.parse.urlencode(params)
        else:
            url = endpoint

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                await asyncio.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        requests.adapters.DEFAULT_RETRIES = 5
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(self.server + url, headers=hdrs) as request:
                    if request.status == 200:
                        if hdrs['Accept'] == 'application/json':
                            data = await request.json()
                        elif 'text' in hdrs['Accept']:
                            data = await request.text()
                    elif request.status_code == 429:
                        # check if we are being rate limited by the server
                        if 'Retry-After' in request.headers:
                            retry = request.headers['Retry-After']
                            time.sleep(float(retry))
                            self.perform_rest_action(endpoint, hdrs, params)
                    else:
                        request.raise_for_status()
                    self.req_count += 1
        except (aiohttp.client_exceptions.ClientConnectorError,
                asyncio.TimeoutError,
                aiohttp.client_exceptions.ServerDisconnectedError,
                aiohttp.client_exceptions.InvalidURL) as exc:
            try:
                code = exc.code
            except AttributeError:
                code = ''
            sys.stderr.write(
                'Request failed for {0}: Message: {1} Status code: {2}.\n'.format(endpoint, exc, code))

        return data

    async def get_gene_inf(self, gene_id):
        gene_obj = await self.perform_rest_action(
            f'/proteins/api/proteins/{gene_id}',
            params={'size': -1}
        )
        return gene_obj

    async def get_protein_inf(self, uniprot_id):
        prot_obj = await self.perform_rest_action(
            '/proteins/api/proteins/{ui}'.format(ui=uniprot_id),
        )
        return prot_obj

    async def get_go_inf(self, uniprot_id, page=1):
        go_obj = await self.perform_rest_action(
            '/QuickGO/services/annotation/search',
            params={'geneProductType': 'protein',
                    'geneProductId': uniprot_id,
                    'page': page}
        )
        return go_obj

    async def get_go_anno(self, go_id):
        go_anno = await self.perform_rest_action(
            '/QuickGO/services/ontology/go/search',
            params={'query': go_id}
        )
        return go_anno


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
