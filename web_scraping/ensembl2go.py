import fire
from rest_api import EnsemblRestClient, UniprotClient
import pandas as pd
from pandas import DataFrame
import gtfparse
import time
from tqdm import tqdm
from pathlib import Path
from functools import reduce


ENS_SERVER = "https://rest.ensembl.org"


def ens_unip_map(ensembl_id):
    client = EnsemblRestClient(server=ENS_SERVER)
    decoded = client.get_uniprot_id(ensembl_id)
    uniprot_ids = [each['primary_id']
                   for each in decoded if 'primary_id' in each]
    if uniprot_ids:
        unprot_df = DataFrame(uniprot_ids, columns=['uniprot_id'])
        unprot_df.loc[:, 'gene_id'] = ensembl_id
        return unprot_df
    else:
        return None


def uniprot_go_map(uniprot_id):

    decoded_list = []

    def download_page_inf(uniprot_id, page=1):
        client = UniprotClient()
        nonlocal decoded_list
        decoded = client.get_go_inf(uniprot_id, page=page)
        decoded_list.append(decoded)
        total_page = decoded['pageInfo']['total']
        cur_page = decoded['pageInfo']['current']
        if cur_page < total_page:
            return download_page_inf(uniprot_id, cur_page + 1)

    download_page_inf(uniprot_id)
    total_goids = []
    for decoded in decoded_list:
        if 'results' in decoded:
            go_ids = [each['goId'] for each in decoded['results']
                      if 'goId' in each]
            total_goids.extend(go_ids)
    if total_goids:
        go_df = DataFrame(total_goids, columns=['go_id'])
        go_df.loc[:, 'uniprot_id'] = uniprot_id
        return go_df
    else:
        return None


def idlist2df(id_list, col_name, identity_map=None):
    if id_list:
        list_df = DataFrame(id_list, columns=[col_name])
        if identity_map is not None:
            for key, val in identity_map.items():
                list_df.loc[:, key] = val
        return list_df


def uniprot_anno_map(uniprot_id):
    client = UniprotClient()
    decoded = client.get_protein_inf(uniprot_id)
    anno_dfs = list()
    identity_map = {'uniprot_id': uniprot_id}
    if 'gene' in decoded:
        gene_db = decoded['gene']
        gene_names = [each['name']['value'] for each in gene_db
                      if 'name' in each]
        anno_dfs.append(idlist2df(gene_names, 'gene_names',
                                  identity_map=identity_map))
    if 'protein' in decoded:
        protein_db = decoded['protein']
        products = None
        if 'submittedName' in protein_db:
            products = [each['fullName']['value'] for
                        each in protein_db['submittedName']
                        if 'value' in each['fullName']]
        elif 'recommendedName' in protein_db:
            products = [protein_db['recommendedName']['fullName']['value']]
        if products is not None:
            anno_dfs.append(idlist2df(products, 'protein_names',
                                      identity_map=identity_map))
    if 'dbReferences' in decoded:
        anno_db = decoded['dbReferences']
        interpro_ids = [each['id'] for each in anno_db
                        if each['type'] == 'InterPro']
        anno_dfs.append(idlist2df(interpro_ids, 'interpro_id',
                                  identity_map=identity_map))
    if 'references' in decoded:
        citation_db = decoded['references']
        pubmed_ids = [each['citation']['dbReferences'][0]['id']
                      for each in citation_db
                      if 'dbReferences' in each['citation']]
        anno_dfs.append(idlist2df(pubmed_ids, 'pubmed_id',
                                  identity_map=identity_map))
    if anno_dfs:
        anno_dfs = [each for each in anno_dfs
                    if each is not None]
        anno_df = reduce(pd.merge, anno_dfs)
        return anno_df


def map2df(map_func, query_ids, msg):
    map_dfs = map(map_func,
                  tqdm(query_ids,
                       ncols=100,
                       desc=f'{msg:<40}'))
    map_dfs = [each for each in map_dfs
               if each is not None]
    if map_dfs:
        map_df = pd.concat(map_dfs)
        return map_df


def format_df(gene_df, sep=','):

    def my_unique(x, sep=','):
        unique_x = list(pd.unique(x))
        return sep.join(unique_x)

    return gene_df.groupby('gene_id').agg(my_unique, sep)


def swissprot_id_from_fa(fa_id):
    return fa_id.split('|')[1]


def go_annotation(input_file):
    input_file = Path(input_file)
    if input_file.suffix == '.gtf':
        # ensembl genes from gtf
        gtf_df = gtfparse.read_gtf(input_file)
        gtf_genes = gtf_df.gene_id.unique()
        map_gene_id_msg = 'Mapping ensembl id <-> uniprot id'
        ens_uni_map_df = map2df(ens_unip_map, gtf_genes, map_gene_id_msg)
    elif input_file.suffix == '.blasttab':
        # uniprot blast out
        blast_df = pd.read_table(input_file, header=None)
        ens_uni_map_df = blast_df.loc[:, [0, 1]].drop_duplicates()
        ens_uni_map_df.columns = ['gene_id', 'uniprot_id']
        ens_uni_map_df.loc[:, 'uniprot_id'] = list(ens_uni_map_df.uniprot_id.map(
            swissprot_id_from_fa))
    # uniprot to go
    go_file = input_file.with_suffix('.go.txt')
    if not go_file.exists():
        map_uni_to_go_msg = 'Retriving GO ids'
        uni_go_map_df = map2df(uniprot_go_map,
                               ens_uni_map_df.uniprot_id.unique(),
                               map_uni_to_go_msg)
        ens_go_df = ens_uni_map_df.merge(
            uni_go_map_df).drop('uniprot_id', axis=1)
        # format go result
        formated_go = format_df(ens_go_df)
        formated_go.to_csv(go_file, sep='\t',
                           header=False)
    # uniprot to annotation
    anno_file = input_file.with_suffix('.anno.txt')
    if not anno_file.exists():
        uni_anno_map_msg = 'Retriving UniProt annotation'
        uni_anno_map_df = map2df(uniprot_anno_map,
                                 ens_uni_map_df.uniprot_id.unique(),
                                 uni_anno_map_msg)
        if uni_anno_map_df is not None:
            ens_anno_df = ens_uni_map_df.merge(
                uni_anno_map_df)
            ens_anno_df.fillna('--', inplace=True)
            formated_anno = format_df(ens_anno_df, sep='|')
            formated_anno.to_csv(anno_file, sep='\t')
        else:
            print('No annotation found.')


if __name__ == '__main__':
    fire.Fire(go_annotation)
