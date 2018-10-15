import fire
from rest_api import EnsemblRestClient, UniprotClient
import pandas as pd
from pandas import DataFrame
import gtfparse
import time
from tqdm import tqdm
from pathlib import Path
from functools import reduce
import glob
import sys


ENS_SERVER = "https://rest.ensembl.org"
OUT_HEADER_BASE = [
    'gene_id',
    'gene_names',
    'uniprot_id',
    'protein_names',
    'protein_existence']

OUT_HEADER_REF = [
    'interpro_ids',
    'interpro_names',
    'pfam_ids',
    'pfam_names',
    'feature_names',
    'feature_dbs',
    'pubmed_ids'
]


def save_download(df, middle_file):
    middle_file = Path(middle_file)
    if middle_file.exists():
        write_mode = 'a'
        write_header = False
    else:
        write_header = True
        write_mode = 'w'
    df.to_csv(middle_file, header=write_header,
              mode=write_mode, index=False,
              sep='\t')


def save_download_to_dir(df, middle_file_dir, file_id):
    middle_file_dir = Path(middle_file_dir)
    middle_file_dir.mkdir(parents=True, exist_ok=True)
    middle_file = middle_file_dir / f'{file_id}.txt'
    df.to_csv(middle_file, index=False, sep='\t')


def ens_unip_map(ensembl_id, middle_file):
    client = EnsemblRestClient(server=ENS_SERVER)
    decoded = client.get_uniprot_id(ensembl_id)
    uniprot_ids = [each['primary_id']
                   for each in decoded if 'primary_id' in each]
    if uniprot_ids:
        unprot_df = DataFrame(uniprot_ids, columns=['uniprot_id'])
    else:
        unprot_df = DataFrame([None], columns=['uniprot_id'])
    unprot_df.loc[:, 'gene_id'] = ensembl_id
    save_download(unprot_df, middle_file)
    return unprot_df


def uniprot_go_map(uniprot_id, middle_file):

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
    else:
        go_df = DataFrame([None], columns=['go_id'])
    go_df.loc[:, 'uniprot_id'] = uniprot_id
    save_download(go_df, middle_file)
    return go_df


def idlist2df(id_list, col_name, identity_map=None):
    if id_list:
        list_df = DataFrame(id_list, columns=[col_name])
        if identity_map is not None:
            for key, val in identity_map.items():
                list_df.loc[:, key] = val
        return list_df


def refdb_anno(anno_db, db_name):
    db_ids = [each['id'] for each in anno_db
              if each['type'] == db_name]
    db_names = [each['properties']['entry name'] for each in anno_db
                if each['type'] == db_name]
    return db_ids, db_names


def featuredb_anno(featuredb):
    no_use_ft = ('', 'DSL')
    feature_db_dict = dict()
    for each in featuredb:
        if 'description' not in each:
            continue
        elif each['description'] in no_use_ft:
            continue
        elif 'evidences' not in each:
            continue
        elif 'source' not in each['evidences'][0]:
            continue
        else:
            feature_db_name = each['evidences'][0]['source']['name']
            feature_db_id = each['evidences'][0]['source']['id']
            if feature_db_name == 'Pfam':
                continue
            feature_db_dict.setdefault(
                'feature_names', []).append(each['description'])
            feature_db_f_name = f'{feature_db_name}:{feature_db_id}'
            feature_db_dict.setdefault(
                'feature_dbs', []
            ).append(feature_db_f_name)
    return feature_db_dict


def commentdb_anno(commentdb):
    comment_db_dict = dict()
    for each in commentdb:
        each_type = f'uniprot_comments({each["type"]})'
        if 'text' not in each:
            continue
        for each_type_cm in each['text']:
            if 'value' in each_type_cm:
                comment_db_dict.setdefault(each_type, []).append(
                    each_type_cm['value']
                )
    return comment_db_dict


def uniprot_anno_map(uniprot_id, middle_file):
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
    if 'proteinExistence' in decoded:
        anno_dfs.append(idlist2df([decoded['proteinExistence']],
                                  'protein_existence',
                                  identity_map=identity_map))
    if 'comments' in decoded:
        comment_db = decoded['comments']
        comment_db_dict = commentdb_anno(comment_db)
        if comment_db_dict:
            for key, val in comment_db_dict.items():
                anno_dfs.append(idlist2df(val, key,
                                          identity_map=identity_map))
    if 'dbReferences' in decoded:
        anno_db = decoded['dbReferences']
        interpro_ids, interpro_names = refdb_anno(anno_db, 'InterPro')
        if interpro_ids:
            anno_dfs.append(idlist2df(interpro_ids, 'interpro_ids',
                                      identity_map=identity_map))
            anno_dfs.append(idlist2df(interpro_names, 'interpro_names',
                                      identity_map=identity_map))
        pfam_ids, pfam_names = refdb_anno(anno_db, 'Pfam')
        if pfam_ids:
            anno_dfs.append(idlist2df(pfam_ids, 'pfam_ids',
                                      identity_map=identity_map))
            anno_dfs.append(idlist2df(pfam_names, 'pfam_names',
                                      identity_map=identity_map))
    if 'features' in decoded:
        feature_db = decoded['features']
        feature_db_dict = featuredb_anno(feature_db)
        if feature_db_dict:
            for key, val in feature_db_dict.items():
                anno_dfs.append(idlist2df(val, key,
                                          identity_map=identity_map))
    if 'references' in decoded:
        citation_db = decoded['references']
        pubmed_ids = [each['citation']['dbReferences'][0]['id']
                      for each in citation_db
                      if 'dbReferences' in each['citation']]
        if pubmed_ids:
            anno_dfs.append(idlist2df(pubmed_ids, 'pubmed_ids',
                                      identity_map=identity_map))
    if anno_dfs:
        anno_dfs = [each for each in anno_dfs
                    if each is not None]
        anno_df = reduce(pd.merge, anno_dfs)
    else:
        anno_df = DataFrame([None], columns=['protein_existence'])
        anno_df.loc[:, 'uniprot_id'] = uniprot_id
    save_download_to_dir(anno_df, middle_file, uniprot_id)
    return anno_df


def map2df(map_func, query_ids, msg,
           middle_file, id_col_name='uniprot_id'):
    if middle_file.exists():
        if middle_file.is_file():
            finished_df = pd.read_table(middle_file)
        elif middle_file.is_dir():
            middle_files = glob.glob(f'{middle_file}/*txt')
            finished_dfs = [pd.read_table(each_file) for each_file
                            in middle_files]
            finished_df = pd.concat(finished_dfs, sort=False)
        else:
            sys.exit(f'unsupported file type of {middle_file}!')
        left_ids = set(query_ids).difference(
            set(finished_df.loc[:, id_col_name]))
    else:
        finished_df = DataFrame([])
        left_ids = query_ids

    map_dfs = list()
    for each_id in tqdm(left_ids,
                        ncols=100,
                        desc=f'{msg:<40}'):
        each_df = map_func(each_id, middle_file)
        if each_df is not None:
            map_dfs.append(each_df)
    if map_dfs:
        map_df = pd.concat(map_dfs, sort=False)
        map_df = pd.concat([finished_df, map_df], sort=False)
    else:
        map_df = finished_df
    return map_df


def format_df(gene_df, sep=',', empty_rep='--'):

    def my_unique(x, sep=','):
        unique_x = list(pd.unique(x.dropna()))
        if not unique_x:
            unique_x = ['--']
        return sep.join(unique_x)

    return gene_df.groupby('gene_id').agg(my_unique, sep)


def swissprot_id_from_fa(fa_id):
    return fa_id.split('|')[1]


def go_annotation(input_file):
    input_file = Path(input_file)
    if input_file.suffix == '.gtf':
        # ensembl genes from gtf
        gtf_df = gtfparse.read_gtf(input_file)
        gtf_genes = set(gtf_df.gene_id.unique())
        ens_uni_map_file = input_file.with_suffix('.uni_id.map')
        map_gene_id_msg = 'Mapping ensembl id <-> uniprot id'
        ens_uni_map_df = map2df(ens_unip_map, gtf_genes,
                                map_gene_id_msg, ens_uni_map_file,
                                id_col_name='gene_id')
        ens_uni_map_df.dropna(inplace=True)
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
        uni_go_map_file = input_file.with_suffix('.uni_go.map')
        uni_go_map_df = map2df(uniprot_go_map,
                               ens_uni_map_df.uniprot_id.unique(),
                               map_uni_to_go_msg,
                               uni_go_map_file)
        uni_go_map_df.dropna(inplace=True)
        ens_go_df = ens_uni_map_df.merge(
            uni_go_map_df).drop('uniprot_id', axis=1)
        # format go result
        formated_go = format_df(ens_go_df)
        formated_go.to_csv(go_file, sep='\t',
                           header=False)
    # uniprot to annotation
    anno_file = input_file.with_suffix('.anno.txt')
    # if not anno_file.exists():
    if True:
        uni_anno_map_msg = 'Retriving UniProt annotation'
        uni_anno_map_file = input_file.with_suffix('.uni_anno.map')
        uni_anno_map_df = map2df(uniprot_anno_map,
                                 ens_uni_map_df.uniprot_id.unique(),
                                 uni_anno_map_msg,
                                 uni_anno_map_file)
        if uni_anno_map_df is not None:
            ens_anno_df = ens_uni_map_df.merge(
                uni_anno_map_df)
            formated_anno = format_df(ens_anno_df, sep='|')
            uniprot_comment_header = [each for each in formated_anno.columns
                                      if each.startswith('uniprot_comments')]
            out_header = OUT_HEADER_BASE + uniprot_comment_header + OUT_HEADER_REF
            out_header = [each for each in out_header
                          if each in formated_anno.columns]
            formated_anno.to_csv(anno_file, sep='\t',
                                 columns=out_header)
        else:
            print('No annotation found.')


if __name__ == '__main__':
    fire.Fire(go_annotation)
