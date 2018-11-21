import fire
from rest_api_asyncio import UniprotClient, get_db
import pandas as pd
from pandas import DataFrame
import gtfparse
import time
from tqdm import tqdm
from pathlib import Path
from functools import reduce
import glob
import sys
import urllib3
import asyncio


OUT_HEADER_BASE = [
    'gene_id',
    'gene_names',
    'uniprot_id',
    'protein_names',
    'protein_existence']

OUT_HEADER_REF = [
    'go_id',
    'go_term',
    'interpro_ids',
    'interpro_names',
    'pfam_ids',
    'pfam_names',
    'feature_names',
    'feature_dbs',
    'pubmed_ids'
]


DEFAULT_CONCUR_REQ = 10
DEFAULT_RETRY_TIMES = 2


def format_df(gene_df, sep=',', empty_rep=None, by='gene_id'):

    def my_unique(x, sep=','):
        unique_x = pd.unique(x.dropna())
        if str(unique_x.dtype) == 'float64':
            unique_x = unique_x.astype('int')
        unique_x = [str(each) for each in unique_x]
        if not unique_x:
            if empty_rep is None:
                return None
            else:
                unique_x = [empty_rep]
        return sep.join(unique_x)

    gene_df = gene_df.groupby(by).agg(my_unique, sep)
    return gene_df.reset_index()


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


def extract_anno_inf(decoded, uniprot_id, middle_file):
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
    return anno_df


async def go_anno_map(go_id, middle_file, client, semaphore, skip):
    with await semaphore:
        decoded = await client.get_go_anno(go_id)
    go_df = DataFrame([go_id, None, None],
                      index=['go_id', 'go_term', 'go_ontology']).T
    if decoded is None:
        if not skip:
            return None
    for each_result in decoded['results']:
        if each_result['id'] == go_id:
            go_df = DataFrame([go_id,
                               each_result['name'],
                               each_result['aspect']],
                              index=['go_id',
                                     'go_term',
                                     'go_ontology']).T
    save_download(go_df, middle_file)
    return go_df


async def uniprot_go_map(uniprot_id, middle_file, client,
                         semaphore, skip=False):
    '''
    function to download go ids using uniport id
    '''
    decoded_list = []

    async def download_page_inf(uniprot_id, semaphore, page=1):
        nonlocal client
        nonlocal decoded_list
        with await semaphore:
            decoded = await client.get_go_inf(uniprot_id, page=page)
        decoded_list.append(decoded)
        if decoded is None:
            return None
        total_page = decoded['pageInfo']['total']
        cur_page = decoded['pageInfo']['current']
        if cur_page < total_page:
            return await download_page_inf(uniprot_id, semaphore, cur_page + 1)

    await download_page_inf(uniprot_id, semaphore)
    total_goids = []
    for decoded in decoded_list:
        if decoded is None:
            if not skip:
                return None
            else:
                continue
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


async def uniprot_anno_map(uniprot_id, middle_file, client,
                           semaphore, skip=False):
    with await semaphore:
        decoded = await client.get_protein_inf(uniprot_id)
    if decoded is None:
        if not skip:
            return None
        else:
            decoded = {}
    anno_df = extract_anno_inf(decoded, uniprot_id, middle_file)
    anno_df = format_df(anno_df, sep='|', by='uniprot_id')
    save_download_to_dir(anno_df, middle_file, uniprot_id)
    return anno_df


async def ens_anno_map(ensembl_id, middle_file, client,
                       semaphore, skip=False):
    with await semaphore:
        decodeds = await client.get_gene_inf(ensembl_id)
    anno_dfs = list()
    if decodeds is not None:
        for decoded in decodeds:
            uniprot_id = decoded['accession']
            anno_dfs.append(
                extract_anno_inf(decoded, uniprot_id,
                                 middle_file))
    else:
        if not skip:
            return None
    if anno_dfs:
        anno_df = pd.concat(anno_dfs, sort=False)
        anno_df = format_df(anno_df, sep='|', by='uniprot_id')
    else:
        anno_df = DataFrame([None], columns=['uniprot_id'])
    anno_df.loc[:, 'gene_id'] = ensembl_id
    out_ensembl_id = ensembl_id.split(':')[-1]
    save_download_to_dir(anno_df, middle_file, out_ensembl_id)
    return anno_df


async def map2df(map_func, query_ids, msg,
                 middle_file, concur_req,
                 retry_limits,
                 id_col_name='uniprot_id',
                 skip=False):
    client = UniprotClient()
    semaphore = asyncio.Semaphore(concur_req)
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
    download_fail = 0
    to_do = [map_func(each_id, middle_file, client, semaphore, skip)
             for each_id in left_ids]
    to_do_iter = asyncio.as_completed(to_do)
    to_do_iter = tqdm(to_do_iter, total=len(left_ids),
                      ncols=100, desc=f'{msg:<40}')
    for future in to_do_iter:
        try:
            each_df = await future
        except FetchError as exc:
            country_code = exc.country_code
            try:
                error_msg = exc.__cause__.args[0]
            except IndexError:
                error_msg = exc.__cause__.__class__.__name__
                msg = '*** Error for {}: {}'
                print(msg.format(country_code, error_msg))
                download_fail += 1
        else:
            if each_df is None:
                download_fail += 1
            else:
                map_dfs.append(each_df)
    if download_fail > 0 and retry_limits > 0:
        print(f'{download_fail} items failed to download.')
        print(f'last {retry_limits} try!')
        retry_limits = retry_limits - 1
        if retry_limits == 1:
            skip = True
        return await map2df(map_func, query_ids, msg,
                            middle_file, concur_req,
                            retry_limits,
                            id_col_name=id_col_name,
                            skip=skip)
    if map_dfs:
        map_df = pd.concat(map_dfs, sort=False)
        map_df = pd.concat([finished_df, map_df], sort=False)
    else:
        map_df = finished_df
    return map_df


def map2df_summary(map_func, query_ids, msg,
                   middle_file, concur_req,
                   retry_limits, loop,
                   id_col_name='uniprot_id'):
    coro = map2df(map_func, query_ids, msg,
                  middle_file, concur_req,
                  retry_limits,
                  id_col_name=id_col_name)
    df = loop.run_until_complete(coro)
    return df


def swissprot_id_from_fa(fa_id):
    return fa_id.split('|')[1]


def rm_db_name(ens_id):
    return ens_id.split(':')[-1]


def format_uniprot_anno(ens_anno_df, anno_file):
    uniprot_comment_header = [each for each in ens_anno_df.columns
                              if each.startswith('uniprot_comments')]
    out_header = OUT_HEADER_BASE + uniprot_comment_header + OUT_HEADER_REF
    out_header = [each for each in out_header
                  if each in ens_anno_df.columns]
    ens_anno_df.to_csv(anno_file, sep='\t',
                       columns=out_header,
                       na_rep='--', index=False)


def go_annotation(input_file, species=None, workers=DEFAULT_CONCUR_REQ,
                  retry=DEFAULT_RETRY_TIMES):
    input_file = Path(input_file)
    anno_file = input_file.with_suffix('.anno.txt')
    gene_list_file = input_file.with_suffix('.pcg.gene.list')
    ens_uni_map_file = input_file.with_suffix('.uni_id.map')
    ens_anno_map_df = None
    loop = asyncio.get_event_loop()

    if input_file.suffix == '.gtf':
        # need species info to choose ensembl sever
        if species is None:
            sys.exit('species info is needed to choose ensembl server.')
        else:
            ens_db = get_db(species)
        # ensembl genes from gtf
        if gene_list_file.is_file():
            gene_list_df = pd.read_table(gene_list_file, header=None,
                                         names=['gene_id'])
            gtf_genes = gene_list_df.gene_id.unique()
        else:
            gtf_df = gtfparse.read_gtf(input_file)
            if 'gene_biotype' in gtf_df.columns:
                gtf_df = gtf_df[gtf_df.gene_biotype == 'protein_coding']
            gtf_df = gtf_df.drop_duplicates(subset=['gene_id'])
            gtf_df.to_csv(gene_list_file, sep='\t', index=False,
                          columns=['gene_id'], header=False)
            gtf_genes = set(gtf_df.gene_id.unique())
        gtf_genes = [f'{ens_db}:{each}' for each in gtf_genes]
        # ensembl to uniprot annotation
        map_gene_id_msg = 'Mapping ensembl id <-> uniprot annotation'
        ens_anno_df = map2df_summary(ens_anno_map, gtf_genes,
                                     map_gene_id_msg, ens_uni_map_file,
                                     workers, retry, loop,
                                     id_col_name='gene_id')
        ens_anno_df.loc[:, 'gene_id'] = ens_anno_df.gene_id.map(rm_db_name)
        ens_uni_map_df = ens_anno_df.loc[
            :, ['gene_id', 'uniprot_id']].dropna()
    elif input_file.suffix == '.blasttab':
        # uniprot blast out
        blast_df = pd.read_table(input_file, header=None)
        ens_uni_map_df = blast_df.loc[:, [0, 1]].drop_duplicates()
        ens_uni_map_df.columns = ['gene_id', 'uniprot_id']
        ens_uni_map_df.loc[:, 'uniprot_id'] = list(ens_uni_map_df.uniprot_id.map(
            swissprot_id_from_fa))
        # uniprot to annotation
        uni_anno_map_msg = 'Retriving UniProt annotation'
        uni_anno_map_file = input_file.with_suffix('.uni_anno.map')
        uni_anno_map_df = map2df_summary(uniprot_anno_map,
                                         ens_uni_map_df.uniprot_id.unique(),
                                         uni_anno_map_msg,
                                         ens_uni_map_file,
                                         workers, retry,
                                         loop)
        ens_anno_df = ens_uni_map_df.merge(
            uni_anno_map_df)
    # uniprot to go
    go_file = input_file.with_suffix('.go.txt')
    map_uni_to_go_msg = 'Retriving GO ids'
    uni_go_map_file = input_file.with_suffix('.uni_go.map')
    uni_go_map_df = map2df_summary(uniprot_go_map,
                                   ens_uni_map_df.uniprot_id.unique(),
                                   map_uni_to_go_msg,
                                   uni_go_map_file,
                                   workers, retry, loop)
    if not go_file.exists():
        uni_go_map_df.dropna(inplace=True)
        ens_go_df = ens_uni_map_df.merge(
            uni_go_map_df).drop('uniprot_id', axis=1)
        # format go result
        formated_go = format_df(ens_go_df)
        formated_go.to_csv(go_file, sep='\t',
                           header=False, index=False)
    go_anno_map_file = input_file.with_suffix('.go_anno.map')
    map_go_to_anno_msg = 'Retriving GO information'
    go_anno_map_df = map2df_summary(go_anno_map,
                                    uni_go_map_df.go_id.unique(),
                                    map_go_to_anno_msg,
                                    go_anno_map_file,
                                    workers, retry, loop,
                                    id_col_name='go_id')
    go_anno_map_df = go_anno_map_df.merge(uni_go_map_df)
    ens_anno_df = ens_anno_df.merge(go_anno_map_df)
    ens_anno_df = format_df(ens_anno_df, empty_rep='--', sep='|')
    format_uniprot_anno(ens_anno_df, anno_file)
    loop.close()


if __name__ == '__main__':
    fire.Fire(go_annotation)
