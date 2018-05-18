import click
import requests
import os
import time
from bs4 import BeautifulSoup
import pandas as pd


ENSEMBL_SEARCH_URL = 'https://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q={s_id};site=ensemblunit'

SERVER = "https://rest.ensemblgenomes.org"
EXT = "/sequence/id/{ensembl_id}?type=genomic"


def get_ensembl_seq(ensembl_id):
    id_ext = EXT.format(ensembl_id=ensembl_id)
    url = '{sev}{ext}'.format(sev=SERVER,
                              ext=id_ext)
    r = requests.get(url, headers={"Content-Type": "text/x-fasta"})
    if not r.ok:
        return None
    else:
        return r.text


def get_ensembl_id(query_id):
    each_url = ENSEMBL_SEARCH_URL.format(s_id=query_id)
    try:
        r = requests.get(each_url)
    except Exception:
        return 'connection_failed'
    time.sleep(10)
    soup = BeautifulSoup(r.content, 'lxml')
    search_result = soup.find('div', {'class': 'searchresults'})
    ensembl_id = None
    if search_result:
        ensembl_id = search_result.a.text
    return ensembl_id


@click.command()
@click.argument(
    'id_list',
    type=click.File('r'),
    required=True
)
@click.argument(
    'out_dir',
    type=click.Path(file_okay=False),
    default=os.getcwd()
)
@click.option(
    '-m',
    '--name_map',
    type=click.Path(dir_okay=False),
    default=None
)
def main(id_list, out_dir, name_map):
    # load downloaded info
    id_map_file = os.path.join(out_dir, 'id.map.txt')
    if os.path.exists(id_map_file) and os.stat(id_map_file).st_size:
        id_map_df = pd.read_table(id_map_file, index_col=0,
                                  header=None, names=['ensembl_id'])
        id_map_df = id_map_df.drop_duplicates()
    else:
        id_map_df = pd.DataFrame({'ensembl_id': []})
    # load failed info
    failed_id_file = os.path.join(out_dir, 'failed.id.txt')
    if os.path.exists(failed_id_file) and os.stat(failed_id_file).st_size:
        failed_id_df = pd.read_table(
            failed_id_file, header=None, names=['failed_id'])
    else:
        failed_id_df = pd.DataFrame({'failed_id': []})
    # get prepared name_map df to accelarate
    if name_map is None:
        name_map_df = pd.DataFrame({'ensembl_id': []})
    else:
        name_map_df = pd.read_table(name_map, index_col=0)
        name_map_df.dropna(inplace=True)
    # make sure out_dir exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # prepare for output file
    genomic_seq_file = os.path.join(out_dir, 'ensembl.genomic.seq.fa')
    genomic_seq_inf = open(genomic_seq_file, 'a')
    # read input ids
    id_list = [each.strip() for each in id_list]
    for each_id in id_list:
        # skip downloaded ids
        if each_id in id_map_df.index:
            continue
        name_map_success = 1
        if each_id in name_map_df.index:
            ensembl_id = str(name_map_df.loc[each_id, 'ensembl_id'])
            ensembl_id = ensembl_id.split()[-1]
            ensembl_seq = get_ensembl_seq(ensembl_id)
            if ensembl_seq is None:
                name_map_success = 0
        else:
            name_map_success = 0
        if not name_map_success:
            ensembl_id = get_ensembl_id(each_id)
            if ensembl_id:
                if ensembl_id == 'connection_failed':
                    ensembl_seq = None
                else:
                    ensembl_seq = get_ensembl_seq(ensembl_id)
                if ensembl_seq is None:
                    failed_id_df.loc[len(failed_id_df),
                                     'failed_id'] = ensembl_id
            else:
                ensembl_seq = None
                ensembl_id = 'Not_found'
        if ensembl_seq and (ensembl_id not in id_map_df.ensembl_id):
            genomic_seq_inf.write(ensembl_seq)
        id_map_df.loc[each_id, 'ensembl_id'] = ensembl_id
    genomic_seq_inf.close()
    id_map_df.to_csv(id_map_file, sep='\t', header=False)
    failed_id_df.to_csv(failed_id_file, sep='\t', header=False,
                        index=False)


if __name__ == '__main__':
    main()
