#import sys
import click
import os
import pandas as pd
import numpy as np
from collections import OrderedDict
from Bio import SeqIO
from rest_api import EnsemblRestClient


SERVER = "https://rest.ensemblgenomes.org"


def get_ensembl_seq(ensembl_id):
    client = EnsemblRestClient(server=SERVER)
    return client.get_seq(ensembl_id)


def homology_json_to_df(ensembl_id, homo_inf):
    homo_dict = OrderedDict()
    if homo_inf:
        for each_homo in homo_inf:
            homo_dict.setdefault('Ensembl_id', []).append(ensembl_id)
            homo_dict.setdefault('Orthologue', []).append(
                each_homo['target']['id'])
            homo_dict.setdefault('Target_percent', []).append(
                each_homo['target']['perc_id'])
            homo_dict.setdefault('Query_percent', []).append(
                each_homo['source']['perc_id'])
    else:
        homo_dict.setdefault('Ensembl_id', []).append(ensembl_id)
        homo_dict.setdefault('Orthologue', []).append(np.nan)
        homo_dict.setdefault('Target_percent', []).append(0)
        homo_dict.setdefault('Query_percent', []).append(0)
    homo_df = pd.DataFrame(homo_dict)
    homo_df = homo_df.set_index('Ensembl_id')
    return homo_df


def get_ensembl_orthologues(ensembl_id):
    client = EnsemblRestClient(server=SERVER)
    decoded = client.get_orthologues(ensembl_id)
    try:
        homo_inf = decoded['data'][0]['homologies']
    except IndexError:
        homo_inf = None
    return homology_json_to_df(ensembl_id, homo_inf)


@click.command()
@click.argument(
    'id_map',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'out_dir',
    type=click.Path(file_okay=False),
    default=os.getcwd()
)
def main(id_map, out_dir):
    # check out_dir existence
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load downloaded info
    meta_file = os.path.join(out_dir, 'meta_table.txt')
    download_table = os.path.join(out_dir, 'othologues_table.txt')
    if os.path.exists(meta_file) and os.stat(meta_file).st_size:
        meta_df = pd.read_table(meta_file, index_col=0)
    else:
        meta_df = pd.DataFrame([])

    # read ensembl ids
    id_map_df = pd.read_table(id_map, header=None, index_col=1,
                              names=['Protein_id'])
    ens_id_obj = id_map_df.index.difference(meta_df.index)

    # check if othologues download is finished
    if ens_id_obj.empty:
        print('Othologues download finished!')
        orthologues_dl_df = pd.read_table(download_table, index_col=0)
    else:
        # start to fetch othologues information
        for each_id in ens_id_obj:
            each_id_orth_df = get_ensembl_orthologues(each_id)
            each_id_orth_df = pd.merge(
                id_map_df, each_id_orth_df,
                left_index=True, right_index=True)
            meta_df = pd.concat([meta_df, each_id_orth_df])
            with open(meta_file, 'a') as meta_file_inf:
                for i in range(len(each_id_orth_df)):
                    meta_file_inf.write('{ei}\t{pi}\t{ot}\t{tp}\t{qp}\n'.format(
                        ei=each_id, pi=each_id_orth_df.ix[i][0],
                        ot=each_id_orth_df.ix[i][1],
                        tp=each_id_orth_df.ix[i][2],
                        qp=each_id_orth_df.ix[i][3]))
        orthologues_df = meta_df
        orthologues_df.index.name = 'Ensembl_id'
        orthologues_df.to_csv(meta_file, sep='\t', na_rep='None')

        # filter identify cutoff > 60
        orthologues_dl_df = orthologues_df.dropna()
        orthologues_dl_df = orthologues_dl_df[
            (orthologues_dl_df.Target_percent > 60) &
            (orthologues_dl_df.Query_percent > 60)]
        orthologues_dl_df.to_csv(download_table, sep='\t')

    # download othologues sequence
    # step check if sequence download is finished
    genomic_seq_file = os.path.join(
        out_dir, 'ensembl.othologues.genomic.seq.fa')
    downloaded_seq = list()
    if os.path.exists(genomic_seq_file):
        for each_rd in SeqIO.parse(genomic_seq_file, 'fasta'):
            downloaded_seq.append(each_rd.id)
    downloaded_seq_index = pd.Index(downloaded_seq)
    otho_idx = pd.Index(orthologues_dl_df.Orthologue)
    left_seq_obj = otho_idx.difference(downloaded_seq_index)
    # download sequences
    if left_seq_obj.empty:
        print('Othologues sequence download finished!')
    else:
        genomic_seq_inf = open(genomic_seq_file, 'a')
        for each_id in left_seq_obj.unique():
            ensembl_seq = get_ensembl_seq(each_id)
            if ensembl_seq:
                genomic_seq_inf.write(ensembl_seq)
        genomic_seq_inf.close()


if __name__ == '__main__':
    main()
