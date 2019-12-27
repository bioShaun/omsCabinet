import gtfparse
import click
import pandas as pd
import os
from collections import OrderedDict


def gene_tr_count_from_gtf(gtf_file):
    gtf_df = gtfparse.read_gtf(gtf_file)
    gene_num = len(gtf_df.gene_id.dropna().unique())
    tr_num = len(gtf_df.transcript_id.dropna().unique())
    return gene_num, tr_num


@click.command()
@click.argument(
    'assemblyline_dir',
    type=click.Path(file_okay=False, exists=True),
    required=True
)
@click.argument(
    'summary_file',
    type=click.Path(exists=False),
    required=True
)
@click.argument(
    'species',
    type=click.STRING,
    required=True
)
def main(assemblyline_dir, summary_file, species):
    lib_file = os.path.join(assemblyline_dir, 'library_id.map')
    lib_df = pd.read_table(lib_file, header=None, index_col=0,
                           names=['sample_id'])
    exp_num_dict = OrderedDict()
    for each_lib in lib_df.index:
        each_lib_exp_gtf = os.path.join(assemblyline_dir,
                                        'classify',
                                        '{lb}.expr.gtf'.format(
                                            lb=each_lib
                                        ))
        each_lib_bgr_gtf = os.path.join(assemblyline_dir,
                                        'classify',
                                        '{lb}.bkgd.gtf'.format(
                                            lb=each_lib
                                        ))
        exp_gene, exp_tr = gene_tr_count_from_gtf(each_lib_exp_gtf)
        bgr_gene, bgr_tr = gene_tr_count_from_gtf(each_lib_bgr_gtf)
        exp_num_dict.setdefault(
            'expressed_gene_number', []).append(exp_gene)
        exp_num_dict.setdefault(
            'expressed_transcript_number', []).append(exp_tr)
        exp_num_dict.setdefault(
            'background_gene_number', []).append(bgr_gene)
        exp_num_dict.setdefault(
            'background_transcript_number', []).append(bgr_tr)
    exp_num_df = pd.DataFrame(exp_num_dict, index=lib_df.sample_id)
    exp_num_df.index.name = 'sample_id'
    exp_num_df.loc[:, 'species'] = species
    exp_num_df.to_csv(summary_file, sep='\t', header=False)


if __name__ == '__main__':
    main()
