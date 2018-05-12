from __future__ import division
import pandas as pd
import numpy as np
import click
import glob
import os
import sys


READ_CUTOFF = 2
SAMPLE_CUTOFF = 1

CIRC_HEADER = [
    'chrom',
    'start',
    'end',
    'name',
    'score',
    'strand',
    'thickStart',
    'thickEnd',
    'itemRgb',
    'exonCount',
    'exonSizes',
    'exonOffsets',
    'readNumber',
    'circType',
    'geneName',
    'isoformName',
    'index',
    'flankIntron'
]


OUT_COL = [
    'circ_name',
    'chrom',
    'start',
    'end',
    'strand',
    'exonCount',
    'length',
    'flankIntron',
    'flankIntronSizeUP',
    'flankIntronSizeDOWN',
    'circType',
    'isoformName',
    'geneName',
    'gene_name',
    'transcript_biotype',
    'sample_id',
    'readNumber'
]

OUT_COL_NAME = [
    'circRNAID',
    'chrom',
    'start',
    'end',
    'strand',
    'exonCount',
    'length',
    'flankIntron',
    'flankIntronSizeUP',
    'flankIntronSizeDOWN',
    'circType',
    'isoformName',
    'geneID',
    'geneSymbol',
    'transcriptBiotype',
    'sampleID',
    'readNumber'
]


SP_EN_NAME = {
    'mus_musculus': 'Mouse',
    'canis_familiaris': 'Dog',
    'cavia_porcellus': 'Guinea_pig',
    'felis_catus': 'Cat',
    'macaca_mulatta': 'Monkey',
    'oryctolagus_cuniculus': 'Rabbit',
    'rattus_norvegicus': 'Rat',
    'ovis_aries': 'Sheep',
    'gallus_gallus': 'Chicken',
    'sus_scrofa': 'Pig',
}


def read_circ_table(circ_dir, suffix='circularRNA_known.txt'):
    circ_tables = glob.glob('{d}/*{s}'.format(d=circ_dir, s=suffix))
    circ_df_list = list()
    for each_file in circ_tables:
        each_df = pd.read_table(each_file, header=None, names=CIRC_HEADER)
        sp_name = os.path.basename(each_file).split('.')[0]
        each_df.loc[:, 'sample_id'] = sp_name
        circ_df_list.append(each_df)
    circ_df = pd.concat(circ_df_list)
    return circ_df


def get_circ_stat(stats_df, col_name, method='mean'):
    all_num = stats_df.loc[:, col_name].apply(method)
    by_tr_type = stats_df.groupby(
        ['transcript_biotype']).agg({col_name: method})
    by_tr_type.loc['Total', col_name] = all_num
    return by_tr_type


def get_circ_num(stats_df):
    all_num = pd.Series(len(stats_df), index=['Total'])
    by_tr_type = stats_df.transcript_biotype.value_counts()
    # circ_gene_type_num = stats_df.groupby(['gene_biotype', 'circType']).size()
    # circ_gene_type_num.name = 'circ_count'
    # circ_gene_type_num = circ_gene_type_num.reset_index()
    # circ_gene_type_num.loc[:, 'name'] = circ_gene_type_num.circType + \
    #     '-' + circ_gene_type_num.gene_biotype
    # by_gene_circ_type = pd.Series(
    #     circ_gene_type_num.circ_count.values, index=circ_gene_type_num.name)
    circ_num_stats = pd.concat(
        [all_num, by_tr_type]
    )
    circ_num_stats.name = 'number'
    return circ_num_stats


def get_host_gene_num(stats_df):
    all_num = pd.Series(len(stats_df.geneName.unique()), index=['Total'])
    ids_by_circ_type = stats_df.groupby('transcript_biotype')[
        'geneName'].unique()
    by_tr_type = pd.Series(map(len, ids_by_circ_type),
                           index=ids_by_circ_type.index)
    # circ_gene_type_num = stats_df.groupby(['gene_biotype', 'circType'])[
    #     'geneName'].unique()
    # num_by_gene_circ_type = map(len, circ_gene_type_num)
    # circ_gene_type_num = circ_gene_type_num.reset_index()
    # circ_gene_type_num.loc[:, 'name'] = circ_gene_type_num.circType + \
    #     '-' + circ_gene_type_num.gene_biotype
    # by_gene_circ_type = pd.Series(
    #     num_by_gene_circ_type, index=circ_gene_type_num.name
    # )
    host_gene_stats = pd.concat(
        [all_num, by_tr_type]
    )
    host_gene_stats.name = 'hostGene'
    return host_gene_stats


def exonSizes_to_len(exonsizes):
    exons = [int(each) for each in exonsizes.split(',')]
    return sum(exons)


def circ_flank_intron(circ_df):

    def flankIntron2size(flankIntron, strand):
        # 1:7089216-7120193|1:7163371-7169514
        intron_list = flankIntron.split('|')
        if len(intron_list) == 1:
            return [np.nan, np.nan]
        try:
            intron_cor_list = [each.split(':')[1].split('-')
                               if each != 'None'
                               else np.nan
                               for each in intron_list]
        except IndexError:
            print intron_list
            sys.exit(1)
        try:
            intron_size = [int(each[1]) - int(each[0])
                           if isinstance(each, list)
                           else each
                           for each in intron_cor_list]
        except TypeError:
            print intron_cor_list
            sys.exit(1)
        if strand == '-':
            intron_size = intron_size[::-1]
        return intron_size
    try:
        tmp = map(flankIntron2size,
                  circ_df.flankIntron,
                  circ_df.strand)
        circ_df.loc[:, 'flankIntronSizeUP'] = [each[0] for each
                                               in tmp]
        circ_df.loc[:, 'flankIntronSizeDOWN'] = [each[1] for each
                                                 in tmp]
        # circ_df.loc[:, 'flankIntronSize'] = tmp
        # circ_df.loc[:, 'flankIntronSize'] = map(flankIntron2size,
        #                                         circ_df.flankIntron,
        #                                         circ_df.strand)
    except ValueError:
        print len(circ_df.flankIntron), len(circ_df.strand), len(circ_df)
        sys.exit(1)
    return circ_df
    # circ_df.loc[:, 'flankIntronSizeUP'] = [each[0] for each
    #                                        in circ_df.flankIntronSize]
    # circ_df.loc[:, 'flankIntronSizeDOWN'] = [each[1] for each
    #                                          in circ_df.flankIntronSize]


@click.command()
@click.argument(
    'circ_dir',
    type=click.Path(file_okay=False, exists=True),
    required=True
)
@click.argument(
    'gene_type',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'tissue_sample',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'mapping_summary',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'exp_table',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'out_dir',
    type=click.Path(file_okay=False, exists=True),
    default=os.getcwd()
)
@click.option(
    '-t',
    '--circ_type',
    type=click.Choice(['all_circ', 'circRNA', 'ciRNA']),
    default='all_circ'
)
@click.option(
    '-s',
    '--species',
    type=click.STRING,
    help='species latin',
    required=True
)
@click.option(
    '-a',
    '--abbr',
    type=click.STRING,
    help='species abbr',
    required=True
)
def main(circ_dir, gene_type, out_dir, species, exp_table,
         tissue_sample, mapping_summary, circ_type, abbr):
    sp_en_name = SP_EN_NAME[species]
    combined_circ_out = os.path.join(circ_dir, 'circ.combined.txt')
    if not os.path.isfile(combined_circ_out):
        circ_df = read_circ_table(circ_dir)
        circ_df.to_csv(combined_circ_out, index=False, sep='\t')
    else:
        circ_df = pd.read_table(combined_circ_out)

    circ_df = circ_df[circ_df.readNumber >= READ_CUTOFF]
    if circ_type != 'all_circ':
        circ_df = circ_df[circ_df.circType == circ_type]
    gene_type_df = pd.read_table(gene_type)
    circ_type_df = pd.merge(circ_df, gene_type_df,
                            left_on='isoformName', right_on='transcript_id',
                            how='left')
    circ_type_df = circ_flank_intron(circ_type_df)

    def get_circ_basic_stats(circ_type_df):
        stats_type_df = circ_type_df.loc[:, ['chrom', 'start', 'end',
                                             'circType', 'exonSizes',
                                             'exonCount', 'transcript_biotype',
                                             'geneName', 'flankIntronSizeUP',
                                             'flankIntronSizeDOWN']].drop_duplicates()
        stats_type_df.loc[:, 'length'] = stats_type_df.exonSizes.map(
            exonSizes_to_len)
        circ_num_stats = get_circ_num(stats_type_df)
        exon_stats = get_circ_stat(stats_type_df,
                                   'exonCount')
        circ_len_stats = get_circ_stat(stats_type_df,
                                       'length')
        f_intron_up_stats = get_circ_stat(stats_type_df,
                                          'flankIntronSizeUP')
        f_intron_down_stats = get_circ_stat(stats_type_df,
                                            'flankIntronSizeDOWN')
        host_gene_stats = get_host_gene_num(stats_type_df)
        junc_reads_stats = get_circ_stat(circ_type_df,
                                         'readNumber', 'sum')
        circ_merged_stats = pd.concat(
            [circ_num_stats, exon_stats, circ_len_stats,
             f_intron_up_stats, f_intron_down_stats,
             host_gene_stats, junc_reads_stats], axis=1)
        circ_merged_stats.index.name = 'Category'
        return circ_merged_stats
    circ_merged_stats = get_circ_basic_stats(circ_type_df)
    stats_out_file = os.path.join(out_dir, '{sp}.{t}.stats.txt'.format(
        t=circ_type,
        sp=sp_en_name))
    circ_merged_stats.to_csv(stats_out_file, sep='\t', float_format='%.2f',
                             na_rep='None')

    # sample summary
    sample_stats_file = os.path.join(
        out_dir, '{t}.stats.sample.txt'.format(t=circ_type))
    sample_stats_list = list()
    for each_sample in circ_type_df.sample_id.unique():
        each_sample_df = circ_type_df[circ_type_df.sample_id == each_sample]
        each_sample_stats = get_circ_basic_stats(each_sample_df)
        each_sample_stats.loc[:, 'sample_id'] = each_sample
        sample_stats_list.append(each_sample_stats)
    sample_df = pd.concat(sample_stats_list)
    # sample_df = sample_df.reset_index().set_index(['sample_id', 'Category'])
    # sample_out_df = sample_df.unstack('Category')

    # mapping reads info
    mapping_summary_df = pd.read_table(mapping_summary)
    mapping_summary_df = mapping_summary_df.loc[:, ['Sample', 'total']]
    mapping_summary_df.columns = ['sample_id', 'Total_reads']
    sample_out_df = sample_df.reset_index().set_index(
        ['sample_id', 'Category'])
    sample_out_df = sample_out_df.unstack(
        'Category').fillna(0).stack('Category')
    sample_out_df = pd.merge(sample_out_df.reset_index(), mapping_summary_df,
                             how='left')
    sample_out_df.loc[:, 'readNumberPortion(%)'] = sample_out_df.readNumber / \
        sample_out_df.Total_reads * 100
    sample_out_df.loc[:, 'Species'] = species
    sample_out_df = sample_out_df.set_index('Species')
    sample_out_df.to_csv(sample_stats_file, sep='\t',
                         float_format='%.3f', na_rep=0,
                         header=False)

    # tissue summary
    tissue_df = pd.read_table(tissue_sample, names=[
                              'tissue', 'sample_id'], header=None)
    tissue_sample_num = tissue_df.groupby('tissue').size()
    circ_tissue_type_df = pd.merge(circ_type_df, tissue_df,
                                   how='left')
    # circ_ids = ['{c}:{s}-{e}'.format(c=circ_tissue_type_df.loc[each, 'chrom'],
    #                                  s=circ_tissue_type_df.loc[each, 'start'],
    #                                  e=circ_tissue_type_df.loc[each, 'end'])
    #             for each in circ_tissue_type_df.index]
    # circ_tissue_type_df.loc[:, 'circ_id'] = circ_ids
    # multi_sample_tissues = tissue_sample_num[tissue_sample_num > 1].index
    # if multi_sample_tissues.empty:
    #     pass
    # else:
    #     circ_by_tissue = circ_tissue_type_df.groupby(
    #         ['tissue', 'circ_id']).size()
    #     circ_by_tissue_part1 = circ_by_tissue[circ_by_tissue > READ_CUTOFF]
    #     circ_by_tissue_part1_df = circ_tissue_type_df.set_index(
    #         ['tissue', 'circ_id'])
    #     circ_by_tissue_part1_df = circ_by_tissue_part1_df.loc[
    #         circ_by_tissue_part1.index].reset_index()
    #     single_sample_tissues = tissue_sample_num[tissue_sample_num == 1].index
    #     circ_by_tissue_part2_df = circ_tissue_type_df[
    #         circ_tissue_type_df.tissue.isin(single_sample_tissues)]
    #     circ_tissue_type_df = pd.concat(
    #         [circ_by_tissue_part1_df, circ_by_tissue_part2_df])
    # tissue_stats_list = list()
    # tissue_stats_file = os.path.join(
    #     out_dir, '{t}.stats.tissue.txt'.format(t=circ_type))
    # for each_tissue in tissue_sample_num.index:
    #     each_tissue_df = circ_tissue_type_df[
    #         circ_tissue_type_df.tissue == each_tissue]
    #     each_tissue_stats = get_circ_basic_stats(each_tissue_df)
    #     each_tissue_stats.loc[:, 'tissue'] = each_tissue
    #     tissue_stats_list.append(each_tissue_stats)
    # tissue_stats_df = pd.concat(tissue_stats_list)
    # tissue_stats_df = tissue_stats_df.reset_index(
    # ).set_index(['tissue', 'Category'])
    # tissue_out_df = tissue_stats_df.unstack('Category')
    # tissue_out_df.to_csv(tissue_stats_file, sep='\t',
    #                      float_format='%.2f', na_rep=0)

    # tissue circ_table
    circ_tissue_type_df = circ_tissue_type_df.sort_values(
        ['chrom', 'start', 'end'])
    circ_name = circ_tissue_type_df.groupby(['chrom', 'start', 'end']).size()
    circ_name.name = 'circ_count'
    circ_name = circ_name.reset_index()
    circ_name.loc[:, 'circ_name'] = [
        '{sp}_circ_{num:0>6}'.format(sp=abbr,
                                     num=each + 1)
        for each in circ_name.index]
    circ_name_df = pd.merge(circ_tissue_type_df, circ_name)
    circ_name_df.loc[:, 'length'] = circ_name_df.exonSizes.map(
        exonSizes_to_len)
    exp_table_df = pd.read_table(exp_table, index_col=0)

    def get_circ_exp_df(circ_df):
        each_tissue_df = circ_df.loc[:, OUT_COL]
        each_tissue_df.columns = OUT_COL_NAME
        each_tissue_df = each_tissue_df.set_index(OUT_COL_NAME[:-1])
        each_tissue_out_df = each_tissue_df.unstack('sampleID')
        each_tissue_out_df.columns = each_tissue_out_df.columns.droplevel()
        each_tissue_out_df.columns.name = ''
        each_tissue_out_df = each_tissue_out_df.fillna(0)
        each_tissue_out_df = each_tissue_out_df.astype('int')
        return each_tissue_out_df

    def combine_gene_exp(circ_df, exp_table_df, samples):
        each_tissue_host_df = exp_table_df.loc[:, samples]
        each_tissue_host_df.columns = ['{cn}(host gene)'.format(cn=each)
                                       for each in each_tissue_host_df.columns]
        each_tissue_out_df = pd.merge(circ_df.reset_index(),
                                      each_tissue_host_df,
                                      left_on='geneID', right_index=True,
                                      how='left')
        return each_tissue_out_df

    for each_tissue in tissue_sample_num.index:
        each_tissue_df = circ_name_df[circ_name_df.tissue == each_tissue]
        each_tissue_out_df = get_circ_exp_df(each_tissue_df)
        each_tissue_samples = sorted(tissue_df[tissue_df.tissue ==
                                               each_tissue].sample_id)
        each_tissue_out_df = combine_gene_exp(
            each_tissue_out_df, exp_table_df, each_tissue_samples)
        tissue_stats_file = os.path.join(
            out_dir, '{sp}.{ts}.{tp}.detail.txt'.format(tp=circ_type,
                                                        ts=each_tissue,
                                                        sp=sp_en_name))
        each_tissue_out_df.to_csv(tissue_stats_file, sep='\t',
                                  index=False, float_format='%.3f',
                                  na_rep='None')
    # merged information
    merged_out_file = os.path.join(
        out_dir, '{sp}.{tp}.detail.txt'.format(tp=circ_type,
                                               sp=sp_en_name))
    merged_out_df = get_circ_exp_df(circ_name_df)
    samples = sorted(circ_name_df.sample_id.unique())
    merged_out_df = combine_gene_exp(merged_out_df, exp_table_df, samples)
    tissue_num_df = pd.DataFrame(
        circ_name_df.groupby('circ_name')['tissue'].unique())
    tissue_num_df.loc[:, 'tissueCount'] = [
        len(each) for each in tissue_num_df.tissue]
    tissue_num_df.loc[:, 'tissue'] = [
        ','.join(each) for each in tissue_num_df.tissue]
    merged_out_df = pd.merge(merged_out_df, tissue_num_df,
                             left_on='circRNAID',
                             right_index=True,
                             how='left')
    merged_out_df.to_csv(merged_out_file, sep='\t',
                         index=False, float_format='%.3f',
                         na_rep='None')


if __name__ == '__main__':
    main()
