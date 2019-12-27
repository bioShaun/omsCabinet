import pandas as pd
import numpy as np
import os
import fire
from pathlib import PurePath, Path


CUTOFF = 0.2


def min_df_value_except0(exp_df):
    tmp_df = exp_df.replace({0: np.inf})
    return tmp_df.min().min()


def filter_stats(matrix, cut_range=range(1, 11),
                 out_dir=None, group_file=None,
                 multi=1):
    exp_df = pd.read_table(matrix, index_col=0)
    if group_file is not None:
        group_df = exp_sample2group(matrix, group_file, by='min')
    for each_cut in cut_range:
        each_cut = each_cut * multi
        if group_file is not None:
            passed_genes = group_df[group_df.T.max() > each_cut].index
            exp_df = exp_df.loc[passed_genes]
        else:
            exp_df = exp_df[exp_df.T.max() > each_cut]
        gene_num = len(exp_df)
        print(f'cutoff {each_cut}: {gene_num}')
        if out_dir is not None:
            out_dir = Path(out_dir)
            out_matrix_name = PurePath(
                matrix).with_suffix(f'.cutoff{each_cut}.txt').name
            out_matrix = out_dir / out_matrix_name
            if not out_dir.exists():
                out_dir.mkdir()
            exp_df.to_csv(out_matrix, sep='\t')


def filter_exp_by_por(matrix, outfile, exp_cutoff, prop_cutoff):
    exp_df = pd.read_table(matrix, index_col=0)
    if exp_cutoff == 0:
        exp_cutoff = min_df_value_except0(exp_df)
    f_exp_df = exp_df > exp_cutoff
    f_exp_prop = f_exp_df.sum(1) / exp_df.shape[1]
    p_exp_df = exp_df.loc[f_exp_prop > prop_cutoff]
    p_exp_df.to_csv(outfile, sep='\t', float_format='%.3f')


def gene_group_long_table(matrix, outfile, cutoff=0, group_file=None):
    if group_file:
        exp_df = exp_sample2group(matrix, group_file, by='max')
        exp_df = exp_df.reset_index()
    else:
        exp_df = pd.read_table(matrix)
    if cutoff == 0:
        cutoff = min_df_value_except0(exp_df)
    m_exp_df = exp_df.melt(id_vars=exp_df.columns[0],
                           value_name='tpm', var_name='tissue')
    m_exp_df = m_exp_df[m_exp_df.tpm > cutoff]
    m_exp_df.to_csv(outfile, sep='\t', index=False)


def exp_sample2group(matrix, group, outfile=None, cutoff=0, by='mean'):
    exp_df = pd.read_table(matrix, index_col=0)
    if cutoff == 0:
        cutoff = min_df_value_except0(exp_df)
    exp_df = exp_df[exp_df.T.max() > cutoff]
    group_df = pd.read_table(group, index_col=1, header=None)
    group_df.columns = ['tissue']
    group_exp_df = pd.merge(exp_df.T, group_df,
                            left_index=True, right_index=True,
                            how='left')
    group_mean_exp_df = group_exp_df.groupby('tissue').agg(by).T
    group_mean_exp_df.index.name = 'Gene_id'
    if outfile is None:
        return group_mean_exp_df
    else:
        group_mean_exp_df.to_csv(outfile, sep='\t',
                                 float_format='%.3f')


def cal_fc_ts(matrix, group, gene_classify,
              out_dir, cutoff=0.1):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    gene_df = pd.read_table(gene_classify, index_col=0)
    group_df = pd.read_table(group, index_col=1, header=None)
    tissue_num = len(group_df.tissue.unique())
    ts_score_cutoff = tissue_num * CUTOFF
    # mean_exp_value = group_mean_exp_df.T.mean()
    group_mean_exp_df = exp_sample2group(matrix, group, cutoff=cutoff)
    sum_exp_value = group_mean_exp_df.T.sum()
    max_exp_value = group_mean_exp_df.T.max()
    ts_score = max_exp_value / sum_exp_value
    max_exp_tissue = group_mean_exp_df.T.idxmax()
    ts_df = pd.concat([ts_score, max_exp_tissue], axis=1)
    ts_df.columns = ['ts_score', 'tissue']
    ts_df.index.name = 'Gene_id'
    ts_df.loc[:, 'ts'] = [
        each > ts_score_cutoff for each in ts_df.loc[:, 'ts_score']]
    ts_by_gene_df = pd.merge(ts_df, gene_df,
                             left_index=True, right_index=True, how='left')
    ts_out = os.path.join(out_dir, 'tissue_specific.score.txt')
    ts_num_out = os.path.join(out_dir, 'tissue_specific.number.txt')
    ts_num_summary_out = os.path.join(
        out_dir, 'tissue_specific.number.summary.txt')
    ts_gene_num = ts_by_gene_df.groupby(['tissue', 'gene_biotype'])['ts'].sum()
    ts_gene_num_df = ts_gene_num.reset_index()
    ts_gene_num_df.to_csv(ts_num_out, sep='\t', index=False)
    ts_by_gene_df.to_csv(ts_out, sep='\t')
    gene_num = ts_by_gene_df.gene_biotype.value_counts()
    ts_num = ts_by_gene_df.groupby('gene_biotype')['ts'].sum()
    ts_over_all_summary = pd.concat([gene_num, ts_num], axis=1)
    ts_over_all_summary.columns = ['detected_genes', 'ts_genes']
    ts_over_all_summary.loc[:, 'ts_portion'] = ts_over_all_summary.ts_genes / \
        ts_over_all_summary.detected_genes
    ts_over_all_summary.index.name = 'Gene_type'
    ts_over_all_summary.to_csv(ts_num_summary_out, sep='\t',
                               float_format='%.3f')


if __name__ == '__main__':
    fire.Fire()
