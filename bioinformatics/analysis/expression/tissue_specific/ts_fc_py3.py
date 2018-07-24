import pandas as pd
import os
import fire


CUTOFF = 0.2


def mean_exp(matrix, group, outfile=None, cutoff=0):
    exp_df = pd.read_table(matrix, index_col=0)
    exp_df = exp_df[exp_df.T.max() > cutoff]
    group_df = pd.read_table(group, index_col=1, header=None)
    group_df.columns = ['tissue']
    group_exp_df = pd.merge(exp_df.T, group_df,
                            left_index=True, right_index=True,
                            how='left')
    group_mean_exp_df = group_exp_df.groupby('tissue').mean().T
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
    group_mean_exp_df = mean_exp(matrix, group, cutoff=cutoff)
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
