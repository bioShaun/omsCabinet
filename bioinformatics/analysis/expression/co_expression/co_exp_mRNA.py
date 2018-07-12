import pandas as pd
import click


@click.command()
@click.option(
    '-g',
    '--gene_module',
    help='gene module table.',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '-c',
    '--cor_table',
    help='correlation table.',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '-o',
    '--output',
    help='module correlated mRNA results.',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '-n',
    '--col_num',
    help='target gene column number',
    default=2,
    type=click.INT
)
def main(gene_module, cor_table, output, col_num):
    index_col = col_num - 1
    gm_df = pd.read_table(gene_module)
    cor_df = pd.read_table(cor_table, header=None,
                           index_col=index_col,
                           names=['mRNA', 'cor', 'pvalue', 'padj'])
    cor_gm_df = pd.merge(cor_df, gm_df,
                         left_index=True, right_on='gene_id')
    cor_gm_mean = cor_gm_df.groupby(['mRNA', 'module_name'])['cor'].mean()
    cor_gm_sum = cor_gm_df.groupby(['mRNA', 'module_name'])['cor'].sum()
    cor_total = cor_gm_df.groupby(['mRNA'])['cor'].sum()
    cor_gm_size = cor_gm_df.groupby(['mRNA', 'module_name'])['cor'].size()
    cor_total_size = cor_gm_df.groupby(['mRNA']).size()
    cor_gm_other = cor_total - cor_gm_sum
    cor_gm_other_size = cor_total_size - cor_gm_size
    cor_gm_other_cor = cor_gm_other / cor_gm_other_size
    cor_gm_diff = cor_gm_mean - cor_gm_other_cor
    cor_gm_summary = pd.concat(
        [cor_gm_mean, cor_gm_other_cor, cor_gm_diff], axis=1)
    cor_gm_summary.columns = ['module_cor', 'non_module_cor', 'diff']
    cor_gm_summary.to_csv(output, sep='\t')


if __name__ == '__main__':
    main()
