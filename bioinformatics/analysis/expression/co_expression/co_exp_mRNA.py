import pandas as pd
import click

@click.command()
@click.option(
    '-t',
    '--ts_score',
    help='tissue specific score table.',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '-c',
    '--cor_table',
    help='mRNA lncRNA correlation table.',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '-o',
    '--output',
    help='tissue correlated mRNA results.',
    type=click.Path(dir_okay=False),
    required=True
)
def main(ts_score, cor_table, output):
    ts_score_df = pd.read_table(ts_score, index_col=0)
    cor_df = pd.read_table(cor_table, header=None)
    cor_df.columns = ['lncRNA', 'mRNA', 'cor', 'pvalue', 'padj']
    mrnas = cor_df.mRNA.unique()
    cor_df = cor_df.set_index(['mRNA', 'lncRNA'])
    tissues = ts_score_df.tissue.unique()
    mrna_cor_dict = {}
    for each_ts in tissues:
        each_ts_lnc = ts_score_df[ts_score_df.cluster_tissue == each_ts].index
        each_non_ts_lnc = ts_score_df[ts_score_df.cluster_tissue != each_ts].index
        for each_mrna in mrnas:
            each_ts_cor = cor_df.loc[each_mrna].loc[each_ts_lnc].cor.mean()
            each_non_ts_cor = cor_df.loc[each_mrna].loc[each_non_ts_lnc].cor.mean()
            mrna_cor_dict.setdefault('mRNA', []).append(each_mrna)
            mrna_cor_dict.setdefault('tissue', []).append(each_ts)
            mrna_cor_dict.setdefault('tissue_cor', []).append(each_ts_cor)
            mrna_cor_dict.setdefault('non_tissue_cor', []).append(each_non_ts_cor)
            mrna_cor_dict.setdefault('diff', []).append(each_ts_cor - each_non_ts_cor)
    mrna_cor_df = pd.DataFrame(mrna_cor_dict)
    mrna_cor_out_df = mrna_cor_df.loc[:, ['mRNA', 'tissue', 'tissue_cor', 'non_tissue_cor', 'diff']]
    mrna_cor_out_df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    main()

