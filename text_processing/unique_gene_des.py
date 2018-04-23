import pandas as pd
import numpy as np
import click


def join_col(obj, sep=',', ignore_char='--'):
    if isinstance(obj, str):
        return obj
    else:
        clean_obj = [each for each in obj if each != ignore_char]
        return sep.join(clean_obj)


@click.command()
@click.argument(
    'gene_des',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'output',
    type=click.Path(dir_okay=False, exists=False),
    required=True
)
def main(gene_des, output):
    gene_des_df = pd.read_table(gene_des)
    gene_col = gene_des_df.columns[0]
    groupdf = gene_des_df.groupby([gene_col])
    des_df = groupdf.aggregate(np.unique)
    des_df = des_df.applymap(join_col)
    des_df.to_csv(output, sep='\t')


if __name__ == '__main__':
    main()
