import click
import pandas as pd


@click.command()
@click.argument(
    'exp_table',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'cutoff',
    type=click.FLOAT,
    default=1
)
def main(exp_table, cutoff):
    exp_df = pd.read_table(exp_table, index_col=0)
    ubi_genes = exp_df.T.min() > cutoff
    ubi_genes = ubi_genes[ubi_genes].index
    for each in ubi_genes:
        print each
    

if __name__ == '__main__':
    main()
