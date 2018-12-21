import pandas as pd
import gtfparse
import click


@click.command()
@click.option(
    '-g',
    '--gtf',
    help='gtf file.',
    required=True,
    type=click.Path(dir_okay=False)
)
@click.option(
    '-o',
    '--output',
    help='output gene length file.',
    required=True,
    type=click.Path(dir_okay=False)
)
def main(gtf, output):
    gtf_df = gtfparse.read_gtf(gtf)
    gtf_exon_df = gtf_df[gtf_df.feature == 'exon']
    gtf_exon_df.loc[:, 'exon_len'] = gtf_exon_df.end - gtf_exon_df.start + 1
    tr_len = gtf_exon_df.groupby(['transcript_id'])['exon_len'].sum()
    tr_gene = gtf_exon_df.loc[:, [
        'transcript_id', 'gene_id']].drop_duplicates()
    tr_gene = tr_gene.set_index('transcript_id')
    tr_gene_len = pd.concat([tr_len, tr_gene], axis=1)
    gene_len = tr_gene_len.groupby(['gene_id'])['exon_len'].median()
    gene_len.to_csv(output, header=False, sep='\t')


if __name__ == '__main__':
    main()
