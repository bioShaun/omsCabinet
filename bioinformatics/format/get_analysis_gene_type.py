import click
from gtfparse import read_gtf
from omstools.utils.config import gtf_tools
import os


@click.command()
@click.option(
    '--linc',
    is_flag=True,
    help='seperate lncRNA to lincRNA/genic_lncRNA'
)
@click.argument(
    'gtf',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'out_dir',
    type=click.Path(exists=True, file_okay=False),
    required=True
)
def main(gtf, out_dir, linc):
    gtf_df = read_gtf(gtf)
    linc_genes = gtf_df[gtf_df.transcript_biotype ==
                        'lincRNA'].gene_id.unique()
    gtf_df.gene_biotype.replace(gtf_tools['dict_GENCODE_CATEGORY_MAP'],
                                inplace=True)
    gtf_df.transcript_biotype.replace(
        gtf_tools['dict_GENCODE_CATEGORY_MAP'], inplace=True
    )
    if 'gene_name' in gtf_df.columns:
        mask = (gtf_df.gene_name == "")
        gtf_df.loc[mask, 'gene_name'] = gtf_df.loc[mask, 'gene_id']
    else:
        gtf_df.loc[:, 'gene_name'] = gtf_df.loc[:, 'gene_id']
    if linc:
        gtf_df = gtf_df.set_index('gene_id')
        gtf_df.loc[linc_genes, 'gene_biotype'] = 'lincRNA'
        gtf_df.gene_biotype.replace(
            {'lncRNA': 'genic_lncRNA'}, inplace=True)
        gtf_df = gtf_df.reset_index()
    gene_df = gtf_df[gtf_df.gene_id != ""]
    gene_type_df = gene_df.loc[:, [
        'gene_id', 'gene_name', 'gene_biotype']].drop_duplicates()
    gene_type_file = os.path.join(out_dir, 'gene_type.txt')
    gene_type_df.to_csv(gene_type_file, sep='\t', index=False)
    tr_df = gtf_df[gtf_df.transcript_id != ""]
    tr_type_df = tr_df.loc[:, ['transcript_id', 'gene_id',
                               'gene_name', 'transcript_biotype',
                               'gene_biotype']].drop_duplicates()
    tr_type_file = os.path.join(out_dir, 'transcript_type.txt')
    tr_type_df.to_csv(tr_type_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
