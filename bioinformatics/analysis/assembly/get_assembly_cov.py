import pandas as pd
import click
import gtfparse
from omstools.utils.config import gtf_tools


@click.command()
@click.argument(
    'gtf_file',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'compare_table',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'species',
    type=click.STRING,
    required=True
)
@click.argument(
    'cov_file',
    type=click.Path(exists=False),
    required=True
)
def main(gtf_file, compare_table, cov_file, species):
    gtf_df = gtfparse.read_gtf(gtf_file)
    gene_df = gtf_df[gtf_df.feature == 'gene']
    gene_type_df = gene_df.loc[:, [
        'gene_id', 'gene_biotype']].drop_duplicates()
    gene_type_df = gene_type_df.set_index('gene_id')
    gene_type_df.gene_biotype.replace(gtf_tools['dict_GENCODE_CATEGORY_MAP'],
                                      inplace=True)
    gene_type_counts = gene_type_df.gene_biotype.value_counts()
    compare_table_df = pd.read_table(compare_table)
    assembly_genes = list()
    ref_assembly_df = compare_table_df[compare_table_df.category_relative ==
                                       'exonic_overlap']
    for each in ref_assembly_df.ref_gene_id:
        assembly_genes.extend(each.split(','))
    assembly_genes = list(set(assembly_genes))
    assembly_gene_type_df = gene_type_df.loc[assembly_genes]
    assembly_gene_type_counts = assembly_gene_type_df.gene_biotype.value_counts()
    merged_df = pd.concat(
        [gene_type_counts, assembly_gene_type_counts], axis=1)
    merged_df.columns = ['referrence', 'assembly']
    merged_df.loc[:, 'coverage'] = merged_df.assembly / merged_df.referrence
    merged_df.loc[:, 'species'] = species
    merged_df.to_csv(cov_file, sep='\t', header=False)


if __name__ == '__main__':
    main()
