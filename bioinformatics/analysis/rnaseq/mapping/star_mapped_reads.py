import click
import pandas as pd
import os


STAR_COLUMNS = ['Uniquely mapped reads number',
                'Number of reads mapped to multiple loci',
                'Uniquely mapped reads %',
                '% of reads mapped to multiple loci']


def read_star_mapping_log(star_log_file):
    star_df = pd.read_table(star_log_file, header=None, sep='|', index_col=0)
    star_df = star_df.dropna()
    star_df = star_df.ix[4:]
    star_df.loc[:, 1] = [each.strip() for each in star_df.loc[:, 1]]
    return star_df


@click.command()
@click.option('-s', '--sample_inf', type=click.Path(exists=True),
              help='sample information, second column is sample id.',
              required=True,)
@click.option('-d', '--mapping_dir', type=click.Path(exists=True),
              help='mapping directory.', required=True)
def main(sample_inf, mapping_dir):
    sample_df = pd.read_table(sample_inf, header=None, index_col=1)
    out_file = os.path.join(mapping_dir, 'star_mapping.number.txt')
    star_log_file_list = [os.path.join(
        mapping_dir, each_sample, 'Log.final.out')
        for each_sample in sample_df.index]
    star_log_df_list = map(read_star_mapping_log, star_log_file_list)
    star_log_df = pd.concat(star_log_df_list, axis=1)
    star_log_df.columns = sample_df.index
    star_log_out_df = star_log_df.T
    star_log_out_df.index.name = 'Sample'
    star_log_out_df.columns = [each.strip()
                               for each in star_log_out_df.columns]
    out_data = star_log_out_df.loc[:, STAR_COLUMNS]
    out_data.columns = ['unique_mapped', 'multi_mapped',
                        'unique_mapped_rate', 'multi_mapped_rate']
    out_data.loc[:, 'unique_mapped'] = out_data.unique_mapped.astype('int')
    out_data.loc[:, 'multi_mapped'] = out_data.multi_mapped.astype('int')
    out_data.loc[:, 'total_mapped_reads'] = out_data.unique_mapped + \
        out_data.multi_mapped

    def rm_percent_lab(percent_num):
        return float(percent_num.rstrip('%'))

    out_data.loc[:, 'unique_mapped_rate'] = out_data.unique_mapped_rate.map(
        rm_percent_lab)
    out_data.loc[:, 'multi_mapped_rate'] = out_data.multi_mapped_rate.map(
        rm_percent_lab)
    out_data.loc[:, 'total_mapped_rate'] = out_data.unique_mapped_rate + \
        out_data.multi_mapped_rate

    out_data.to_csv(out_file, sep='\t', float_format='%.2f')


if __name__ == '__main__':
    main()
