import os
import click
import pandas as pd


REPEAT_HEADER = [
    'chrom',
    'start',
    'end',
    'circRNAID',
    'score',
    'strand',
    'region',
    'repeat_chrom',
    'repeat_start',
    'repeat_end',
    'repeat_id',
    'repeat_score',
    'repeat_strand',
    'repeat_type',
    'repeat_class',
    'overlap'
]


REGION_MAP = {
    'up': 'flankIntronUpSINE',
    'down': 'flankIntronDownSINE',
}


def reapeat_type_stats(repeat_df):
    repeat_type_df = repeat_df.loc[:, ['region', 'repeat_class']]
    t_num = repeat_type_df.groupby(['region']).size()
    repeat_type_num = repeat_type_df.groupby(
        ['region'])['repeat_class'].value_counts()
    rp_portion = pd.DataFrame(repeat_type_num / t_num)
    rp_portion.columns = ['portion']
    return rp_portion


def get_sine_content(repeat_df):
    repeat_df.region.replace(REGION_MAP, inplace=True)
    sine_df = repeat_df[repeat_df.repeat_class == 'Type I Transposons/SINE']
    sine_counts = sine_df.groupby(
        ['circRNAID', 'region', 'repeat_class']).size()
    sine_counts = pd.DataFrame(sine_counts)
    sine_counts.columns = ['counts']
    sine_counts.index = sine_counts.index.droplevel('repeat_class')
    sine_counts = sine_counts.unstack('region')
    sine_counts.columns = sine_counts.columns.droplevel()
    return sine_counts


@click.command()
@click.argument(
    'repeat_overlap',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'name',
    type=click.STRING,
    required=True
)
@click.argument(
    'out_dir',
    type=click.Path(file_okay=False, exists=True),
    required=True
)
def main(repeat_overlap, name, out_dir):
    repeat_df = pd.read_table(repeat_overlap, header=None,
                              names=REPEAT_HEADER)
    # get repeat class portion
    rp_class_portion = reapeat_type_stats(repeat_df)
    rp_class_file = os.path.join(
        out_dir, '{n}.repeat.class.txt'.format(n=name))
    rp_class_portion.to_csv(rp_class_file, sep='\t')
    # get SINE content for each circRNA up/down stream flank intron
    sine_content_file = os.path.join(
        out_dir, '{n}.SINE.content.txt'.format(n=name)
    )
    sine_content_df = get_sine_content(repeat_df)
    sine_content_df.to_csv(sine_content_file, sep='\t', na_rep=0)


if __name__ == '__main__':
    main()
