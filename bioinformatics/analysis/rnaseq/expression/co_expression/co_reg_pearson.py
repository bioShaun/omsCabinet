import click
import pandas as pd
import random
import numpy as np


@click.command()
@click.option(
    '-l',
    '--lnc_pc_pair',
    help='lncRNA vs protein_coding gene pair',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '-p',
    '--pc_pc_pair',
    help='protein_coding vs protein_coding gene pair',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '--pc_cor',
    help='protein_coding vs protein_coding pearson correlation matrix.',
    required=True,
    type=click.Path(dir_okay=False)
)
@click.option(
    '--lnc_cor',
    help='lncRNA vs protein_coding pearson correlation matrix.',
    required=True,
    type=click.Path(dir_okay=False)
)
def main(lnc_pc_pair, pc_pc_pair, pc_cor, lnc_cor):
    lnc_pc_pair_df = pd.read_table(lnc_pc_pair, index_col=[0, 1],
                                   header=None)
    lnc_pc_pair_df.columns = ['direction']
    pc_pc_pair_df = pd.read_table(pc_pc_pair, index_col=[0, 1],
                                  header=None)
    pc_pc_pair_df.columns = ['direction']
    lnc_pc_cor_df = pd.read_table(lnc_cor, index_col=[0, 1],
                                  header=None)
    lnc_pc_cor_df.columns = ['cor', 'pvalue', 'qvalue']
    pc_pc_cor_df = pd.read_table(pc_cor, index_col=[0, 1],
                                 header=None)
    pc_pc_cor_df.columns = ['cor', 'pvalue', 'qvalue']
    cor_list = list()
    for each_dir in lnc_pc_pair_df.direction.unique():
        each_dir_lnc_pc = lnc_pc_pair_df[
            lnc_pc_pair_df.direction == each_dir].index.values
        each_dir_pc_pc = pc_pc_pair_df[
            pc_pc_pair_df.direction == each_dir].index.values
        np.random.seed(0)
        each_random_lnc_pc = np.random.choice(
            lnc_pc_cor_df.index.values, len(each_dir_lnc_pc))
        np.random.seed(0)
        each_random_pc_pc = np.random.choice(
            each_dir_pc_pc, len(each_dir_lnc_pc))
        each_dir_lnc_pc_cor_mean = pd.DataFrame(
            lnc_pc_cor_df.loc[each_dir_lnc_pc].cor)
        each_dir_lnc_pc_cor_mean.loc[:, 'direction'] = each_dir
        each_dir_lnc_pc_cor_mean.loc[:, 'type'] = 'lncRNA-protein_coding'
        each_dir_pc_pc_cor_mean = pd.DataFrame(
            pc_pc_cor_df.loc[each_random_pc_pc].cor)
        each_dir_pc_pc_cor_mean.loc[:, 'direction'] = each_dir
        each_dir_pc_pc_cor_mean.loc[
            :, 'type'] = 'protein_coding-protein_coding'
        each_random_lnc_pc_cor_mean = pd.DataFrame(
            lnc_pc_cor_df.loc[each_random_lnc_pc].cor)
        each_random_lnc_pc_cor_mean.loc[:, 'direction'] = each_dir
        each_random_lnc_pc_cor_mean.loc[:, 'type'] = 'random'
        cor_list.extend([each_dir_lnc_pc_cor_mean,
                         each_dir_pc_pc_cor_mean,
                         each_random_lnc_pc_cor_mean])
    co_reg_cor_df = pd.concat(cor_list)
    co_reg_cor_df.to_csv('co_reg_correlation.txt', sep='\t')

if __name__ == '__main__':
    main()
