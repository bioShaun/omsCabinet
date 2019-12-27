import click
import pandas as pd
import numpy as np


def get_cor(gene_pair, exp_df):
    cor_list = list()
    for gene1, gene2 in gene_pair:
        gene1_exp = exp_df.loc[gene1]
        gene2_exp = exp_df.loc[gene2]
        cor_val = np.corrcoef(gene1_exp, gene2_exp)[0][1]
        cor_list.append(cor_val)
    return cor_list


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
    '--pc_exp',
    help='protein coding expression table.',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '--lnc_exp',
    help='lncRNA expression table',
    type=click.Path(dir_okay=False),
    required=True
)
def main(lnc_pc_pair, pc_pc_pair, pc_exp, lnc_exp):
    lnc_exp_df = pd.read_table(lnc_exp, index_col=0)
    pc_exp_df = pd.read_table(pc_exp, index_col=0)
    merged_exp_df = pd.concat([lnc_exp_df, pc_exp_df])
    lnc_pc_pair_df = pd.read_table(lnc_pc_pair, index_col=[0, 1],
                                   header=None)
    lnc_pc_pair_df.columns = ['direction']
    pc_pc_pair_df = pd.read_table(pc_pc_pair, index_col=[0, 1],
                                  header=None)
    pc_pc_pair_df.columns = ['direction']
    cor_list = list()
    for each_dir in lnc_pc_pair_df.direction.unique():
        each_dir_lnc_pc = lnc_pc_pair_df[
            lnc_pc_pair_df.direction == each_dir].index.values
        each_dir_pc_pc = pc_pc_pair_df[
            pc_pc_pair_df.direction == each_dir].index.values
        each_random_lnc_pc_dict = dict()
        np.random.seed(0)
        each_random_lnc = np.random.choice(
            lnc_exp_df.index, 1000)
        np.random.seed(0)
        each_random_pc = np.random.choice(
            pc_exp_df.index, 1000)
        each_random_lnc_pc_dict['lncRNA'] = each_random_lnc
        each_random_lnc_pc_dict['mRNA'] = each_random_pc
        each_random_lnc_pc_dict['direction'] = [each_dir for each in range(1000)]
        each_random_lnc_pc_dict['type'] = ['random' for each in range(1000)]
        each_random_lnc_pc_df = pd.DataFrame(each_random_lnc_pc_dict)
        each_random_lnc_pc_df = each_random_lnc_pc_df.set_index(['lncRNA', 'mRNA'])
        np.random.seed(0)
        each_random_pc_pc = np.random.choice(
            each_dir_pc_pc, len(each_dir_lnc_pc))
        each_random_pc_pc_df = pc_pc_pair_df.loc[each_random_pc_pc]
        each_random_pc_pc_df.loc[:, 'type'] = 'protein_coding-protein_coding'
        each_lnc_pc_df = lnc_pc_pair_df.loc[each_dir_lnc_pc]
        each_lnc_pc_df.loc[:, 'type'] = 'lncRNA-protein_coding'
        cor_list.extend([each_lnc_pc_df, each_random_pc_pc_df, each_random_lnc_pc_df])
    co_reg_cor_df = pd.concat(cor_list)
    co_reg_cor_df.index.names = ['gene1', 'gene2']
    co_reg_cor_df.loc[:, 'cor'] = get_cor(co_reg_cor_df.index.values,
                                          merged_exp_df)
    co_reg_cor_df.to_csv('co_reg_correlation.txt', sep='\t')

if __name__ == '__main__':
    main()
