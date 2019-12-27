#! /usr/bin/env python3

import fire
import pandas as pd
import numpy as np
import sys
import pathlib


class Correlation(object):

    def __init__(self,
                 sample_tissue_file, out_dir=pathlib.Path.cwd(),
                 method='pearson', exp_table=None,
                 cor_table=None, gene_type_file=None):
        self.exp_table = exp_table
        self.gene_type_file = gene_type_file
        self.sample_tissue_file = sample_tissue_file
        self.out_dir = out_dir
        self.method = method
        self.cor_table = cor_table

    def cal_cor(self):
        # exp_df = pd.read_table(self.exp_table,
        #                        index_col=0)
        # mask = exp_df.T.mean()
        # log_exp_df = np.log10()
        return self.cor_df

    def cor_classify(self, gene_type):
        cor_classify_file = pathlib.PurePath(self.out_dir) / \
            '{gt}.cor.classify.txt'.format(gt=gene_type)

        if self.cor_table:
            self.cor_df = pd.read_table(self.cor_table,
                                        index_col=0)
        elif self.exp_table:
            self.cal_cor()
        else:
            print('Epression table or Correlation table must be provided.')
            sys.exit(1)
        self.cor_df.loc[:, 'sample_id1'] = self.cor_df.index
        melt_cor_df = pd.melt(self.cor_df, id_vars=[
                              'sample_id1'], value_vars=self.cor_df.index,
                              value_name='cor_value',
                              var_name='sample_id2')
        melt_cor_df = melt_cor_df[melt_cor_df.sample_id1 !=
                                  melt_cor_df.sample_id2]
        melt_cor_df.loc[:, 'gene_type'] = gene_type

        def get_sample_index(sample, sep='_'):
            sample_idx = sample.split(sep)[-1]
            if sep == '_':
                try:
                    int(sample_idx)
                except ValueError:
                    return get_sample_index(sample, '-')
            else:
                try:
                    int(sample_idx)
                except ValueError:
                    sample_idx = np.nan
            return sample_idx

        melt_cor_df.loc[:, 'sample_idx1'] = melt_cor_df.sample_id1.map(
            get_sample_index)
        melt_cor_df.loc[:, 'sample_idx2'] = melt_cor_df.sample_id2.map(
            get_sample_index)

        self.tissue_df = pd.read_table(self.sample_tissue_file, index_col=0,
                                       header=None, names=['tissue'])
        melt_cor_df = melt_cor_df.merge(self.tissue_df,
                                        left_on='sample_id1',
                                        right_index=True)
        melt_cor_df = melt_cor_df.merge(self.tissue_df,
                                        left_on='sample_id2',
                                        right_index=True,
                                        suffixes=('1', '2'))
        melt_cor_df.loc[:, 'category'] = 'None'
        mask = melt_cor_df.sample_idx1 == melt_cor_df.sample_idx2
        melt_cor_df.loc[mask, 'category'] = 'within_tissue'
        mask = melt_cor_df.tissue1 == melt_cor_df.tissue2
        melt_cor_df.loc[mask, 'category'] = 'within_sample'
        melt_cor_df.to_csv(cor_classify_file, sep='\t', header=None,
                           index=False)


if __name__ == '__main__':
    fire.Fire(Correlation)
