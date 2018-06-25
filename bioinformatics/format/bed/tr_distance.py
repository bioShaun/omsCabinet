#! /usr/bin/env python3

import pandas as pd
import pathlib
import fire
import numpy as np


BEDHEADER = [
    'chrom',
    'start',
    'end',
    'tr_id',
    'score',
    'strand',
    'thickStart',
    'thickEnd',
    'itemRgb',
    'blockCount',
    'blockSizes',
    'blockStarts'
]


class Bed(object):

    def __init__(self, bedfile):
        self.bed = pathlib.PurePath(bedfile)
        self.bed_df = pd.read_table(bedfile,
                                    header=None,
                                    names=BEDHEADER)

    def nearby_tr(self, bed_df,
                  tr_idx, num=1, distance=None,
                  gene_type=None):
        tr_inf = bed_df.loc[tr_idx]
        same_chr_tr = bed_df[bed_df.chrom == tr_inf.chrom]
        same_chr_tr = bed_df[bed_df.gene_id != tr_inf.gene_id]
        if gene_type:
            same_chr_tr = same_chr_tr[
                same_chr_tr.gene_biotype == gene_type]
        up_tr_inf = same_chr_tr[same_chr_tr.index < tr_idx]
        up_tr_inf = up_tr_inf.loc[up_tr_inf.index[::-1][:num]]
        down_tr_inf = same_chr_tr[same_chr_tr.index > tr_idx]
        down_tr_inf = down_tr_inf.loc[down_tr_inf.index[:num]]
        return up_tr_inf, down_tr_inf

    def get_distance(self, tr_a, tr_b):
        try:
            distance = max(tr_a.start - tr_b.end,
                           tr_b.start - tr_a.end,
                           0)
        except ValueError:
            print('tra:{a}\ntrb:{b}'.format(a=tr_a.start, b=tr_b.start))
            return None
        else:
            return distance

    def neareast_tr(self, tr_type_file, nearby_type, outfile=None):
        tr_type_df = pd.read_table(tr_type_file)
        bed_df = pd.merge(self.bed_df, tr_type_df,
                          left_on='tr_id',
                          right_on='transcript_id',
                          how='left')
        neareast_tr_dict = dict()
        for each_idx in bed_df.index:
            up_tr, down_tr = self.nearby_tr(bed_df, each_idx,
                                            gene_type=nearby_type)
            each_tr = bed_df.loc[each_idx]
            if each_tr.strand == '-':
                up_tr, down_tr = down_tr, up_tr
            up_dis = down_dist = np.inf
            up_id = down_id = 'None'
            if not up_tr.empty:
                up_dis = self.get_distance(each_tr,
                                           up_tr.iloc[0])
                up_id = up_tr.iloc[0].gene_id
            if not down_tr.empty:
                down_dist = self.get_distance(each_tr,
                                              down_tr.iloc[0])
                down_id = down_tr.iloc[0].gene_id
            neareast_tr_dict.setdefault(
                'trancript_id', []).append(each_tr.tr_id)
            neareast_tr_dict.setdefault(
                'gene_id', []).append(each_tr.gene_id)
            neareast_tr_dict.setdefault(
                'upstream_gene', []).append(up_id)
            neareast_tr_dict.setdefault(
                'upstream_distance', []).append(up_dis)
            neareast_tr_dict.setdefault(
                'downstream_gene', []).append(down_id)
            neareast_tr_dict.setdefault(
                'downstream_distance', []).append(down_dist)
            neareast_tr_dict.setdefault(
                'gene_biotype', []
            ).append(each_tr.transcript_biotype)
        self.neareast_tr_df = pd.DataFrame(neareast_tr_dict)
        if outfile:
            self.neareast_tr_df.to_csv(outfile, sep='\t', index=False,
                                       columns=['trancript_id',
                                                'gene_id',
                                                'upstream_gene',
                                                'upstream_distance',
                                                'downstream_gene',
                                                'downstream_distance',
                                                'gene_biotype'])

    def neareast_gene(self, tr_type_file, nearby_type, outfile=None):
        self.neareast_tr(tr_type_file, nearby_type)
        up_gene_df = self.neareast_tr_df.loc[:, ['gene_id',
                                                 'upstream_gene',
                                                 'upstream_distance',
                                                 'gene_biotype']]
        nearest_up_gene_idx = up_gene_df.groupby(
            ['gene_id'])['upstream_distance'].idxmin()
        nearest_up_df = up_gene_df.loc[nearest_up_gene_idx]
        down_gene_df = self.neareast_tr_df.loc[:, ['gene_id',
                                                   'downstream_gene',
                                                   'downstream_distance',
                                                   'gene_biotype']]
        nearest_dn_gene_idx = down_gene_df.groupby(
            ['gene_id'])['downstream_distance'].idxmin()
        nearest_down_df = down_gene_df.loc[nearest_dn_gene_idx]
        self.nearest_gene_df = pd.merge(nearest_up_df, nearest_down_df)
        if outfile:
            self.nearest_gene_df.to_csv(outfile, sep='\t', index=False,
                                        columns=['gene_id',
                                                 'gene_biotype',
                                                 'upstream_gene',
                                                 'upstream_distance',
                                                 'downstream_gene',
                                                 'downstream_distance'])


if __name__ == '__main__':
    fire.Fire(Bed)
