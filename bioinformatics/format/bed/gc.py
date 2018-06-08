#! /usr/bin/env python3

import fire
from Bio import SeqIO
import pandas as pd
import envoy
import pathlib
import shutil
import sys


CURRENT_DIR = pathlib.Path.cwd()
BED_HEADER = [
    'chrom',
    'start',
    'end',
    'feature_id',
    'score',
    'strand',
    'gene_id',
    'feature'
]


def cmd_exists(cmd):
    return shutil.which(cmd) is not None


class GCcontent(object):

    def __init__(self, genome, feature_bed, out_dir=CURRENT_DIR):
        self.genome = genome
        self.bed = pathlib.PurePath(feature_bed)
        self.fa = self.bed.with_suffix('.fa')
        self.gc_df = None
        self.gc_file = out_dir / 'gc.content.txt'

    def get_feature_fa(self):
        '''
        using bedtools to extract fasta from bedfile
        '''
        if pathlib.Path(self.fa).is_file():
            print('Fasta file exists.')
            return 1
        if cmd_exists('bedtools'):
            bed2fa_cmd = 'bedtools getfasta -fi {t.genome} \
-fo {t.fa} -name -bed {t.bed} -s'.format(
                t=self
            )
            envoy.run(bed2fa_cmd)
        else:
            sys.exit('Not find bedtools!')
        return 1

    def get_gc_count(self):
        '''
        calculate GC content for fasta file
        '''
        gc_dict = dict()
        for each_rd in SeqIO.parse(self.fa, 'fasta'):
            g_count = each_rd.seq.lower().count('g')
            c_count = each_rd.seq.lower().count('c')
            each_rd_id = each_rd.id.split('(')[0]
            each_tr_id = each_rd.id.split('|')[0]
            gc_dict.setdefault('feature_id', []).append(each_rd_id)
            gc_dict.setdefault('gc_count', []).append(g_count + c_count)
            gc_dict.setdefault('transcript_id', []).append(each_tr_id)
        self.gc_df = pd.DataFrame(gc_dict)

    def get_gc_stats(self):
        '''
        calculate GC portion for each region of transcripts
        '''
        bed_df = pd.read_table(self.bed, header=None,
                               names=BED_HEADER)
        merged_df = pd.merge(bed_df, self.gc_df)
        merged_df.loc[:, 'length'] = merged_df.end - merged_df.start
        agg_dict = {
            'length': 'sum',
            'gc_count': 'sum',
        }
        merged_df_group = merged_df.groupby(
            ['transcript_id', 'gene_id', 'feature']).agg(agg_dict)
        merged_df_group.loc[:, 'GC'] = merged_df_group.gc_count / \
            merged_df_group.length
        merged_df_group.to_csv(self.gc_file, sep='\t', columns=['GC'])

    def pipeline(self):
        print('Generating fasta file.')
        self.get_feature_fa()
        print('Calculating GC content.')
        self.get_gc_count()
        print('Calculating GC portion for each transcripts.')
        self.get_gc_stats()
        print('Analysis finished.')


if __name__ == '__main__':
    fire.Fire(GCcontent)
