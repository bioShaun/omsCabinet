#! /usr/bin/env python3

import fire
import pandas as pd
import shutil
import pathlib
import envoy
import sys


CURRENT_DIR = pathlib.Path.cwd()


def cmd_exists(cmd):
    return shutil.which(cmd) is not None


class Repeat(object):

    def __init__(self, bedfile, repeat_file, out_dir=CURRENT_DIR):
        self.bed = pathlib.PurePath(bedfile)
        self.repeat = repeat_file
        self.out_dir = out_dir
        self.repeat_ovlp = self.bed.with_suffix('.repeat.overlap')
        self.repeat_file = out_dir / 'repeat.content.txt'

    def get_repeat_overlap(self):
        if pathlib.Path(self.repeat_ovlp).is_file():
            print('Finished repeat intersection analysis.')
            return 1
        if cmd_exists('bedtools'):
            repeat_is_cmd = 'bedtools intersect -a {t.bed} \
-b {t.repeat} -wao'.format(t=self)
            run_inf = envoy.run(repeat_is_cmd)
            with open(self.repeat_ovlp, 'w') as ovlp_inf:
                ovlp_inf.write(run_inf.std_out)
            if run_inf.std_err:
                sys.exit(run_inf.std_err)
        else:
            sys.exit('bedtools is required!')

    def get_repeat_stats(self):
        repeat_df = pd.read_table(self.repeat_ovlp, header=None)
        repeat_df = repeat_df.loc[:, [0, 1, 2, 3, 6, 7, 12]]
        repeat_df.columns = ['chrom', 'start', 'end',
                             'feature_id', 'gene_id',
                             'feature', 'repeat_len']

        def get_tr_id(feature_id):
            return feature_id.split('|')[0]

        repeat_df.loc[:, 'length'] = repeat_df.end - repeat_df.start
        repeat_df.loc[:, 'transcript_id'] = repeat_df.feature_id.map(get_tr_id)
        agg_dict = {
            'repeat_len': 'sum',
            'length': 'sum',
        }
        repeat_df_group = repeat_df.groupby(
            ['transcript_id', 'gene_id',
             'feature', 'feature_id', 'length'])['repeat_len'].sum()
        repeat_df_group = repeat_df_group.reset_index()
        repeat_df_group = repeat_df_group.groupby(
            ['transcript_id', 'gene_id', 'feature']).agg(agg_dict)
        repeat_df_group.loc[:, 'repeat'] = repeat_df_group.repeat_len / \
            repeat_df_group.length
        repeat_df_group.to_csv(self.repeat_file, sep='\t', columns=['repeat'])

    def pipeline(self):
        print('Get repeat/transcript overlap region.')
        self.get_repeat_overlap()
        print('Calculating repeat portion for each feature in transcripts.')
        self.get_repeat_stats()
        print('Repeat analysis finished.')


if __name__ == '__main__':
    fire.Fire(Repeat)
