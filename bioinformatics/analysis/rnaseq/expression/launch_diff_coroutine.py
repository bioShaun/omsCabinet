#! /usr/bin/env python3

import fire
import asyncio
import pandas as pd
import itertools
import pathlib
from asyncio.subprocess import PIPE, STDOUT


SEMA = asyncio.Semaphore(8)
DIFF_R = '/project0/OM-mRNA-pig-limingzhou-P160901/sus_scrofa/diff/diff_analysis.R'


def get_compare_names(sample_inf, method='pairwise',
                      contrasts=None, replicates=True):
    if contrasts is None:
        group_sample_df = pd.read_table(
            sample_inf, header=None, index_col=0)
        compare_list = list(itertools.combinations(
            group_sample_df.index.unique(), 2))
        tmp_compare_list = compare_list[:]
        if replicates:
            for each_pair in tmp_compare_list:
                shape1 = len(group_sample_df.loc[each_pair[0]].shape)
                shape2 = len(group_sample_df.loc[each_pair[1]].shape)
                if shape1 < 2 or shape2 < 2:
                    compare_list.remove(each_pair)
        if method == 'pairwise':
            compare_name_list = ['{0}_vs_{1}'.format(
                each_compare[0], each_compare[1])
                for each_compare in compare_list]
        elif method == 'preferential':
            unique_group_list = set([item for sublist in compare_list
                                     for item in sublist])
            compare_name_list = ['{g}_vs_other'.format(g=each)
                                 for each in unique_group_list]
        else:
            pass
    else:
        contrasts_df = pd.read_table(
            contrasts, header=None)
        compare_name_list = ['{0}_vs_{1}'.format(contrasts_df.loc[i, 0],
                                                 contrasts_df.loc[i, 1])
                             for i in contrasts_df.index]
    return compare_name_list


async def launch_cmd(cmd):
    with (await SEMA):
        p = await asyncio.create_subprocess_shell(cmd,
                                                  stdout=PIPE,
                                                  stderr=STDOUT)
        return (await p.communicate())[0].splitlines()


def check_dir(diff_dir):
    diff_dir = pathlib.Path(diff_dir)
    if diff_dir.exists() and list(diff_dir.iterdir()):
        return True
    else:
        return False


class Diff(object):

    def __init__(self, sample_inf, counts,
                 tpm_table, out_dir, gene_class,
                 qvalue=0.05, thread=8,
                 fc=1, method='pairwise', lib_size=''):
        self.sample_inf = sample_inf
        self.counts = counts
        self.tpm_table = tpm_table
        self.out_dir = out_dir
        self.gene_class = gene_class
        self.qvalue = qvalue
        self.fc = fc
        self.method = method
        self.lib_size = lib_size
        self.compares = get_compare_names(sample_inf,
                                          method=method)
        self.thread = thread

    def launch_diff(self):
        if not pathlib.Path(self.out_dir).exists():
            pathlib.Path(self.out_dir).mkdir()
        cmd_list = list()
        for each_compare in self.compares:
            each_comp_dir = pathlib.PurePath(self.out_dir) / each_compare
            if check_dir(each_comp_dir):
                print('{comp} analysis is finished.'.format(comp=each_compare))
            else:
                each_cmd = 'Rscript {diff_r} --counts {t.counts} \
--compare {comp} --sample_inf {t.sample_inf} --tpm_table {t.tpm_table} \
-o {comp_dir} --qvalue {t.qvalue} --logfc {t.fc}'.format(
                    diff_r=DIFF_R, t=self,
                    comp=each_compare, comp_dir=each_comp_dir
                )
                if self.lib_size:
                    each_cmd = '{cmd} --lib_size {lib}'.format(
                        lib=self.lib_size, cmd=each_cmd
                    )
                cmd_list.append(each_cmd)
        if cmd_list:
            loop = asyncio.get_event_loop()
            f = asyncio.wait([launch_cmd(each) for each in cmd_list])
            loop.run_until_complete(f)
        else:
            print('all compare analysis is finished.')


if __name__ == '__main__':
    fire.Fire(Diff)
