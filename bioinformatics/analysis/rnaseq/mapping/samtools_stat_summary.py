import fire
import pandas as pd
from io import StringIO
from pathlib import PurePath


def extract_label_content(samtools_stats, label):
    content_list = []
    with open(samtools_stats) as stats_inf:
        for eachline in stats_inf:
            eachline = eachline.strip()
            if eachline.startswith(label):
                if '#' in eachline:
                    eachline_inf = eachline.strip().split('\t')
                    content_list.append('\t'.join(eachline_inf[:-1]))
                else:
                    content_list.append(eachline)
    return content_list


def extract_mapping_stats(samtools_stats, fai):
    sn_list = extract_label_content(samtools_stats, 'SN')
    sn_stringio = StringIO('\n'.join(sn_list))
    sn_df = pd.read_table(sn_stringio, header=None,
                          index_col=1, names=['lable', 'stats_value'])
    mapped_reads = int(sn_df.loc['reads mapped:'].stats_value)
    total_reads = int(sn_df.loc['raw total sequences:'].stats_value)
    mapped_rate = mapped_reads / total_reads
    mapped_bases = sn_df.loc['bases mapped (cigar):'].stats_value
    fai_df = pd.read_table(fai, header=None)
    genome_length = fai_df.loc[:, 1].sum()
    average_depth = mapped_bases / genome_length

    cov_list = extract_label_content(samtools_stats, 'COV')
    cov_stringio = StringIO('\n'.join(cov_list))
    cov_df = pd.read_table(cov_stringio, header=None, names=[
                           'label', 'cov_lab', 'cov_count', 'base_count'])
    cov_df.loc[:, 'base_num'] = cov_df.cov_count * cov_df.base_count
    cov4 = 1 - (cov_df[cov_df.cov_count < 4].base_num.sum() /
                mapped_bases)
    return mapped_reads, mapped_rate, average_depth, cov4


def merge_mapping_stats(mapping_dir, sample_inf, fa_idx,
                        outfile, suffix='realn.bam.stat'):
    sample_df = pd.read_table(sample_inf,
                              header=None,
                              index_col=0)
    mapping_dir = PurePath(mapping_dir)
    mapping_files = [mapping_dir / f'{each}/{each}.{suffix}'
                     for each in sample_df.index]
    fa_idxs = [fa_idx for each in mapping_files]
    mapping_stats = list(
        map(extract_mapping_stats, mapping_files, fa_idxs))
    mapping_df = pd.DataFrame(mapping_stats,
                              columns=['Mapped_Reads', 'Mapping_Rate',
                                       'Average_Depth', 'Coverage_4X'])
    mapping_df.index = sample_df.index
    mapping_df.index.name = 'Sample_ID'
    mapping_df.to_csv(outfile, sep='\t',
                      float_format='%.3f')


if __name__ == '__main__':
    fire.Fire(merge_mapping_stats)
