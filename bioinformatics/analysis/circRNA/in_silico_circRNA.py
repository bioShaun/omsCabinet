import click
import pandas as pd
import numpy as np
import gtfparse
import itertools
from collections import OrderedDict


def get_introns(exon_df):
    if len(exon_df.shape) == 1:
        # single exon transcript has no intron
        return None
    intron_list = list()
    e1 = exon_df.iloc[0]
    for j in xrange(1, len(exon_df)):
        e2 = exon_df.iloc[j]
        intron_list.append([e1.end, e2.start_0base])
        e1 = e2
    return intron_list


def get_flank_intron(chrom, exon_combine, introns):
    if exon_combine[0] > 0:
        li = '{c}:{s}-{e}'.format(c=chrom,
                                  s=introns[exon_combine[0] - 1][0],
                                  e=introns[exon_combine[0] - 1][1])
    else:
        li = 'None'
    try:
        ri = '{c}:{s}-{e}'.format(c=chrom,
                                  s=introns[exon_combine[1]][0],
                                  e=introns[exon_combine[1]][1])
    except IndexError:
        ri = 'None'
    return '{li}|{ri}'.format(li=li, ri=ri)


@click.command()
@click.argument(
    'input_gtf',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'in_silico_circ',
    type=click.Path(exists=False),
    required=True
)
@click.option(
    '-r',
    '--real_circ_table',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
def main(input_gtf, in_silico_circ, real_circ_table):
    real_circ_df = pd.read_table(real_circ_table)
    real_circ_df.loc[:, 'chrom'] = real_circ_df.chrom.astype(str)
    real_circ_df = real_circ_df.set_index(['chrom', 'start', 'end'])
    circ_num = len(real_circ_df)
    is_circ_dict = OrderedDict()
    gtf_df = gtfparse.read_gtf(input_gtf)
    exon_df = gtf_df[gtf_df.feature == 'exon']
    exon_df = exon_df.set_index('transcript_id')
    exon_df.loc[:, 'start_0base'] = exon_df.start - 1
    for each_tr in exon_df.index.unique():
        # in silico circRNA from same transcript set
        if each_tr not in real_circ_df.isoformName.unique():
            continue
        each_tr_exons = exon_df.loc[each_tr]
        each_tr_introns = get_introns(each_tr_exons)
        if each_tr_introns is None:
            continue
        exon_num = len(each_tr_exons)
        chrom = str(each_tr_exons.iloc[0].seqname)
        strand = each_tr_exons.iloc[0].strand
        gene = each_tr_exons.iloc[0].gene_id
        for each_com in itertools.combinations_with_replacement(
                range(exon_num), 2):
            start = each_tr_exons.iloc[each_com[0]].start_0base
            end = each_tr_exons.iloc[each_com[1]].end
            # filter real circRNA from in silico circRNA
            if (chrom, start, end) in real_circ_df.index:
                continue
            flank_intron = get_flank_intron(chrom, each_com, each_tr_introns)
            is_circ_dict.setdefault('chrom', []).append(chrom)
            is_circ_dict.setdefault('start', []).append(start)
            is_circ_dict.setdefault('end', []).append(end)
            is_circ_dict.setdefault('strand', []).append(strand)
            is_circ_dict.setdefault('flankIntron', []).append(flank_intron)
            is_circ_dict.setdefault('isoformName', []).append(each_tr)
            is_circ_dict.setdefault('geneID', []).append(gene)
    is_circ_df = pd.DataFrame(is_circ_dict)
    np.random.seed(0)
    selected_circ = np.random.choice(is_circ_df.index.values, circ_num)
    is_circ_df = is_circ_df.loc[selected_circ]
    is_circ_df.loc[:, 'circRNAID'] = ['in_silico_circ_{num:0>10}'.format(
        num=each + 1
    ) for each in range(len(is_circ_df))]
    is_circ_df = is_circ_df.set_index('circRNAID')
    is_circ_df.to_csv(in_silico_circ, sep='\t')


if __name__ == '__main__':
    main()
