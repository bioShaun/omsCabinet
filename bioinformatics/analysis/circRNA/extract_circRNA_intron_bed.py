#!/usr/bin/env python

########################################################
# Extract flank intron of circRNA to bed format
#
# ------------------------------------------------------
# output upstream and downstrean intron to same file
# using "up", "down" to label them
# 2015.5.15
########################################################

import pandas as pd
import click
import os


def get_bed_line(intron_line):
    if intron_line != 'None':
        chrom, pos = intron_line.split(':')
        start, end = pos.split('-')
        return '{c}\t{s}\t{e}'.format(
            c=chrom, s=start, e=end
        )
    else:
        return None


@click.command()
@click.argument(
    'circ_bed',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'out_file',
    type=click.Path(dir_okay=False, exists=False),
    required=True
)
def main(circ_bed, out_file):
    circ_df = pd.read_table(circ_bed)
    circ_intron_list = [each.split('|')
                        for each in circ_df.flankIntron]
    intron_inf = open(out_file, 'w')
    for n, each in enumerate(circ_intron_list):
        if len(each) == 2:
            strand = circ_df.strand[n]
            circ_id = circ_df.circRNAID[n]
            if strand == '+':
                up_intron_line, down_intron_line = each
            else:
                up_intron_line, down_intron_line = list(reversed(each))
            up_intron_line = get_bed_line(up_intron_line)
            down_intron_line = get_bed_line(down_intron_line)
            if up_intron_line:
                intron_inf.write('{il}\t{ci}\t.\t{st}\tup\n'.format(
                    il=up_intron_line,
                    ci=circ_id,
                    st=strand))
            if down_intron_line:
                intron_inf.write('{il}\t{ci}\t.\t{st}\tdown\n'.format(
                    il=down_intron_line,
                    ci=circ_id,
                    st=strand))
    intron_inf.close()


if __name__ == '__main__':
    main()
