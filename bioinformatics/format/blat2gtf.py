import click
import os
import pandas as pd


def bedline2gtf(bedline_split, chr_len):
    gtf_lines = list()
    eachline = list(bedline_split)
    chrom = eachline[0]
    start = int(eachline[1])
    end = int(eachline[2])
    f_id = eachline[3]
    strand = eachline[5]
    coverage = eachline[12]
    # if strand == '-':
    #     start, end = end, start
    blocksize = [int(each) for each in eachline[10].split(',')
                 if each]
    blockoffset = [int(each) for each in eachline[11].split(',')
                   if each]
    if strand == '-':
        blockoffset = list(reversed([(chr_len - each)
                                     for each in blockoffset]))
        blocksize = reversed(blocksize)
    feature = 'transcript_id "{tr_id}"; coverage "{cov}";'.format(
        tr_id=f_id, cov=coverage)
    # write transcript line
    transcript_line = '{ch}\tbed2gtf\ttranscript\t{st}\t{ed}\t.\t{std}\t.\t{ft}\n'.format(
        ch=chrom, st=start + 1, ed=end, std=strand, ft=feature
    )
    gtf_lines.append(transcript_line)
    for m, each_size in enumerate(blocksize):
        each_offset = blockoffset[m]
        each_start = each_offset
        each_end = each_start + each_size
        if strand == '-':
            each_end = each_offset
            each_start = each_end - each_size
        each_exon = '{ch}\tbed2gtf\texon\t{st}\t{ed}\t.\t{std}\t.\t{ft}\n'.format(
            ch=chrom, st=each_start + 1, ed=each_end, std=strand, ft=feature
        )
        gtf_lines.append(each_exon)
    return gtf_lines


def get_gtf_lines(blat_df):
    gtf_lines = list()
    for each_index in blat_df.index:
        bedlines = list(blat_df.loc[each_index])
        chr_len = bedlines[-1]
        bed_inf = bedlines[:-1]
        gtf_lines.extend(bedline2gtf(bed_inf, chr_len))
    return gtf_lines


def write_obj(obj, outfile):
    if obj:
        with open(outfile, 'w') as out_inf:
            for each in obj:
                out_inf.write(each)


@click.command()
@click.argument(
    'blat_out',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.option(
    '-c',
    '--cutoff',
    default=0,
    type=click.FLOAT,
    help='Blat coverage cutoff (0-1), default is 0 : no limit.'
)
def main(blat_out, cutoff):
    blat_out_prefix = os.path.splitext(blat_out)[0]
    blat_gtf_file = '{pre}.cutoff_{cut}.pass.gtf'.format(
        pre=blat_out_prefix, cut=cutoff)
    blat_failed_gtf_file = '{pre}.cutoff_{cut}.failed.gtf'.format(
        pre=blat_out_prefix, cut=cutoff)
    blat_df = pd.read_table(blat_out, skiprows=5, header=None)
    max_map = blat_df.groupby([13])[0].idxmax()
    blat_best_df = blat_df.loc[max_map]
    blat_best_df.loc[:, 21] = 0
    # coverage
    blat_best_df.loc[:, 22] = blat_best_df.loc[:, 0] / blat_best_df.loc[:, 14]
    blat_best_df.loc[:, 22] = blat_best_df.loc[:, 22].round(3)
    bed_col = [9, 11, 12, 13, 21, 8, 11, 12, 21, 17, 18, 19, 22, 10]
    blat_bed_df = blat_best_df.loc[:, bed_col]
    # output gtf file
    mask = blat_bed_df.loc[:, 22] > cutoff
    pass_blat_gtfs = get_gtf_lines(blat_bed_df.loc[mask])
    failed_blat_gtfs = get_gtf_lines(blat_bed_df.loc[~mask])
    write_obj(pass_blat_gtfs, blat_gtf_file)
    write_obj(failed_blat_gtfs, blat_failed_gtf_file)


if __name__ == '__main__':
    main()
