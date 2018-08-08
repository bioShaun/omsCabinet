import click
import os
import pandas as pd
import envoy


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


def blat(database, query, blatout):
    print 'Running Blat, this may take a while...'
    r = envoy.run('/public/software/blat/blat {d} {q} {o}'.format(
        q=query, d=database, o=blatout
    ))
    return r.std_err


def df2gtf(blat_df, prefix, name, cutoff):
    blat_df.loc[:, 21] = 0
    # coverage
    blat_df.loc[:, 22] = blat_df.loc[:, 0] / blat_df.loc[:, 10]
    blat_df = blat_df.sort_values([9, 0])
    blat_df.loc[:, 22] = blat_df.loc[:, 22].round(3)
    bed_col = [13, 15, 16, 9, 21, 8, 15, 16, 21, 17, 18, 20, 22, 10]
    blat_bed_df = blat_df.loc[:, bed_col]
    # output gtf file
    mask = blat_bed_df.loc[:, 22] > cutoff
    pass_blat_gtfs = get_gtf_lines(blat_bed_df.loc[mask])
    failed_blat_gtfs = get_gtf_lines(blat_bed_df.loc[~mask])
    blat_gtf_file = '{pre}.{name}.cutoff_{cut}.passed.gtf'.format(
        pre=prefix, cut=cutoff, name=name)
    blat_failed_gtf_file = '{pre}.{name}.cutoff_{cut}.failed.gtf'.format(
        pre=prefix, cut=cutoff, name=name)
    write_obj(pass_blat_gtfs, blat_gtf_file)
    write_obj(failed_blat_gtfs, blat_failed_gtf_file)


@click.command()
@click.argument(
    'database',
    type=click.Path(dir_okay=False, exists=True),
    required=True,
)
@click.argument(
    'query',
    type=click.Path(dir_okay=False, exists=True),
    required=True,
)
@click.option(
    '-c',
    '--cutoff',
    default=0.8,
    type=click.FLOAT,
    help='Blat coverage cutoff (0-1), default is 0.8.'
)
def main(database, query, cutoff):
    blat_out = '{query}.blatout'.format(query=query)
    if not os.path.exists(blat_out):
        blat(database, query, blat_out)
    blat_out_prefix = os.path.splitext(blat_out)[0]
    blat_df = pd.read_table(blat_out, skiprows=5, header=None)
    blat_df = blat_df.drop_duplicates()
    max_map = blat_df.groupby([9])[0].idxmax()
    blat_best_df = blat_df.loc[max_map]
    df2gtf(blat_df, blat_out_prefix, 'all', cutoff)
    df2gtf(blat_best_df, blat_out_prefix, 'best', cutoff)


if __name__ == '__main__':
    main()
