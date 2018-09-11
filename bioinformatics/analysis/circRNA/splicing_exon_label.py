import fire
import pandas as pd


def splicing_exon_position(circ_bed, transcritp_bed, outfile):
    circ_df = pd.read_table(circ_bed)
    tr_df = pd.read_table(transcritp_bed, header=None)
    tr_total_exon = tr_df.loc[:, [3, 9]]
    tr_total_exon.columns = ['isoformName', 'isoform_exon']
    circ_df = circ_df.merge(tr_total_exon)
    splicing_exon_df = circ_df.loc[:, ['chrom', 'start', 'end',
                                       'strand', 'index', 'isoform_exon']
                                   ].drop_duplicates()

    def splicing_exon(x):
        return int(x.split(',')[0])

    splicing_exon_df.loc[:, 'splicing_exon'] = splicing_exon_df.loc[
        :, 'index'].map(splicing_exon)

    def splicing_exon_pos(exon_num, total_exon, strand):
        if strand == '-':
            exon_num = total_exon - exon_num + 1
        if exon_num == 1:
            return 'first'
        elif exon_num == total_exon:
            return 'last'
        else:
            return 'middle'

    splicing_exon_df.loc[:, 'splicing_exon_name'] = list(
        map(splicing_exon_pos, splicing_exon_df.splicing_exon,
            splicing_exon_df.isoform_exon, splicing_exon_df.strand))

    splicing_exon_df.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    fire.Fire(splicing_exon_position)
