import fire
import gtfparse
import pandas as pd


def gtf2genebed(gtf, bed):
    gtf_df = gtfparse.read_gtf(gtf)
    gene_start = gtf_df.groupby(['gene_id'])['start'].min()
    gene_end = gtf_df.groupby(['gene_id'])['end'].max()
    gene_chr = gtf_df.groupby(['gene_id'])['seqname'].first()
    gene_strand = gtf_df.groupby(['gene_id'])['strand'].first()
    gene_bed_df = pd.concat([gene_chr, gene_start, gene_end, gene_strand],
                            axis=1)
    gene_bed_df = gene_bed_df.reset_index()
    gene_bed_df.loc[:, 'score'] = '.'
    gene_bed_df.to_csv(
        bed,
        sep='\t',
        columns=['seqname', 'start', 'end', 'gene_id', 'score', 'strand'],
        index=False,
        header=False)


if __name__ == "__main__":
    fire.Fire(gtf2genebed)
