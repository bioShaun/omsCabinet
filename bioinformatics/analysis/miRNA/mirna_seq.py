import os
import click
import pandas as pd
from Bio import SeqIO


def get_sp_miseq(sp_obj, all_seq):
    '''get a subset of miRNA fasta belong specific species
       from miRbase mature/hairpin fasta
    '''
    out_list = list()
    for each_rd in SeqIO.parse(all_seq, 'fasta'):
        each_seq_sp = each_rd.id.split('-')[0]
        if each_seq_sp in sp_obj:
            out_list.append(each_rd)
    return out_list


def output_seq(seq_list, outfile):
    if seq_list:
        SeqIO.write(seq_list, outfile, "fasta")
    else:
        print 'No sequence for file {out}.'.format(out=outfile)


@click.command()
@click.argument(
    'all_sp_table',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'mirbase_dir',
    type=click.Path(file_okay=False, exists=True),
    required=True
)
@click.argument(
    'abbr',
    type=click.STRING,
    required=True
)
@click.argument(
    'out_dir',
    type=click.Path(file_okay=False),
    default=os.getcwd()
)
def main(mirbase_dir, abbr, out_dir, all_sp_table):
    all_sp_df = pd.read_table(all_sp_table, index_col=0)
    mature_fa = os.path.join(mirbase_dir, 'mature.fa')
    hairpin_fa = os.path.join(mirbase_dir, 'hairpin.fa')
    if abbr in all_sp_df.index:
        other_sp = all_sp_df.index.drop(abbr)
    else:
        other_sp = all_sp_df.index
    sp_mature_seqs = get_sp_miseq([abbr], mature_fa)
    other_mature_seqs = get_sp_miseq(other_sp, mature_fa)
    sp_hairpin_seqs = get_sp_miseq([abbr], hairpin_fa)
    sp_mature_seq_file = os.path.join(
        out_dir, '{sp}.mature.fa'.format(sp=abbr))
    other_mature_seq_file = os.path.join(
        out_dir, '{sp}-other.mature.fa'.format(sp=abbr))
    sp_hairpin_seq_file = os.path.join(
        out_dir, '{sp}.hairpin.fa'.format(sp=abbr))
    output_seq(sp_mature_seqs, sp_mature_seq_file)
    output_seq(other_mature_seqs, other_mature_seq_file)
    output_seq(sp_hairpin_seqs, sp_hairpin_seq_file)


if __name__ == '__main__':
    main()
