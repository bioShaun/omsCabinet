from Bio import SeqIO
import click
from collections import Counter


@click.command()
@click.argument(
    'input_seq',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'dedup_seq',
    type=click.Path(exists=False),
    required=True
)
def main(input_seq, dedup_seq):
    my_records = []
    seq_counter = Counter()
    for seq_record in SeqIO.parse(input_seq, "fasta"):
        seq_id = seq_record.id
        seq_counter.update([seq_id])
        if seq_counter[seq_id] == 1:
            my_records.append(seq_record)
    SeqIO.write(my_records, dedup_seq, "fasta")


if __name__ == '__main__':
    main()
