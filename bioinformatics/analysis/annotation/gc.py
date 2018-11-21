from __future__ import division, print_function
import click
from Bio import SeqIO


class Bed(object):

    @staticmethod
    def from_string(line):
        f = Bed()
        # read bed line
        fields = line.strip().split('\t')
        f.seqid = fields[0]
        f.start = int(fields[1])
        f.end = int(fields[2])
        f.fieldid = fields[3]
        f.score = 0 if (fields[4] == '.') else float(fields[4])
        f.strand = fields[5]
        f.line = line.strip()
        return f

    @property
    def length(self):
        return self.end - self.start

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            # read the bed line
            if not line:
                continue
            if not line.strip():
                continue
            if line.startswith('#'):
                continue
            yield Bed.from_string(line)


@click.command()
@click.option(
    '-b',
    '--bed6',
    type=click.File('r'),
    help='bed6 file.',
    required=True,
)
@click.option(
    '-g',
    '--genome',
    type=click.Path(dir_okay=False),
    help='genome fasta.',
    required=True,
)
def main(bed6, genome):
    bed_dict = dict()
    for each_f in Bed.parse(bed6):
        bed_dict.setdefault(each_f.seqid, []).append(each_f)
    for seq_record in SeqIO.parse(genome, "fasta"):
        if seq_record.id not in bed_dict:
            continue
        for each_f in bed_dict[seq_record.id]:
            g_count = seq_record.seq[each_f.start:each_f.end].lower().count(
                'g')
            c_count = seq_record.seq[each_f.start:each_f.end].lower().count(
                'c')
            gc_count = g_count + c_count
            gc_por = gc_count / each_f.length
            # gc_count = g_count + c_count
            print('{line}\t{cnt}\t{por}'.format(
                line=each_f.line, por=gc_por,
                cnt=gc_count))


if __name__ == '__main__':
    main()
