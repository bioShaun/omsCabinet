import attr
import click
import shutil
import delegator
import itertools
import click_completion
import pandas as pd
from Bio import SeqIO
from loguru import logger
from pathlib import Path

click_completion.init()
CURRENT_DIR = Path().cwd()


class BlastFormatError(Exception):
    pass


class BlastNotFound(Exception):
    pass


DB_SUFFIX_DICT = {
    'nucl': set(['.nin', '.nsq', '.nhr']),
    'prot': set(['.pin', '.psq', '.phr']),
}

BLAST_PROGRAM_DICT = {
    ('nucl', 'nucl'): 'blastn',
    ('nucl', 'prot'): 'blastx',
    ('prot', 'prot'): 'blastp',
    ('prot', 'nucl'): 'tblastn',
}

DEFAULT_HIT_FIELDS = [
    'qseqid',
    'sseqid',
    'pident',
    'length',
    'mismatch',
    'gapopen',
    'qstart',
    'qend',
    'sstart',
    'send',
    'evalue',
    'bitscore',
]


def is_prot(seq_record):
    nucl_set = {'A', 'T', 'G', 'C', 'N'}
    seq_set = {seq_i for seq_i in seq_record.seq.upper() if seq_i.isalpha()}
    if seq_set.difference(nucl_set):
        return True
    return False


def filter_by_identity(blasttab, identity):
    blast_df = pd.read_csv(blasttab,
                           sep='\t',
                           header=None,
                           names=DEFAULT_HIT_FIELDS)
    f_blast_df = blast_df[blast_df.pident >= identity]
    f_blast_file = blasttab.with_suffix(f'.pident{identity}.blasttab')
    f_blast_df.to_csv(f_blast_file, sep='\t', header=False, index=False)
    return f_blast_file


def bidirectional_hit(blasttab1, blasttab2, outfile):
    blast1_df = pd.read_csv(blasttab1,
                            sep='\t',
                            header=None,
                            names=DEFAULT_HIT_FIELDS)
    blast2_df = pd.read_csv(blasttab2,
                            sep='\t',
                            header=None,
                            names=DEFAULT_HIT_FIELDS)
    bbh_df = blast1_df.merge(blast2_df,
                             left_on=['qseqid', 'sseqid'],
                             right_on=['sseqid', 'qseqid'])
    bbh_df.to_csv(outfile,
                  sep='\t',
                  index=False,
                  columns=['qseqid_x', 'sseqid_x', 'pident_x', 'pident_y'],
                  header=False)


@attr.s
class Blastobj:

    blast_file = attr.ib()
    blast_path = None
    _db_type = None
    _fa_file = None

    @property
    def _is_fasta(self):
        if self.blast_file.is_file():
            fasta = SeqIO.parse(self.blast_file, "fasta")
            if any(fasta):
                self._fa_file = self.blast_file
                prot_seqs = itertools.takewhile(is_prot,
                                                itertools.islice(fasta, 10))
                if any(prot_seqs):
                    self._db_type = 'prot'
                else:
                    self._db_type = 'nucl'
                return True
        return False

    @property
    def _is_db(self):
        suffix_set = set()
        for file_i in self.blast_file.parent.glob(f'{self.blast_file.name}.*'):
            suffix_set.add(file_i.suffix)
        for dtype in DB_SUFFIX_DICT:
            if len(suffix_set.intersection(DB_SUFFIX_DICT[dtype])) == 3:
                self._db_type = dtype
                return True
        return False

    @property
    def db_type(self):
        if self._db_type is None:
            if self._is_fasta:
                pass
            elif self._is_db:
                pass
            else:
                raise BlastFormatError(
                    f'Unkown Blast format : {self.blast_file}, '
                    'fasta or blast database is alowed.')
        return self._db_type

    @property
    def _make_db(self):
        """
        make blast database
        """
        logger.info(f'Making blast db for {self.fasta}')
        delegator.run((f'{self.blast_path}/makeblastdb '
                       f'-in {self.fasta} -dbtype {self._db_type}'))

    @property
    def db(self):
        if not self._is_db:
            self._make_db
        return self.blast_file

    @property
    def _extract_fasta(self):
        """
        extract fasta from blast database
        """
        logger.info(f'Extracting sequence from db {self.db}')
        delegator.run((f'{self.blast_path}/blastdbcmd '
                       f'-entry all -db {self.db} -out {self.db}.fasta'))
        self._fa_file = Path(f'{self.db}.fasta')

    @property
    def fasta(self):
        if self._is_fasta:
            pass
        elif self._is_db:
            self._extract_fasta
        else:
            raise BlastFormatError(f'Unkown Blast format : {self.blast_file}, '
                                   'fasta or blast database is alowed.')
        return self._fa_file

    def run_blast(self, other, blast_params, outdir):
        blast_program = BLAST_PROGRAM_DICT.get((self.db_type, other.db_type))
        blastout_name = f'{self.fasta.name}.{other.db.name}.blasttab'
        self.blastout_file = outdir / blastout_name
        logger.info(
            f'Running {blast_program} on {self.fasta} against {other.db}')
        delegator.run((f'{self.blast_path}/{blast_program} '
                       f'-query {self.fasta} '
                       f'-db {other.db} '
                       f'-evalue {blast_params.evalue} '
                       f'-outfmt {blast_params.outfmt} '
                       f'-max_target_seqs {blast_params.max_target_seqs} '
                       f'-num_threads {blast_params.num_threads} '
                       f'{blast_params.optionals} '
                       f'-out {self.blastout_file}'))
        if blast_params.perc_identity:
            logger.info(('Filtering blast results percent identity '
                         f'less than {blast_params.perc_identity}'))
            self.blastout_file = filter_by_identity(self.blastout_file,
                                                    blast_params.perc_identity)


@attr.s
class BlastParams:

    evalue = attr.ib()
    outfmt = attr.ib()
    perc_identity = attr.ib()
    max_target_seqs = attr.ib()
    num_threads = attr.ib()
    optionals = attr.ib()


@click.command()
@click.option(
    '-i1',
    '--input1',
    help='first input file.',
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    '-i2',
    '--input2',
    help='second input file, a blast database or a fasta file.',
    type=click.Path(),
    required=True,
)
@click.option(
    '-o',
    '--outdir',
    help='result directory.',
    default=CURRENT_DIR,
)
@click.option(
    '-evalue',
    help='blast evalue cutoff',
    type=click.STRING,
    default='1e-5',
)
@click.option(
    '-outfmt',
    help='blast alignment view options',
    type=click.INT,
    default=6,
)
@click.option(
    '-perc_identity',
    help='Percent identity cutoff',
    type=click.INT,
    default=None,
)
@click.option(
    '-max_target_seqs',
    help=('Maximum number of aligned sequences to keep, '
          'Not applicable for outfmt <= 4.'),
    type=click.INT,
    default=500,
)
@click.option(
    '-num_threads',
    help='Number of threads (CPUs) to use in the BLAST search',
    default=4,
)
@click.option(
    '-optional',
    '--optional_blast_options',
    help='other blast optional parameters.',
    default='',
)
@click.option(
    '-bi',
    '--bidirectional',
    is_flag=True,
)
@click.option(
    '--blast_path',
    help='blast executable path.',
    default='',
)
def main(input1, input2, outdir, evalue, outfmt, perc_identity,
         max_target_seqs, num_threads, optional_blast_options, bidirectional,
         blast_path):
    if blast_path:
        blast_path = Path(blast_path)
    elif shutil.which('blastn'):
        blast_path = Path(shutil.which('blastn')).parent
    else:
        raise BlastNotFound

    Blastobj.blast_path = blast_path

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if bidirectional:
        max_target_seqs = 1

    blast_params = BlastParams(
        evalue,
        outfmt,
        perc_identity,
        max_target_seqs,
        num_threads,
        optional_blast_options,
    )

    input1_blast_obj = Blastobj(blast_file=Path(input1))
    input2_blast_obj = Blastobj(blast_file=Path(input2))
    input1_blast_obj.run_blast(input2_blast_obj, blast_params, outdir)

    if bidirectional:
        input2_blast_obj.run_blast(input1_blast_obj, blast_params, outdir)
        bidirectional_file = input1_blast_obj.blastout_file.with_suffix('.bbh')
        bidirectional_hit(input1_blast_obj.blastout_file,
                          input2_blast_obj.blastout_file, bidirectional_file)


if __name__ == "__main__":
    main()
