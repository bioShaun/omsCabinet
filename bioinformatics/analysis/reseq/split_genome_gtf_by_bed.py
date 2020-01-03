import fire
import shutil
import delegator
from HTSeq import GFF_Reader
from pathlib import PurePath


class AppNotFound(Exception):
    pass


class FileIsRequired(Exception):
    pass


def gffline(gffreader_item):
    gff_attr = gffreader_item.attr.copy()
    attr_inf = [f'{key}={val}' for key, val
                in gff_attr.items()]
    attr_line = ';'.join(attr_inf)
    basic_inf = gffreader_item.get_gff_line().split('\t')[:-1]
    basic_inf.append(attr_line)
    return '\t'.join(basic_inf)


def check_app(app_name):
    if shutil.which(app_name) is None:
        raise AppNotFound(app_name)


def split_genome_gtf(split_bed, genome_fa=None, genome_gtf=None,
                     split_fa=True, split_gtf=True):

    if split_fa:
        if genome_fa is None:
            raise FileIsRequired(
                'genome fa is required when --split_fa is set!')
        # get split chr fasta
        check_app('bedtools')
        split_chr_fa = PurePath(genome_fa).with_suffix('.splitChr.fa')
        BED2FA_CMD = ('bedtools getfasta '
                      f'-fi {genome_fa} '
                      f'-bed {split_bed} '
                      f'-fo {split_chr_fa} -name')
        delegator.run(BED2FA_CMD)

    if split_gtf:
        if genome_gtf is None:
            raise FileIsRequired(
                'genome gtf is required when --split_gtf is set!')
        # get split chr gtf
        chr_dict = dict()

        with open(split_bed) as bed_inf:
            for eachline in bed_inf:
                eachline_inf = eachline.strip().split()
                chrom = eachline_inf[0]
                start = int(eachline_inf[1])
                end = int(eachline_inf[2])
                split_chr = eachline_inf[3]
                chr_dict.setdefault(chrom, {})[(start, end)] = split_chr

        gtf_suffix = PurePath(genome_gtf).suffix
        split_chr_gtf = PurePath(genome_gtf).with_suffix(f'.splitChr{gtf_suffix}')
        with open(split_chr_gtf, 'w') as gtf_inf:
            for eachline in GFF_Reader(genome_gtf):
                chrom = eachline.iv.chrom
                start = eachline.iv.start
                end = eachline.iv.end
                if chrom in chr_dict:
                    for each_inter in chr_dict[chrom]:
                        if start >= each_inter[0] and end <= each_inter[1]:
                            new_chr = chr_dict[chrom][each_inter]
                            new_start = start - each_inter[0]
                            new_end = end - each_inter[0]
                else:
                    new_chr = chrom
                    new_start = start
                    new_end = end
                eachline.iv.chrom = new_chr
                eachline.iv.start = new_start
                eachline.iv.end = new_end
                if 'gtf' in gtf_suffix:
                    output_line = eachline.get_gff_line().strip() + ';'
                elif 'gff' in gtf_suffix:
                    output_line = gffline(eachline)
                gtf_inf.write(f'{output_line}\n')


if __name__ == '__main__':
    fire.Fire(split_genome_gtf)
