import fire
import gzip
import pandas as pd
from pathlib import PurePath


def extract_vcf_header(vcf):
    if vcf.suffix == '.gz':
        vcf_inf = gzip.open(vcf)
    else:
        vcf_inf = open(vcf)
    prefix = ''
    for eachline in vcf_inf:
        if vcf.suffix == '.gz':
            eachline = eachline.decode()
        prefix += eachline
        if eachline[:6] == '#CHROM':
            return eachline.strip().split('\t'), prefix


def replace_score(vcf, map_file, outfile):
    vcf = PurePath(vcf)
    vcf_header, vcf_prefix = extract_vcf_header(vcf)
    vcf_df = pd.read_csv(vcf, sep='\t', comment='#',
                         header=None, names=vcf_header)
    map_df = pd.read_csv(
        map_file, sep='\t', header=None,
        names=['chrom', 'start', 'ref', 'alt', 'sample_id', ],
        index_col=0)
    vcf_df.loc[:, 'QUAL'] = [
        map_df.loc[each].sample_id if each in map_df.index
        else 'NA'
        for each in vcf_df.loc[:, 'ID']
    ]
    with gzip.open(outfile, 'wt') as out_inf:
        out_inf.write(vcf_prefix)

    vcf_df.to_csv(outfile, sep='\t', header=False,
                  compression='gzip', mode='a', index=False)


if __name__ == '__main__':
    fire.Fire(replace_score)
