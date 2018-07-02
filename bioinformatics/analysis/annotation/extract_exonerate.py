import fire


HEADER = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attributes',
    'query_id',
    'identity',
    'similarity',
    'query_length',
    'alignment_length',
]


class Exonerate(object):

    def __init__(self, exonerate_file, chrom=None, start_point=0):
        self.exonerate_file = exonerate_file
        self.chrom = chrom
        self.start_point = start_point
        self.exonerate_gff = dict()
        self.exonerate_inf = dict()

    def extract_gff(self):
        gff_flag = 0
        query_id = None
        with open(self.exonerate_file) as exon_inf:
            for eachline in exon_inf:
                eachline = eachline.lstrip()
                if eachline.startswith('Query:'):
                    query_id = eachline.split(':')[-1].strip()
                elif eachline.startswith('##gff-version 2'):
                    gff_flag = 1
                elif eachline.startswith('# --- END OF GFF DUMP ---'):
                    gff_flag = 0
                elif gff_flag and (not eachline.startswith('#')):
                    eachline_inf = eachline.rstrip('\n').split('\t')
                    start = int(eachline_inf[3])
                    end = int(eachline_inf[4])
                    if self.start_point:
                        eachline_inf[3] = str(start + self.start_point - 1)
                        eachline_inf[4] = str(end + self.start_point - 1)
                    if self.chrom:
                        eachline_inf[0] = self.chrom
                    outline = '\t'.join(eachline_inf)
                    self.exonerate_gff.setdefault(
                        query_id, []).append(outline)

    def extract_ryo(self):
        ryo_flag = 0
        with open(self.exonerate_file) as exon_inf:
            for eachline in exon_inf:
                eachline = eachline.strip()
                if '# --- END OF GFF DUMP ---' in eachline:
                    ryo_flag = 1
                elif ryo_flag and (not eachline.startswith('#')):
                    ryo_inf = eachline.rstrip(
                        '-- completed exonerate analysis')
                    ryo_inf_list = ryo_inf.split()
                    ryo_id = ryo_inf_list[0]
                    ryo_line = '\t'.join(ryo_inf_list)
                    self.exonerate_inf[ryo_id] = ryo_line
                    ryo_flag = 0

    def merge(self, outfile):
        self.extract_gff()
        self.extract_ryo()
        header = '\t'.join(HEADER)
        with open(outfile, 'w') as outfile_inf:
            outfile_inf.write('{hd}\n'.format(hd=header))
            for each_id in self.exonerate_gff:
                each_id_inf = self.exonerate_inf[each_id]
                for eachline in self.exonerate_gff[each_id]:
                    outfile_inf.write('{gff}\t{si}\n'.format(
                        gff=eachline, si=each_id_inf
                    ))


if __name__ == '__main__':
    fire.Fire(Exonerate)
