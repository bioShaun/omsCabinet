import fire
import pandas as pd
from Bio import SeqIO
import re
import sys

HEDER = [
    'UniqueIdentifier',
    'ClusterName',
    'Members',
    'TaxonName',
    'TaxonIdentifier',
    'RepresentativeMember'
]
PATTERN = '{qid} (.*) n=(\d+) Tax=(.*) TaxID=(\d+) RepID=(\w+)'


def extract_des_from_seq(uniref_fa, outfile):
    des_dict = dict()
    for each_rd in SeqIO.parse(uniref_fa, 'fasta'):
        each_pat = re.compile(PATTERN.format(qid=each_rd.id))
        try:
            des, mem, tax, taxid, repid = each_pat.match(
                each_rd.description).groups()
        except AttributeError:
            print(each_rd.description)
            sys.exit(1)
        des_dict.setdefault('UniqueIdentifier', []).append(each_rd.id)
        des_dict.setdefault('ClusterName', []).append(des)
        des_dict.setdefault('Members', []).append(mem)
        des_dict.setdefault('TaxonName', []).append(tax)
        des_dict.setdefault('TaxonIdentifier', []).append(taxid)
        des_dict.setdefault('RepresentativeMember', []).append(repid)
    des_df = pd.DataFrame(des_dict)
    des_df.to_csv(outfile, sep='\t', columns=HEDER, index=False)


if __name__ == '__main__':
    fire.Fire(extract_des_from_seq)
