
import re
import sys
import copy
import fire
import pandas as pd
from HTSeq import GFF_Reader
from collections import OrderedDict

TR_PATTERN = re.compile(r'(\S+).p\d+$')
CDS_PATTERN = re.compile(r'(\S+).p\d+$')
EXON_PATTERN = re.compile(r'(\S+).p\d+(.exon\w+)')
UTR5_PATTERN = re.compile(r'(\S+).p\d+(.utr5p\w+)')
UTR3_PATTERN = re.compile(r'(\S+).p\d+(.utr3p\w+)')

ID_PATTERN = {
    'mRNA': TR_PATTERN,
    'five_prime_UTR': UTR5_PATTERN,
    'exon': EXON_PATTERN,
    'CDS': CDS_PATTERN,
    'three_prime_UTR': UTR3_PATTERN,
}

ATTR_FILTER = ('fixed', 'tr_id', 'gene_id', 'geneID')


def rename_id(gff_attr, old_pre='MSTRG', new_pre='BMSK10'):
    name_pattern = re.compile(f'(.*){old_pre}.(\d+)(.*)')
    for key in gff_attr:
        val = gff_attr[key]
        if name_pattern.match(str(val)):
            pre, id_idx, sfx = name_pattern.match(val).groups()
            gff_attr[key] = f'{pre}{new_pre}{id_idx:0>6}{sfx}'
    return gff_attr


def gffline(gffreader_item, rename, name_prefix,
            id_filter=ATTR_FILTER):
    gff_attr = gffreader_item.attr.copy()
    if rename:
        gff_attr = rename_id(gff_attr, new_pre=name_prefix)
    attr_inf = [f'{key}={val}' for key, val
                in gff_attr.items()
                if key not in id_filter]
    attr_line = ';'.join(attr_inf)
    basic_inf = gffreader_item.get_gff_line().split('\t')[:-1]
    basic_inf.append(attr_line)
    return '\t'.join(basic_inf)


def gtfline(gff_item, rename, name_prefix):
    gff_attr = gff_item.attr.copy()
    if rename:
        gff_attr = rename_id(gff_attr, new_pre=name_prefix)
    gene_id = gff_attr.get('gene_id')
    if gff_item.type == 'gene':
        gene_id = gff_attr.get('ID')
    tr_id = gff_attr.get('tr_id')
    attr_line = f'gene_id "{gene_id}";'
    if gff_item.type != 'gene':
        attr_line = f'{attr_line} transcript_id "{tr_id}";'
    basic_inf = gff_item.get_gff_line().split('\t')[:-1]
    basic_inf.append(attr_line)
    return '\t'.join(basic_inf)


def update_gene_inf(gene_inf, gff_item):
    if gene_inf is None:
        gene_inf = copy.deepcopy(gff_item)
        attr_dict = {'ID': gff_item.attr['Parent']}
        gene_inf.type = 'gene'
        gene_inf.attr = attr_dict
    else:
        if not gene_inf.attr.get('fixed'):
            if gene_inf.iv.start > gff_item.iv.start:
                gene_inf.iv.start = gff_item.iv.start
            if gene_inf.iv.end < gff_item.iv.end:
                gene_inf.iv.end = gff_item.iv.end
    return gene_inf


def fix_id(input_id, id_type, flag=True):
    if flag:
        id_match = ID_PATTERN.get(id_type).match(input_id)
        return ''.join(id_match.groups())
    else:
        return input_id


def gff2dict(gff, fix_id_flag=False):
    by_gene_dict = OrderedDict()
    by_tr_dict = OrderedDict()
    gene_entry_dict = dict()
    tr2gene = dict()
    for eachline in GFF_Reader(gff):
        if eachline.type == 'gene':
            gene_id = eachline.attr['ID']
            eachline.attr['gene_id'] = gene_id
            gene_entry_dict[gene_id] = eachline
            gene_entry_dict[gene_id].attr['fixed'] = True
            continue
        if 'geneID' in eachline.attr:
            parent = eachline.attr['geneID']
            eachline.attr['Parent'] = parent
        else:
            parent = eachline.attr['Parent']
        if eachline.type in ["transcript", "mRNA"]:
            tr_id = fix_id(eachline.attr['ID'], eachline.type, fix_id_flag)
            eachline.attr['ID'] = tr_id
            gene_id = parent
            tr2gene[tr_id] = parent
        else:
            if 'ID' in eachline.attr:
                eachline.attr['ID'] = fix_id(
                    eachline.attr.get('ID'), eachline.type, fix_id_flag)
            tr_id = fix_id(parent, 'mRNA', fix_id_flag)
            eachline.attr['Parent'] = tr_id
            gene_id = tr2gene[tr_id]
        eachline.attr['tr_id'] = tr_id
        eachline.attr['gene_id'] = gene_id
        by_gene_dict.setdefault(gene_id, []).append(eachline)
        by_tr_dict.setdefault(tr_id, []).append(eachline)
        gene_entry_dict[gene_id] = update_gene_inf(
            gene_entry_dict.get(gene_id), eachline)
    return by_gene_dict, by_tr_dict, gene_entry_dict, tr2gene


def wirte_gff_inf(gff_inf, out_pre, fmt='gff',
                  rename=True, name_prefix='Novel',
                  zero_len=6):
    out_file = f'{out_pre}.{fmt}'
    with open(out_file, 'w') as out_inf:
        for line in gff_inf:
            if fmt == 'gff':
                line_str = gffline(line, rename, name_prefix)
            elif fmt == 'gtf':
                line_str = gtfline(line, rename, name_prefix)
            else:
                raise ValueError(f'Invalid format {fmt}')
            out_inf.write(f'{line_str}\n')


def add_by_tr(raw_tr, f_tr, raw_gene_entry, raw_tr2gene,
              outprefix, rename, name_prefix):
    out_inf_list = []
    for tr_id in raw_tr:
        gene_id = raw_tr2gene[tr_id]
        if gene_id in raw_gene_entry:
            out_inf_list.append(raw_gene_entry.pop(gene_id))
        if tr_id in f_tr:
            out_inf_list.extend(f_tr[tr_id])
        else:
            out_inf_list.extend(raw_tr[tr_id])

    wirte_gff_inf(out_inf_list, outprefix, fmt='gff',
                  rename=rename, name_prefix=name_prefix)
    wirte_gff_inf(out_inf_list, outprefix, fmt='gtf',
                  rename=rename, name_prefix=name_prefix)


def add_by_gene(raw_gene, f_gene, raw_gene_entry,
                f_gene_entry, outprefix, rename,
                name_prefix, rm_gene):
    if rm_gene:
        rm_gene_df = pd.read_csv(rm_gene, sep='\t',
                                 header=None, index_col=0)
    out_inf_list = []
    for gene_id in raw_gene:
        if gene_id in f_gene_entry:
            out_inf_list.append(f_gene_entry[gene_id])
            out_inf_list.extend(f_gene[gene_id])
        else:
            if gene_id not in rm_gene_df.index:
                out_inf_list.append(raw_gene_entry[gene_id])
                out_inf_list.extend(raw_gene[gene_id])
    wirte_gff_inf(out_inf_list, outprefix, fmt='gff',
                  rename=rename, name_prefix=name_prefix)
    wirte_gff_inf(out_inf_list, outprefix, fmt='gtf',
                  rename=rename, name_prefix=name_prefix)


def add_feature_gff(raw_gff, feature_gff, outprefix,
                    by='tr', rename=True, rm_gene=None,
                    name_prefix='Novel'):
    raw_by_gene, raw_by_tr, raw_gene_entry, raw_tr2gene = gff2dict(raw_gff)
    f_by_gene, f_by_tr, f_gene_entry, f_tr2gene = gff2dict(
        feature_gff, fix_id_flag=True)
    if by == 'tr':
        add_by_tr(raw_by_tr, f_by_tr, raw_gene_entry,
                  raw_tr2gene, outprefix, rename,
                  name_prefix)
    else:
        add_by_gene(raw_by_gene, f_by_gene, raw_gene_entry,
                    f_gene_entry, outprefix, rename,
                    name_prefix, rm_gene)


if __name__ == '__main__':
    fire.Fire(add_feature_gff)
