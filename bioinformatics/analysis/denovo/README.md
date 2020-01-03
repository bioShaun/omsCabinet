## Denovo 

存放与基因组 Denovo 相关脚本

### 注释

与基因组注释相关的脚本

#### merge_contig.py

将 fasta 文件和 gtf 文件中的短 contigs 合并

```bash

python merge_contig.py \
    --genome_fa genome.fa \ # 需要合并 contig 的基因组文件
    --contig_list contig.list \ # 需要合并的 contig 列表
    --gtf_file genome.gtf \ # 需要合并 contig 的 gtf 文件，如不需要则省略此参数
    --n_sep 100 \ # contig 之间 N 间隔的数目，默认 100
    --merge_name chrUn # 合并生成的 contig 名称，默认 chrUn

# output
# genome.merge_ctg.fa 合并 contig 后的基因组文件
# genome.ctg.offset.txt 合并前的 contig 在合并的 contig 中的位置
# genome.merge_ctg.gtf 合并 contig 后的 gtf 文件

```

#### add_feature_gff.py

合并通过 TransDecoder 预测编码区域的 GFF 文件和原始 GFF 文件。

```bash

python add_feature_gff.py \
    --raw-gff assembly.gff3 \ # 原始 GFF3 
    --feature-gff transdecode.gff3 \ # TransDecoder 预测 CDS 区域的 GFF3
    --outprefix merge \ # 输出 merge.gff, merge.gtf
    --by gene \ # 使用 transdecode.gff3 中的基因注释替换 assembly.gff3 中的基因注释，无论 assembly.gff3 中该基因是否包含更多转录本信息，选择 tr 模式时，仅替换 assembly.gff3 相同的转录本信息
    --rename True \ # 是否用新的前缀替换默认的 MSTRG， 默认是
    --name-prefix Novel \ # 新前缀， 默认 Novel
    --rm_gene gene.list \ # --by gene 模式下，去除原始 GFF3 中的部分基因

```