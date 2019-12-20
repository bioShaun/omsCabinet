## Denovo 

存放与基因组 Denovo 相关脚本

### 注释

与基因组注释相关的脚本

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