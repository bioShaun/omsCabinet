## Reseq

### 参考基因组文件准备

#### split_genome_gtf_by_bed.py

针对单条染色体超过524M的基因组，建议用此脚本将基因组和 gtf 文件拆分，便于后续使用 gatk 等软件分析。

```bash

python split_genome_gtf_by_bed.py \
    --split-bed split.bed \ # 根据此 bed 文件对基因组进行拆分
    --split-fa True \ # 是否拆分基因组文件，默认 True
    --genome-fa genome.fa \ # 基因组序列文件，默认 None， --split-fa True 时，必填
    --split-gtf True \ # 是否拆分 gtf 文件，默认 True
    --genome-gtf genome.gtf(gff) # 基因组 gtf(gff) 文件，默认 None, --split-gtf True 时，必填


# split bed example
chr1A   0   299292406   1AL
chr1A   299292406   604927367   1AS

# output
# genome.splitChr.fa
# genoem.splitChr.gtf(gff)

```