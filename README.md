### Analysis Pipeline (Workflow) 
## Phase 1: Data Pre-processing & Quality Control (数据质控与清洗) 完全没必要，下次直接用fastp，省去cat和qc
## qc.sh & qc2.sh: 自动化提取所有 .fastq.gz 文件执行 FastQC，并通过 MultiQC 汇总测序质量。区分了主测与补测批次。
# Merge
```
#!/bin/bash
# 1. 定义主测和补测的路径
DIR_MAIN="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/ML150013750"
DIR_TOPUP="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/ML150013680"
OUT_DIR="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis/01_merged_data"

mkdir -p ${OUT_DIR}

echo "正在将 Arg-mock-1 的主测和补测数据进行物理拼接 (cat)..."

# 2. 暴力拼接 R1 文件
cat ${DIR_MAIN}/Arg-mock-1_L1_UDI489.R1.fastq.gz \
    ${DIR_TOPUP}/Arg-mock-1_L1_UDI489.R1.fastq.gz \
    > ${OUT_DIR}/Arg-mock-1_merged.R1.fastq.gz

# 3. 暴力拼接 R2 文件
cat ${DIR_MAIN}/Arg-mock-1_L1_UDI489.R2.fastq.gz \
    ${DIR_TOPUP}/Arg-mock-1_L1_UDI489.R2.fastq.gz \
    > ${OUT_DIR}/Arg-mock-1_merged.R2.fastq.gz

echo "拼接完成！产物为 Arg-mock-1_merged.R1/R2.fastq.gz"
```
# qc.sh
```
#!/bin/bash
set -euo pipefail

# 1. 环境配置 (请根据你的实际 conda 环境路径调整)
export PATH="/mnt/alamo01/users/chenyun730/micromamba/envs/sra-tools/bin:$PATH"

# 2. 路径定义
RAW_DATA_ROOT="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01"
BASE_DIR="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"

mkdir -p "${BASE_DIR}"/{fastqc_results,logs}

# 3. 自动搜寻所有 fastq.gz 文件 (排除 QCFILE 目录)
echo "Scanning for fastq files..."
ALL_FILES=$(find "${RAW_DATA_ROOT}" -name "*.fastq.gz" | grep -v "QCFILE")
SAMPLE_PREFIXES=$(echo "$ALL_FILES" | sed 's/\.R[12]\.fastq\.gz//' | sort | uniq)

# 4. 质控处理函数
process_qc() {
    local prefix="$1"
    local filename=$(basename "$prefix")
    local dirname=$(basename $(dirname "$prefix"))
    # 使用 目录名+文件名 作为唯一 ID，防止补测样本重名覆盖
    local sample_id="${dirname}_${filename}"
    local log_file="${BASE_DIR}/logs/${sample_id}_qc.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START: ${sample_id}" | tee -a "${log_file}"

    fastqc -t 2 \
        "${prefix}.R1.fastq.gz" \
        "${prefix}.R2.fastq.gz" \
        -o "${BASE_DIR}/fastqc_results" \
        >> "${log_file}" 2>&1
}

# 5. 并行执行 (控制同时运行 6 个任务)
MAX_JOBS=6
for prefix in $SAMPLE_PREFIXES; do
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 5
    done
    process_qc "$prefix" &
done

wait

# 6. MultiQC 汇总
echo "Generating MultiQC report..."
multiqc "${BASE_DIR}/fastqc_results" -o "${BASE_DIR}/fastqc_results" --name "IDE8-ARG_QC_Summary"

echo "ALL DONE! Results at: ${BASE_DIR}/fastqc_results/IDE8-ARG_QC_Summary.html"
(rnaseq) yuansongwei7@mgt01:/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/sc
```
# qc2.sh
```
# 定义路径
DIR_MAIN="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/ML150013750"
DIR_TOPUP="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/ML150013680"
COMPARE_DIR="/mnt/alamo01/users/yuansongwei7/analysis/IDE8-ARG/compare_mock1"

mkdir -p $COMPARE_DIR

# 1. 单独对【主测】的 mock-1 跑质控
fastqc ${DIR_MAIN}/Arg-mock-1_L1_UDI489.R1.fastq.gz ${DIR_MAIN}/Arg-mock-1_L1_UDI489.R2.fastq.gz -o $COMPARE_DIR

# 2. 单独对【补测】的 mock-1 跑质控
fastqc ${DIR_TOPUP}/Arg-mock-1_L1_UDI489.R1.fastq.gz ${DIR_TOPUP}/Arg-mock-1_L1_UDI489.R2.fastq.gz -o $COMPARE_DIR
```
# 1_fastp_qc.sh
```
#!/bin/bash
#PBS -N fastp_qc
#PBS -q batch
#PBS -l nodes=1:ppn=24
#PBS -l mem=50gb
#PBS -V
#PBS -o logs/fastp.out
#PBS -e logs/fastp.err

set -euo pipefail

# ⚠️ 注意：因为有了 #PBS -V，我们直接删掉了那些找环境的代码，计算节点会直接使用你现在的 rnaseq 环境！

# 1. 路径定义
BASE_DIR="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
RAW_DATA_DIR="${BASE_DIR}/01_merged_data"
CLEAN_DIR="${BASE_DIR}/02_clean_data"

mkdir -p "${CLEAN_DIR}" "${BASE_DIR}/logs" "${BASE_DIR}/fastp_reports"

# 2. 处理函数
process_fastp() {
    local r1_file="$1"
    local r2_file="${r1_file/.R1./.R2.}"
    local prefix=$(echo "$r1_file" | sed 's/\.R1\.fastq\.gz//')
    local log_file="${BASE_DIR}/logs/${prefix}_fastp.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START fastp: ${prefix}" | tee -a "${log_file}"

    fastp -w 4 \
        -i "${RAW_DATA_DIR}/${r1_file}" -I "${RAW_DATA_DIR}/${r2_file}" \
        -o "${CLEAN_DIR}/${prefix}_clean_R1.fastq.gz" \
        -O "${CLEAN_DIR}/${prefix}_clean_R2.fastq.gz" \
        -h "${BASE_DIR}/fastp_reports/${prefix}_fastp.html" \
        -j "${BASE_DIR}/fastp_reports/${prefix}_fastp.json" \
        >> "${log_file}" 2>&1
}

# 3. 循环处理
cd "${RAW_DATA_DIR}"
echo "找到以下 R1 文件："
ls *.R1.fastq.gz

MAX_JOBS=6
for r1 in *.R1.fastq.gz; do
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 5
    done
    process_fastp "$r1" &
done

wait

# 4. 汇总
echo "Generating MultiQC report..."
multiqc "${BASE_DIR}/fastp_reports" -n "IDE8-ARG_Fastp_Summary" -o "${BASE_DIR}/fastp_reports"

echo "🎉 ALL DONE!"

```

## Phase 2: Host Genome Mapping & Quantification (宿主比对与定量)
# 2_host_mapping.sh
```
#PBS -N STAR_Host_Mapping
#PBS -q batch
#PBS -l nodes=1:ppn=32
#PBS -l mem=100gb
#PBS -V
#PBS -o logs/mapping.out
#PBS -e logs/mapping.err

set -e

# 1. 路径定义
BASE_DIR="/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
DATA_DIR="${BASE_DIR}/02_clean_data"
REF_DIR="${BASE_DIR}/reference"
INDEX_DIR="${REF_DIR}/star_idx_host_only"
OUT_DIR="${BASE_DIR}/03_alignment"
MATRIX_DIR="${BASE_DIR}/04_counts"

mkdir -p ${OUT_DIR} ${MATRIX_DIR}

# 2. 建立索引 (如果你之前建好了，它会自动跳过)
if [ ! -f ${INDEX_DIR}/Genome ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Host Index Construction..."
    STAR --runMode genomeGenerate --runThreadN 32 --genomeDir ${INDEX_DIR} \
         --genomeFastaFiles ${REF_DIR}/tick_genome.fa \
         --sjdbGTFfile ${REF_DIR}/tick_annotation.gtf \
         --sjdbOverhang 149
fi

# 3. 批量比对
cd ${DATA_DIR}
SAMPLES=$(ls *_clean_R1.fastq.gz | sed 's/_clean_R1\.fastq\.gz//' | sort | uniq)

for s in ${SAMPLES}; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Mapping $s ..."
    STAR --runThreadN 24 --genomeDir ${INDEX_DIR} \
         --readFilesIn ${s}_clean_R1.fastq.gz ${s}_clean_R2.fastq.gz \
         --readFilesCommand zcat --outFileNamePrefix ${OUT_DIR}/${s}_ \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 60000000000 \
         --quantMode GeneCounts --outSAMunmapped Within
done

# 4. 提取表达矩阵
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating Count Matrix..."
paste ${OUT_DIR}/*_ReadsPerGene.out.tab | grep -v "N_" | awk '{printf "%s", $1; for(i=4;i<=NF;i+=4) printf "\t%s", $i; printf "\n"}' > ${MATRIX_DIR}/raw_counts_matrix.txt

echo -n "Gene_ID" > ${MATRIX_DIR}/header.txt
ls ${OUT_DIR}/*_ReadsPerGene.out.tab | xargs -n 1 basename | sed 's/_ReadsPerGene.out.tab//' | tr '\n' '\t' | awk '{print "\t" $0}' >> ${MATRIX_DIR}/header.txt
echo "" >> ${MATRIX_DIR}/header.txt

cat ${MATRIX_DIR}/header.txt ${MATRIX_DIR}/raw_counts_matrix.txt > ${MATRIX_DIR}/final_counts.txt

echo "🎉 ALL DONE! Ultimate Matrix saved at ${MATRIX_DIR}/final_counts.txt"
```
## 注释
#  7_annotate_matrix.py
```
import pandas as pd
import re
import os

# ================= 1. 路径配置 =================
BASE_DIR = "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
matrix_file = os.path.join(BASE_DIR, "04_counts/final_counts.txt")
blast_file = os.path.join(BASE_DIR, "reference/tick_vs_swissprot.blast.txt")
gff_file = os.path.join(BASE_DIR, "reference/tick_annotation.gtf")
sample_info_file = os.path.join(BASE_DIR, "04_counts/sample_info.csv")
output_file = os.path.join(BASE_DIR, "04_counts/annotated_comprehensive_counts.txt")

# ================= 2. 解析 BLAST 结果 =================
print("1. 提取 BLAST 详细信息...")
blast_df = pd.read_csv(blast_file, sep='\t', header=None,
                       names=['Protein_ID', 'SwissProt_ID', 'pident', 'length', 'mismatch', 'gapopen',
                              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

def extract_gene_name(sp_id):
    try: return sp_id.split('|')[-1].split('_')[0]
    except: return "-"

def extract_species(sp_id):
    try: return sp_id.split('_')[-1]
    except: return "-"

blast_df['Protein_ID_clean'] = blast_df['Protein_ID'].apply(lambda x: x.split('.')[0])
blast_df['Gene_Symbol'] = blast_df['SwissProt_ID'].apply(extract_gene_name)
blast_df['Species'] = blast_df['SwissProt_ID'].apply(extract_species)
blast_df = blast_df.drop_duplicates(subset=['Protein_ID_clean'], keep='first')
blast_dict = blast_df.set_index('Protein_ID_clean')[['Gene_Symbol', 'Species', 'evalue', 'bitscore']].to_dict('index')

# ================= 3. 解析 GTF =================
print("2. 解析 GTF 文件...")
loc_to_gtf = {}
with open(gff_file, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        if 'CDS' in line and 'protein_id' in line and 'gene_id' in line:
            gene_match = re.search(r'gene_id "([^"]+)"', line)
            prot_match = re.search(r'protein_id "([^"]+)"', line)
            trans_match = re.search(r'transcript_id "([^"]+)"', line)
            if gene_match and prot_match:
                loc_id = gene_match.group(1)
                prot_id_full = prot_match.group(1)
                prot_id_clean = prot_id_full.split('.')[0]
                trans_id = trans_match.group(1) if trans_match else "-"
                loc_to_gtf[loc_id] = {'Transcript_ID': trans_id, 'Protein_ID_full': prot_id_full, 'Protein_ID_clean': prot_id_clean}

# ================= 4. 处理 Sample 映射与重命名 =================
print("3. 处理样本名精简与分组符号化 (+/-)...")
rename_map = {}
final_column_order = []

if os.path.exists(sample_info_file):
    sample_df = pd.read_csv(sample_info_file)
    for _, row in sample_df.iterrows():
        orig_id = row['Sample_ID'].strip()

        # A. 提取精简 ID: 去掉 _L1_UDI 开头的后缀
        # 例如 Arg-mock-2_L1_UDI490 -> Arg-mock-2
        short_id = re.sub(r'_L1_UDI.*', '', orig_id)

        # B. 处理 True_Group: plus -> +, minus -> -
        group_display = row['True_Group'].replace('_plus_', '+_').replace('_minus_', '-_')

        # C. 组合新名字: Group:Short_ID
        new_name = f"{group_display}:{short_id}"

        rename_map[orig_id] = new_name
        final_column_order.append(new_name)
else:
    print(f"错误：未找到 {sample_info_file}")

# ================= 5. 注释并重排矩阵 =================
print("4. 组装并格式化最终矩阵...")
matrix_df = pd.read_csv(matrix_file, sep='\t')
matrix_df.columns = matrix_df.columns.str.strip()
matrix_df = matrix_df.loc[:, ~matrix_df.columns.str.contains('^Unnamed')]

# 插入注释列
t_ids, p_ids, syms, specs, evals, scores = [], [], [], [], [], []
for loc in matrix_df['Gene_ID']:
    gtf_info = loc_to_gtf.get(loc, {})
    p_clean = gtf_info.get('Protein_ID_clean', '-')
    t_ids.append(gtf_info.get('Transcript_ID', '-'))
    p_ids.append(gtf_info.get('Protein_ID_full', '-'))
    b_info = blast_dict.get(p_clean, {})
    syms.append(b_info.get('Gene_Symbol', '-'))
    specs.append(b_info.get('Species', '-'))
    evals.append(b_info.get('evalue', '-'))
    scores.append(b_info.get('bitscore', '-'))

matrix_df.insert(1, 'Transcript_ID', t_ids)
matrix_df.insert(2, 'Protein_ID', p_ids)
matrix_df.insert(3, 'Gene_Symbol', syms)
matrix_df.insert(4, 'Species', specs)
matrix_df.insert(5, 'BLAST_evalue', evals)
matrix_df.insert(6, 'BLAST_bitscore', scores)

# 执行重命名
matrix_df.rename(columns=rename_map, inplace=True)

# 强制排序
base_cols = ['Gene_ID', 'Transcript_ID', 'Protein_ID', 'Gene_Symbol', 'Species', 'BLAST_evalue', 'BLAST_bitscore']
# 检查重命名后的样本是否在矩阵中
valid_samples = [name for name in final_column_order if name in matrix_df.columns]
matrix_df = matrix_df[base_cols + valid_samples]

matrix_df.to_csv(output_file, sep='\t', index=False)
print(f"🎉 任务完成！矩阵已按 '+/-' 格式重命名并排序。")
print(f"示例列名: {valid_samples[0] if valid_samples else 'N/A'}")
```
## Phase 3: Pairwise Differential Expression (基础差异分析)
# 5_annotete_DE_analysis.R
```
#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)

work_dir <- "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
setwd(work_dir)

out_dir <- "05_annotate"
dir.create(out_dir, showWarnings = FALSE)

count_file <- "04_counts/annotated_comprehensive_counts.txt"
meta_file <- "sample_info.csv" 

# =========================================================================
# 2. 读取并拆分数据
# =========================================================================
data <- read.table(count_file, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)

anno <- data[, 1:7]
rownames(anno) <- anno$Gene_ID
counts <- data[, 8:ncol(data)]
rownames(counts) <- data$Gene_ID

# 读取 metadata
metadata <- read.csv(meta_file, row.names=1) 

# 【修复核心 1】：仅修改 metadata 的行名来对齐矩阵，但保持 True_Group 列内容为英文！
short_ids <- gsub("_L1_UDI.*", "", rownames(metadata))
# 模拟 Python 里的显示名字
group_display <- gsub("_plus_", "+_", metadata$True_Group)
group_display <- gsub("_minus_", "-_", group_display)

# 赋予新行名以匹配 count 矩阵
rownames(metadata) <- paste0(group_display, ":", short_ids)

counts <- counts[, rownames(metadata)]
stopifnot(all(colnames(counts) == rownames(metadata)))

# =========================================================================
# 3. 构建 DESeq2 对象
# =========================================================================
# 此时 metadata$True_Group 依然是安全的 "Arg_plus_LGTV"，不会引发报错
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ True_Group)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)

# =========================================================================
# 5. 提取差异结果
# =========================================================================
get_de_results <- function(dds, group1, group2, prefix) {
  res <- results(dds, contrast=c("True_Group", group1, group2))
  res_df <- as.data.frame(res)
  
  res_df <- cbind(anno[rownames(res_df), c("Gene_Symbol", "Species", "BLAST_evalue")], res_df)
  res_ordered <- res_df[order(res_df$padj),]
  
  out_file <- file.path(out_dir, paste0("DE_", prefix, "_", group1, "_vs_", group2, ".csv"))
  write.csv(res_ordered, file=out_file)
  
  sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0.58, na.rm=TRUE)
  sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < -0.58 , na.rm=TRUE)
  cat(sprintf("Comparison: %s vs %s -> Up: %d, Down: %d\n", group1, group2, sig_up, sig_down))
}

# 【修复核心 2】：调用比较组时，用回英文字母
get_de_results(dds, "Arg_plus_LGTV", "Arg_plus_mock", "LGTV_in_ArgPlus")
get_de_results(dds, "Arg_minus_LGTV", "Arg_minus_mock", "LGTV_in_ArgMinus")
get_de_results(dds, "Arg_plus_LGTV", "Arg_minus_LGTV", "ArgPlus_vs_ArgMinus_in_LGTV")

get_de_results(dds, "Arg_plus_SFTSV", "Arg_plus_mock", "SFTSV_in_ArgPlus")
get_de_results(dds, "Arg_minus_SFTSV", "Arg_minus_mock", "SFTSV_in_ArgMinus")
get_de_results(dds, "Arg_plus_SFTSV", "Arg_minus_SFTSV", "ArgPlus_vs_ArgMinus_in_SFTSV")

# =========================================================================
# 6. PCA 降维分析
# =========================================================================
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("True_Group", "Genotype", "Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file.path(out_dir, "PCA_plot.pdf"), width=8, height=6)
print(
  ggplot(pcaData, aes(PC1, PC2, color=True_Group, shape=Genotype)) +
    geom_point(size=4, alpha=0.8) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_bw() + 
    ggtitle("PCA Plot of Arg +/- with LGTV/SFTSV Infection") +
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
)
dev.off()

cat("🎉 差异分析结束！请检查 05_annotate 目录。\n")
```
## KEGG&GO
# enrichmeng.R
```
#!/usr/bin/env Rscript
# =====================================================================
# 非模式生物(蜱虫)借助同源基因(Human)进行 GO & KEGG 自动化富集分析
# =====================================================================

# 安装必要的包 (如果还没装的话)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# 1. 路径设置
work_dir <- "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
setwd(work_dir)

de_dir <- "05_annotate"
out_dir <- "07_enrichment"
dir.create(out_dir, showWarnings = FALSE)

# =====================================================================
# 2. 核心分析函数：输入DE表，输出GO/KEGG结果和全套图表
# =====================================================================
run_enrichment <- function(de_file_name, prefix_name, fc_cutoff=1.5, pvalue_cutoff=0.05) {

  message("\n=======================================================")
  message("🚀 开始处理分析组: ", prefix_name)

  de_file <- file.path(de_dir, de_file_name)
  if (!file.exists(de_file)) {
    message("⚠️ 找不到文件: ", de_file, "，跳过。")
    return(NULL)
  }

  # 读取差异结果并筛选显著 DEG
  res <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)
  sig <- res[!is.na(res$pvalue) & res$pvalue < pvalue_cutoff &
             !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > log2(fc_cutoff), ]

  if (nrow(sig) == 0) {
    message("⚠️ ", prefix_name, " 没有显著差异基因，跳过。")
    return(NULL)
  }

  # 提取合法的 Gene Symbol 并【全部大写】以匹配人类数据库
  symbols <- sig$Gene_Symbol
  symbols <- symbols[!is.na(symbols) & symbols != "" & symbols != "-"]
  symbols <- toupper(symbols)
  symbols <- unique(symbols)

  if (length(symbols) < 5) {
    message("⚠️ ", prefix_name, " 拥有合法 Symbol 的基因太少(<5)，做富集没意义，跳过。")
    return(NULL)
  }

  # ---------------------------------------------------------
  # ID 转换: SYMBOL -> ENTREZID
  # ---------------------------------------------------------
  message("🔄 正在将 Symbol 转换为 Entrez ID...")
  gene_df <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  if (nrow(gene_df) == 0) {
    message("⚠️ ", prefix_name, " 未能成功匹配到任何 Human Entrez ID，跳过。")
    return(NULL)
  }

  target_entrez <- gene_df$ENTREZID

  # =====================================================================
  # 3. GO 富集分析 & 绘图
  # =====================================================================
  message("🧬 正在运行 GO 富集分析...")
  go <- enrichGO(
    gene          = target_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "all",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  if (!is.null(go) && nrow(go@result[go@result$p.adjust < 0.05, ]) > 0) {
    # 保存表格
    write.csv(go@result, file = file.path(out_dir, paste0(prefix_name, "_GO_results.csv")))

    # 画图：为了防止报错，加个简单的 trycatch
    tryCatch({
      # 3.1 柱状图
      p1 <- barplot(go, showCategory = 20, color = "pvalue", label_format = 100) + ggtitle(paste0(prefix_name, " - GO 富集分析柱状图"))
      ggsave(file.path(out_dir, paste0(prefix_name, "_GO_barplot.pdf")), p1, width = 10, height = 8)

      # 3.2 气泡图
      p2 <- dotplot(go, showCategory = 20, color = "pvalue", label_format = 100) + ggtitle(paste0(prefix_name, " - GO 富集分析气泡图"))
      ggsave(file.path(out_dir, paste0(prefix_name, "_GO_dotplot.pdf")), p2, width = 10, height = 8)

      # 3.3 按分类拆分柱状图
      p3 <- barplot(go, drop = TRUE, showCategory = 10, split = "ONTOLOGY", label_format = 100, color = "pvalue") +
            facet_grid(ONTOLOGY~., scale = 'free') + ggtitle(paste0(prefix_name, " - GO 分类柱状图"))
      ggsave(file.path(out_dir, paste0(prefix_name, "_GO_split_barplot.pdf")), p3, width = 12, height = 10)

      # 3.4 按分类拆分气泡图
      p4 <- dotplot(go, showCategory = 10, split = "ONTOLOGY", label_format = 100, color = "pvalue") +
            facet_grid(ONTOLOGY~., scale = 'free') + ggtitle(paste0(prefix_name, " - GO 分类气泡图"))
      ggsave(file.path(out_dir, paste0(prefix_name, "_GO_split_dotplot.pdf")), p4, width = 12, height = 10)

      message("✅ GO 可视化完成！")
    }, error = function(e) { message("⚠️ GO 画图失败: ", e$message) })
  } else {
    message("⚠️ ", prefix_name, " 没有显著富集的 GO 条目。")
  }

  # =====================================================================
  # 4. KEGG 富集分析 & 绘图
  # =====================================================================
  message("🧪 正在运行 KEGG 富集分析...")
  kegg <- enrichKEGG(
    gene          = target_entrez,
    organism      = 'hsa',
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.5, # 注意：KEGG较难富集，你这里设了 0.5 作为阈值，建议如果出图太多可以调回 0.05
    qvalueCutoff  = 0.5
  )

  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
    write.csv(kegg@result, file = file.path(out_dir, paste0(prefix_name, "_KEGG_results.csv")))

    tryCatch({
      # 4.1 KEGG 柱状图
      pk1 <- barplot(kegg, drop = TRUE, showCategory = 15, label_format = 50, color = "pvalue") + ggtitle(paste0(prefix_name, " - KEGG 富集分析柱状图"))
      ggsave(file.path(out_dir, paste0(prefix_name, "_KEGG_barplot.pdf")), pk1, width = 10, height = 8)

      # 4.2 KEGG 气泡图
      pk2 <- dotplot(kegg, showCategory = 20, orderBy = "GeneRatio", label_format = 50, color = "pvalue") + ggtitle(paste0(prefix_name, " - KEGG 富集分析气泡图"))
      ggsave(file.path(out_dir, paste0(prefix_name, "_KEGG_dotplot.pdf")), pk2, width = 10, height = 8)

      message("✅ KEGG 可视化完成！")
    }, error = function(e) { message("⚠️ KEGG 画图失败: ", e$message) })
  } else {
    message("⚠️ ", prefix_name, " 没有显著富集的 KEGG 通路。")
  }
}

# =====================================================================
# 5. 批量执行
# =====================================================================
# 注意这里的 FC 和 P 值过滤标准要和你前一步热图画图时保持一致
run_enrichment("DE_LGTV_in_ArgPlus_Arg_plus_LGTV_vs_Arg_plus_mock.csv", "ArgPlus_LGTV_vs_Mock", fc_cutoff=2)
run_enrichment("DE_SFTSV_in_ArgPlus_Arg_plus_SFTSV_vs_Arg_plus_mock.csv", "ArgPlus_SFTSV_vs_Mock", fc_cutoff=2)
run_enrichment("DE_ArgPlus_vs_ArgMinus_in_LGTV_Arg_plus_LGTV_vs_Arg_minus_LGTV.csv", "LGTV_ArgPlus_vs_ArgMinus", fc_cutoff=1.5)
run_enrichment("DE_ArgPlus_vs_ArgMinus_in_SFTSV_Arg_plus_SFTSV_vs_Arg_minus_SFTSV.csv", "SFTSV_ArgPlus_vs_ArgMinus", fc_cutoff=3)

message("\n🎉 富集分析全部结束！所有的 CSV 表格和精美的 PDF 气泡图/柱状图都保存在了 07_enrichment 目录下！")

```
## venn
# 15_four_way_venn.R
```
#!/usr/bin/env Rscript
# =====================================================================
# 终极经典韦恩图：(四圈独立纯色叠加 + 对称排版 + 无百分比 + 精准四组提取)
# =====================================================================

# 安装经典韦恩图包
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram", repos = "https://cloud.r-project.org")

library(VennDiagram)
library(grid)

# 关闭 VennDiagram 默认的繁琐日志输出
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# 1. 路径设置
work_dir <- "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
setwd(work_dir)

de_dir <- "05_annotate"
out_dir <- "09_four_way_venn"
dir.create(out_dir, showWarnings = FALSE)

# =====================================================================
# 2. 读取四个核心对比组数据
# =====================================================================
message("📦 正在加载四个核心对比组的数据...")

f_lgtv_v  <- file.path(de_dir, "DE_LGTV_in_ArgPlus_Arg_plus_LGTV_vs_Arg_plus_mock.csv")
f_sftsv_v <- file.path(de_dir, "DE_SFTSV_in_ArgPlus_Arg_plus_SFTSV_vs_Arg_plus_mock.csv")
f_lgtv_a  <- file.path(de_dir, "DE_ArgPlus_vs_ArgMinus_in_LGTV_Arg_plus_LGTV_vs_Arg_minus_LGTV.csv")
f_sftsv_a <- file.path(de_dir, "DE_ArgPlus_vs_ArgMinus_in_SFTSV_Arg_plus_SFTSV_vs_Arg_minus_SFTSV.csv")

res_lgtv_v  <- read.csv(f_lgtv_v, row.names = 1, stringsAsFactors = FALSE)
res_sftsv_v <- read.csv(f_sftsv_v, row.names = 1, stringsAsFactors = FALSE)
res_lgtv_a  <- read.csv(f_lgtv_a, row.names = 1, stringsAsFactors = FALSE)
res_sftsv_a <- read.csv(f_sftsv_a, row.names = 1, stringsAsFactors = FALSE)

# =====================================================================
# 3. 提取目标基因集 (全部采用绝对值 abs()，提取所有“差异基因”)
# =====================================================================
message("🎯 正在执行差异基因阈值筛选 (包含上调与下调)...")

genes_lgtv_v <- rownames(res_lgtv_v[!is.na(res_lgtv_v$pvalue) & res_lgtv_v$pvalue < 0.05 &
                                    !is.na(res_lgtv_v$log2FoldChange) & abs(res_lgtv_v$log2FoldChange) > log2(1.5), ])

genes_sftsv_v <- rownames(res_sftsv_v[!is.na(res_sftsv_v$pvalue) & res_sftsv_v$pvalue < 0.05 &
                                      !is.na(res_sftsv_v$log2FoldChange) & abs(res_sftsv_v$log2FoldChange) > log2(1.5), ])

genes_lgtv_a <- rownames(res_lgtv_a[!is.na(res_lgtv_a$pvalue) & res_lgtv_a$pvalue < 0.05 &
                                    !is.na(res_lgtv_a$log2FoldChange) & abs(res_lgtv_a$log2FoldChange) > log2(1.5), ])

genes_sftsv_a <- rownames(res_sftsv_a[!is.na(res_sftsv_a$pvalue) & res_sftsv_a$pvalue < 0.05 &
                                      !is.na(res_sftsv_a$log2FoldChange) & abs(res_sftsv_a$log2FoldChange) > log2(1.5), ])

# 【亮点设计：完美对称排版】
# 按照 1、4 为外圈，2、3 为内圈的绘图逻辑，我们这样排列就能实现左右对称！
venn_list <- list(
  "Arg+(LGTV vs mock)"  = genes_lgtv_v,   # 左外圈
  "LGTV(Arg+ vs Arg-)"    = genes_lgtv_a,   # 左内圈
  "SFTSV(Arg+ vs Arg-)"   = genes_sftsv_v,  # 右内圈
  "Arg+(SFTSV vs mock)" = genes_sftsv_a   # 右外圈
)

# =====================================================================
# 4. 绘制经典四圈韦恩图
# =====================================================================
message("🎨 正在绘制四个大圈独立着色的经典韦恩图...")

# 定义四个大圈的颜色（LGTV用暖色红黄，SFTSV用冷色蓝绿，对比分明）
my_colors <- c("#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A")

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  col = "transparent",          # 隐藏边框线，让颜色更纯粹
  fill = my_colors,             # 四个圈四种独立颜色
  alpha = 0.5,                  # 半透明叠加，中间交集会自动加深
  label.col = "black",          # 数字颜色
  cex = 1.3,                    # 数字大小 (只有count，没有百分比)
  fontface = "bold",
  fontfamily = "sans",
  cat.col = my_colors,          # 标题颜色与圈颜色一致
  cat.cex = 1.1,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  margin = 0.1
)

pdf_file <- file.path(out_dir, "Four_Way_Venn_Classic_Symmetric.pdf")
pdf(pdf_file, width = 9, height = 8)
grid.draw(venn.plot)
dev.off()
message("✅ 经典韦恩图保存成功: ", pdf_file)

# =====================================================================
# 5. 精准提取你要求的 4 个交集
# =====================================================================
extract_detailed_info <- function(target_genes, file_name) {
  if (length(target_genes) == 0) {
    message("⚠️ 交集 [", file_name, "] 无基因，跳过导出。")
    return()
  }

  target_cols <- c("Gene_Symbol", "Species")
  actual_cols <- intersect(target_cols, colnames(res_lgtv_v))
  df <- res_lgtv_v[target_genes, actual_cols, drop = FALSE]

  df$LGTV_Virus_log2FC <- res_lgtv_v[target_genes, "log2FoldChange"]
  df$LGTV_Arg_log2FC   <- res_lgtv_a[target_genes, "log2FoldChange"]
  df$SFTSV_Virus_log2FC <- res_sftsv_v[target_genes, "log2FoldChange"]
  df$SFTSV_Arg_log2FC   <- res_sftsv_a[target_genes, "log2FoldChange"]

  csv_path <- file.path(out_dir, paste0(file_name, ".csv"))
  write.csv(df, csv_path, row.names = TRUE)
  message("  - 💾 成功导出表: ", file_name, ".csv (共 ", length(target_genes), " 个基因)")
}

message("📊 正在严格按照您的 4 个条件提取核心交集数据...")

# 目标 1: LGTV vs MOCK 和 LGTV_Arg 的交集
int_1 <- intersect(genes_lgtv_v, genes_lgtv_a)
extract_detailed_info(int_1, "Target_1_LGTV_Virus_AND_LGTV_Arg")

# 目标 2: SFTSV vs MOCK 和 SFTSV_Arg 的交集
int_2 <- intersect(genes_sftsv_v, genes_sftsv_a)
extract_detailed_info(int_2, "Target_2_SFTSV_Virus_AND_SFTSV_Arg")

# 目标 3: LGTV vs MOCK 和 SFTSV vs MOCK 的交集
int_3 <- intersect(genes_lgtv_v, genes_sftsv_v)
extract_detailed_info(int_3, "Target_3_LGTV_Virus_AND_SFTSV_Virus")

# 目标 4: 四者的终极交集 (即：LGTV_Virus ∩ SFTSV_Virus ∩ LGTV_Arg ∩ SFTSV_Arg)
int_4 <- Reduce(intersect, list(genes_lgtv_v, genes_sftsv_v, genes_lgtv_a, genes_sftsv_a))
extract_detailed_info(int_4, "Target_4_Ultimate_All_Four")

message("\n🎉 逻辑已彻底理顺！请查阅 09_four_way_venn 文件夹！")

```

## 热图
# 4heatmap.R
```
#!/usr/bin/env Rscript
# =====================================================================
# 目标：读取 Venn 交叉基因列表并绘制标准化表达热图
# 升级版：智能三级基因命名 + SCI高级双轨配色 + 自定义列排序
# =====================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap", repos = "https://cloud.r-project.org")
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 package not found. Please install it.")
  library(pheatmap)
  library(DESeq2)
  library(grid)
})

# =====================================================================
# 1. 路径设置
# =====================================================================
analysis_dir <- "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"

count_file <- file.path(analysis_dir, "04_counts/annotated_comprehensive_counts.txt")
meta_file  <- file.path(analysis_dir, "04_counts/sample_info.csv")
gtf_file   <- file.path(analysis_dir, "reference/tick_annotation.gtf")

out_dir <- file.path(analysis_dir, "06_heatmaps/Venn_Targets_Heatmaps")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

venn_dir <- file.path(analysis_dir,"09_four_way_venn")
target_files <- list(
  T1 = file.path(venn_dir,"Target_1_LGTV_Virus_AND_LGTV_Arg.csv"),
  T2 = file.path(venn_dir,"Target_2_SFTSV_Virus_AND_SFTSV_Arg.csv"),
  T3 = file.path(venn_dir,"Target_3_LGTV_Virus_AND_SFTSV_Virus.csv"),
  T4 = file.path(venn_dir,"Target_4_Ultimate_All_Four.csv")
)

# =====================================================================
# 2. 数据预处理核心
# =====================================================================
message("📦 [1/4] 正在读取原始 Counts 和 Metadata 并进行样本名标准化...")
if(!file.exists(count_file)) stop("Count file not found at: ", count_file)
if(!file.exists(meta_file)) stop("Meta file not found at: ", meta_file)

data <- read.table(count_file, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)
counts_raw <- data[, 8:ncol(data)]
rownames(counts_raw) <- data$Gene_ID

metadata <- read.csv(meta_file, stringsAsFactors = FALSE)
metadata$True_Group <- factor(metadata$True_Group)

metadata$Sample_ID <- sub("_L1_UDI.*", "", metadata$Sample_ID)
metadata$Sample_ID <- sub("_merged$", "", metadata$Sample_ID)
metadata$Group_Disp <- gsub("_plus_", "+_",as.character(metadata$True_Group))
metadata$Group_Disp <- gsub("_minus_", "-_", metadata$Group_Disp)

base_colnames <- paste0(metadata$Group_Disp, ":", metadata$Sample_ID)
base_groupnames <- sapply(base_colnames, function(x) {
  g <- strsplit(x, ":")[[1]][1]
  g <- gsub("_", "", g)
  return(g)
})
clean_sample_names <- unname(ave(base_groupnames, base_groupnames, FUN = function(x) paste0(x, ".", seq_along(x))))
rownames(metadata) <- clean_sample_names
colnames(counts_raw) <- rownames(metadata)[match(colnames(counts_raw), base_colnames)]

# --- 重新引入 GTF 解析：提取蛋白功能 (Product) ---
message("🔍 [2/4] 正在扫描 GTF 提取基因功能描述 (用于优化 LOC 注释)...")
if(!file.exists(gtf_file)) stop("GTF file not found at: ", gtf_file)
gtf_lines <- readLines(gtf_file)
gtf_lines <- gtf_lines[grepl("product", gtf_lines)]
locs <- sub('.*gene_id "([^"]+)".*', '\\1', gtf_lines)
prods <- ifelse(grepl('product "', gtf_lines),
                sub('.*product "([^"]+)".*', '\\1', gtf_lines), NA)
loc_to_product <- data.frame(LOC = locs, Product = prods)
loc_to_product <- loc_to_product[!duplicated(loc_to_product$LOC), ]
rownames(loc_to_product) <- loc_to_product$LOC
rm(gtf_lines, locs, prods); gc()

message("🧬 [3/4] 正在进行 DESeq2 VST 转换...")
dds <- DESeqDataSetFromMatrix(countData = counts_raw, colData = metadata, design = ~ True_Group)
vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)

# =====================================================================
# 3. 定义通用绘图函数 (三级智能命名版)
# =====================================================================
plot_gene_list_heatmap <- function(gene_id_list, samples_to_plot, output_prefix, title_suffix) {

  valid_ids <- intersect(gene_id_list, rownames(vst_matrix))
  if(length(valid_ids) < 2) return(NULL)

  valid_samples <- intersect(samples_to_plot, colnames(vst_matrix))
  if(length(valid_samples) == 0) return(NULL)

  expr <- vst_matrix[valid_ids, valid_samples, drop = FALSE]

  # --- 【核心升级：智能三级命名 + 种属标注 + 自动识别 misc_RNA】 ---
  annotation_ref <- data[match(valid_ids, data$Gene_ID), c("Gene_ID", "Gene_Symbol", "Species")]
  rownames(annotation_ref) <- annotation_ref$Gene_ID

  new_rownames <- sapply(valid_ids, function(loc) {
    sym  <- annotation_ref[loc, "Gene_Symbol"]
    spec <- annotation_ref[loc, "Species"]
    prod <- if(loc %in% rownames(loc_to_product)) loc_to_product[loc, "Product"] else ""

    # 1. 处理种属缩写 (例如 Ixodes scapularis -> I. scapularis)
    spec_short <- if(!is.na(spec) && spec != "") {
      parts <- strsplit(spec, " ")[[1]]
      if(length(parts) >= 2) paste0(substr(parts[1], 1, 1), ". ", parts[2]) else spec
    } else { "" }

    # 2. 识别功能描述是否为 misc_RNA 或其他有用类型
    biotype_info <- ""
    if (grepl("misc_RNA|lncRNA|ncRNA|snRNA", prod, ignore.case = TRUE)) {
      biotype_info <- "misc_RNA"
    }

    # --- 判定输出格式 ---
    # 情况 A：有正常的 Symbol (不含 LOC)
    if (!is.na(sym) && sym != "" && sym != "-" && !grepl("^LOC", sym, ignore.case = TRUE)) {
      return(paste0(sym, " (", spec_short, ")"))
    }

    # 情况 B：只有 LOC，尝试通过 Product 提取有用信息
    if (!is.na(prod) && prod != "") {
      # 深度清洗：去掉 uncharacterized LOCxxxx 这种废话
      clean_prod <- gsub("uncharacterized LOC[0-9]+", "", prod, ignore.case = TRUE)
      clean_prod <- gsub("uncharacterized protein", "", clean_prod, ignore.case = TRUE)
      clean_prod <- gsub("isoform X.*", "", clean_prod)
      clean_prod <- gsub(",?\\s*transcript variant.*", "", clean_prod, ignore.case = TRUE)
      clean_prod <- trimws(clean_prod)

      # 如果清洗后有实质意义
      if (clean_prod != "" && !grepl("^LOC", clean_prod) && !grepl("^hypothetical", clean_prod, ignore.case = TRUE)) {
        if (nchar(clean_prod) > 40) clean_prod <- paste0(substr(clean_prod, 1, 37), "...")
        return(paste0(loc, " (", clean_prod, ", ", spec_short, ")"))
      }
    }

    # 情况 C：实在没有功能描述，如果是 misc_RNA，至少标出类型
    if (biotype_info != "") {
      return(paste0(loc, " (", biotype_info, ", ", spec_short, ")"))
    }

    # 情况 D：最后兜底：LOC 编号 + 种属
    return(paste0(loc, " (", spec_short, ")"))
  })

  # ... (前面是 new_rownames <- sapply(...) 的大括号结束位置) ...

  # ===================================================================
  # 【新增】：彻底踢除名为 "LOCxxxx (-)" 的三无废基因
  # ===================================================================
  # grepl 检测是以 LOC 开头、并以 (-) 结尾的字符串
  valid_gene_idx <- !grepl("^LOC[0-9]+\\s*\\(-\\)$", new_rownames)

  expr <- expr[valid_gene_idx, , drop = FALSE]
  new_rownames <- new_rownames[valid_gene_idx]

  n_genes <- nrow(expr)
  if(n_genes < 2) {
    message("⚠️ 剔除无意义 LOC 基因后，剩余基因不足2个，跳过本图。")
    return(NULL)
  }

  rownames(expr) <- make.unique(as.character(new_rownames))

  # ... (下面继续是 强制左右排序 的代码) ...

  rownames(expr) <- make.unique(as.character(new_rownames))
  n_genes <- nrow(expr)

  # --- 强制左右排序 ---
  is_virus <- ifelse(grepl("mock", colnames(expr)), 0, 1)
  is_arg_minus <- ifelse(grepl("Arg\\+", colnames(expr)), 0, 1)
  ordered_cols <- colnames(expr)[order(is_virus, is_arg_minus)]
  expr <- expr[, ordered_cols]

  # --- 高级同色系双轨注释 ---
  annotation_col <- data.frame(
    Virus = ifelse(grepl("mock", colnames(expr)), "Mock",
            ifelse(grepl("SFTSV", colnames(expr)), "SFTSV", "LGTV")),
    Arg   = ifelse(grepl("Arg\\+", colnames(expr)), "Arg+", "Arg-")
  )
  rownames(annotation_col) <- colnames(expr)

  ann_colors <- list(
    Virus = c(SFTSV = "#841F27", LGTV = "#D6604D", Mock = "#2166AC"),
    Arg   = c(`Arg+` = "#F4A582", `Arg-` = "#92C5DE")
  )

  # --- 绘图 ---
  my_color <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  my_breaks <- seq(-2.5, 2.5, length.out = 101)
  plot_height <- min(max(8, n_genes * 0.16), 45)

  ph <- pheatmap(expr,
                 scale = "row",
                 cluster_rows = TRUE, cluster_cols = FALSE,
                 annotation_col = annotation_col, annotation_colors = ann_colors,
                 color = my_color, breaks = my_breaks,
                 show_rownames = TRUE, fontsize_row = 8, fontsize_col = 10,
                 angle_col = 45, cellwidth = 25,
                 main = paste0(title_suffix, "\nZ-score Expression Heatmap (VST scaled)\n | n=", n_genes),
                 silent = TRUE)

  out_pdf <- file.path(out_dir, paste0("VennHeatmap_", output_prefix, ".pdf"))
  pdf(out_pdf, width = 14, height = plot_height)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
  message("  ✅ 成功生成高级热图: ", basename(out_pdf), " (展示 ", n_genes, " 个基因)")
}

# =====================================================================
# 4. 执行绘图 (针对四个 Target 文件)
# =====================================================================
message("🎨 [4/4] 正在读取 Venn 目标文件并绘制热图...")

get_genes_from_csv <- function(file_path) {
  if(!file.exists(file_path)) stop("文件未找到: ", file_path)
  df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  if("Gene_ID" %in% colnames(df)) return(df$Gene_ID) else return(df[[1]])
}

message("\n--- 处理 Target 1 ---")
genes_t1 <- get_genes_from_csv(target_files$T1)
samples_lgtv <- clean_sample_names[grepl("Arg.*(LGTV|mock)", clean_sample_names)]
plot_gene_list_heatmap(genes_t1, samples_lgtv, "T1_LGTV_Intersect", "Target 1 (LGTV response & LGTV Arg-regulated)")

message("\n--- 处理 Target 2 ---")
genes_t2 <- get_genes_from_csv(target_files$T2)
samples_sftsv <- clean_sample_names[grepl("Arg.*(SFTSV|mock)", clean_sample_names)]
plot_gene_list_heatmap(genes_t2, samples_sftsv, "T2_SFTSV_Intersect", "Target 2 (SFTSV response & SFTSV Arg-regulated)")

message("\n--- 处理 Target 3 ---")
genes_t3 <- get_genes_from_csv(target_files$T3)
samples_all <- clean_sample_names[grepl("Arg", clean_sample_names)]
plot_gene_list_heatmap(genes_t3, samples_all, "T3_Common_Virus", "Target 3 (Common LGTV & SFTSV response)")

message("\n--- 处理 Target 4 ---")
genes_t4 <- get_genes_from_csv(target_files$T4)
plot_gene_list_heatmap(genes_t4, samples_all, "T4_Ultimate_Intersect", "Target 4 (Ultimate Intersection of 4 lists)")

message("\n🎉 全部结束！已完美生成极其整洁的高级标准化表达热图。")

```

## 靶点基因箱型图
# 
```
#!/usr/bin/env Rscript
# =====================================================================
# 核心基因表达图 (终极大满贯版：全组别 P 值闭环 + SCI 柳叶刀风格)
# =====================================================================

if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2", repos = "https://cloud.r-project.org")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
if (!requireNamespace("ggsci", quietly = TRUE)) install.packages("ggsci", repos = "https://cloud.r-project.org")

library(DESeq2)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggsci)

# 1. 路径设置
work_dir <- "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
setwd(work_dir)

count_file <- "04_counts/annotated_comprehensive_counts.txt"
meta_file  <- "04_counts/sample_info.csv"
venn_dir   <- "09_four_way_venn"
de_dir     <- "05_annotate"
out_dir    <- "10_expression_plots"
dir.create(out_dir, showWarnings = FALSE)

# =====================================================================
# 2. 读取 Counts 矩阵并进行 VST
# =====================================================================
message("📦 正在加载底层矩阵并进行 VST 归一化...")

data <- read.table(count_file, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)
counts_raw <- data[, 8:ncol(data)]
rownames(counts_raw) <- data$Gene_ID

metadata <- read.csv(meta_file, stringsAsFactors = FALSE)
metadata$True_Group <- factor(metadata$True_Group)

metadata$Sample_ID <- sub("_L1_UDI.*", "", metadata$Sample_ID)
metadata$Sample_ID <- sub("_merged$", "", metadata$Sample_ID)
metadata$Group_Disp <- gsub("_plus_", "+_",as.character(metadata$True_Group))
metadata$Group_Disp <- gsub("_minus_", "-_", metadata$Group_Disp)
rownames(metadata) <- paste0(metadata$Group_Disp, ":", metadata$Sample_ID)
counts_raw <- counts_raw[, rownames(metadata)]

metadata$Arg_Status <- ifelse(grepl("plus", metadata$True_Group), "Arg+", "Arg-")
metadata$Virus_Status <- gsub("Arg_plus_|Arg_minus_", "", metadata$True_Group)
metadata$Virus_Status <- factor(metadata$Virus_Status, levels = c("mock", "LGTV", "SFTSV"), labels = c("Mock", "LGTV", "SFTSV"))

dds <- DESeqDataSetFromMatrix(countData = counts_raw, colData = metadata, design = ~ True_Group)
vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)

# =====================================================================
# 3. 预加载 DESeq2 的原始统计结果
# =====================================================================
message("🔍 正在预加载 DESeq2 原始数据池...")
res_lgtv_v  <- read.csv(file.path(de_dir, "DE_LGTV_in_ArgPlus_Arg_plus_LGTV_vs_Arg_plus_mock.csv"), row.names = 1)
res_sftsv_v <- read.csv(file.path(de_dir, "DE_SFTSV_in_ArgPlus_Arg_plus_SFTSV_vs_Arg_plus_mock.csv"), row.names = 1)
res_lgtv_a  <- read.csv(file.path(de_dir, "DE_ArgPlus_vs_ArgMinus_in_LGTV_Arg_plus_LGTV_vs_Arg_minus_LGTV.csv"), row.names = 1)
res_sftsv_a <- read.csv(file.path(de_dir, "DE_ArgPlus_vs_ArgMinus_in_SFTSV_Arg_plus_SFTSV_vs_Arg_minus_SFTSV.csv"), row.names = 1)

format_p <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("p < 0.001")
  if (p > 0.05) return("ns")
  return(paste0("p = ", signif(p, 2)))
}

# =====================================================================
# 4. 核心绘图函数
# =====================================================================
plot_target_genes <- function(csv_file_name, prefix) {
  target_file <- file.path(venn_dir, csv_file_name)
  if (!file.exists(target_file)) return(NULL)
  
  df_target <- read.csv(target_file, row.names = 1, stringsAsFactors = FALSE)
  genes <- rownames(df_target)
  if (length(genes) == 0) return(NULL)
  
  message("\n🎨 正在为 [", prefix, "] 绘制终极版 SCI 表达图...")
  
  display_names <- sapply(genes, function(g) {
    sym <- df_target[g, "Gene_Symbol"]
    if (!is.na(sym) && sym != "" && sym != "-") return(sym) else return(g)
  })
  
  expr_subset <- vst_matrix[genes, , drop = FALSE]
  rownames(expr_subset) <- display_names
  expr_long <- melt(expr_subset, varnames = c("Gene", "Sample"), value.name = "Expression")
  expr_long$Arg_Status <- metadata[as.character(expr_long$Sample), "Arg_Status"]
  expr_long$Virus_Status <- metadata[as.character(expr_long$Sample), "Virus_Status"]
  
  n_genes <- length(genes)
  n_cols <- min(4, n_genes)
  plot_h <- min(max(6, ceiling(n_genes / 4) * 4.5), 45)
  
  seg_h_list <- list()
  seg_v_list <- list()
  text_list <- list()
  
  for (g in genes) {
    sym <- display_names[[g]]
    y_base <- max(expr_long$Expression[expr_long$Gene == sym], na.rm = TRUE)
    y_min <- min(expr_long$Expression[expr_long$Gene == sym], na.rm = TRUE)
    y_range <- y_base - y_min
    step <- y_range * 0.20  
    tick <- y_range * 0.04
    
    # 提取预先跑好的 DESeq2 p 值
    p_la <- format_p(res_lgtv_a[g, "pvalue"])
    p_sa <- format_p(res_sftsv_a[g, "pvalue"])
    p_lv <- format_p(res_lgtv_v[g, "pvalue"])
    p_sv <- format_p(res_sftsv_v[g, "pvalue"])
    
    # 【重点新增】：实时计算 Mock 组内 Arg+ vs Arg- 的 P 值
    mock_plus <- expr_long$Expression[expr_long$Gene == sym & expr_long$Virus_Status == "Mock" & expr_long$Arg_Status == "Arg+"]
    mock_minus <- expr_long$Expression[expr_long$Gene == sym & expr_long$Virus_Status == "Mock" & expr_long$Arg_Status == "Arg-"]
    p_ma_val <- tryCatch(t.test(mock_plus, mock_minus)$p.value, error = function(e) NA)
    p_ma <- format_p(p_ma_val)
    
    # 构建 5 条精美的连线坐标 (Mock, LGTV, SFTSV 组内连线处于同一高度，互不打架)
    coords <- list(
      list(x1 = 0.8, x2 = 1.2, y = y_base + step * 0.8, label = p_ma), # 【新增】Mock: Arg- vs Arg+
      list(x1 = 1.8, x2 = 2.2, y = y_base + step * 0.8, label = p_la), # LGTV: Arg- vs Arg+
      list(x1 = 2.8, x2 = 3.2, y = y_base + step * 0.8, label = p_sa), # SFTSV: Arg- vs Arg+
      list(x1 = 1.2, x2 = 2.2, y = y_base + step * 2.0, label = p_lv), # Arg+: Mock vs LGTV
      list(x1 = 1.2, x2 = 3.2, y = y_base + step * 3.2, label = p_sv)  # Arg+: Mock vs SFTSV
    )
    
    for (c in coords) {
      seg_h_list[[length(seg_h_list) + 1]] <- data.frame(Gene = sym, x = c$x1, xend = c$x2, y = c$y, yend = c$y)
      seg_v_list[[length(seg_v_list) + 1]] <- data.frame(Gene = sym, x = c$x1, xend = c$x1, y = c$y, yend = c$y - tick)
      seg_v_list[[length(seg_v_list) + 1]] <- data.frame(Gene = sym, x = c$x2, xend = c$x2, y = c$y, yend = c$y - tick)
      text_list[[length(text_list) + 1]] <- data.frame(Gene = sym, x = (c$x1 + c$x2) / 2, y = c$y + tick*2.0, label = c$label)
    }
  }
  
  df_seg_h <- do.call(rbind, seg_h_list)
  df_seg_v <- do.call(rbind, seg_v_list)
  df_text <- do.call(rbind, text_list)
  
  p_box <- ggplot(expr_long, aes(x = Virus_Status, y = Expression, fill = Arg_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 0.8, color = "black", 
                 position = position_dodge(0.8), width = 0.6) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
               size = 1.2, color = "black", alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, 
                 fill = "white", color = "black", stroke = 0.8,
                 position = position_dodge(0.8)) +
    geom_segment(data = df_seg_h, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, size = 0.6) +
    geom_segment(data = df_seg_v, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, size = 0.6) +
    geom_text(data = df_text, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.5, fontface = "italic") +
    facet_wrap(~ Gene, scales = "free_y", ncol = n_cols) +
    scale_fill_lancet() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.45))) +
    theme_classic(base_size = 14) +
    labs(x = "", y = "Normalized Expression (VST)", fill = "Arginine Status") +
    theme(
      strip.text = element_text(face = "bold", size = 13),
      strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.8),
      axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir, paste0(prefix, "_Boxplot_SCI_Style.pdf")), p_box, width = max(8, n_cols * 3.5), height = plot_h)
  message("  ✅ 成功生成全连线 SCI 图谱: ", prefix)
}

# =====================================================================
# 5. 执行批量绘图
# =====================================================================
plot_target_genes("Target_1_LGTV_Virus_AND_LGTV_Arg.csv", "Target_1_LGTV")
plot_target_genes("Target_2_SFTSV_Virus_AND_SFTSV_Arg.csv", "Target_2_SFTSV")
plot_target_genes("Target_3_LGTV_Virus_AND_SFTSV_Virus.csv", "Target_3_Common_Virus")
plot_target_genes("Target_4_Ultimate_All_Four.csv", "Target_4_Ultimate")

message("\n🎉 殿堂级全比较基因表达图完美收官！恭喜完成！")
```
## Phase 4: Advanced Likelihood Ratio Test (高级交互建模 LRT) ⭐核心⭐
# 1LRT.R
```
#!/usr/bin/env Rscript
# =====================================================================
# DESeq2 LRT (似然比检验) 进阶版：条件依赖响应模型 + 效应量双重过滤
# =====================================================================

if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")

library(DESeq2)
library(dplyr)

# 1. 路径设置
work_dir <- "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
setwd(work_dir)

count_file <- "04_counts/annotated_comprehensive_counts.txt"
meta_file  <- "04_counts/sample_info.csv"
out_dir    <- "11_LRT_analysis"
dir.create(out_dir, showWarnings = FALSE)

# =====================================================================
# 2. 读取 Counts 矩阵并重构严谨的 Factor 级别
# =====================================================================
message("📦 正在加载矩阵并重构双因素 Metadata...")

data <- read.table(count_file, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)
counts_raw <- data[, 8:ncol(data)]
rownames(counts_raw) <- data$Gene_ID

gene_anno <- data[, 1:7] 
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)

metadata$Sample_ID <- sub("_L1_UDI.*", "", metadata$Sample_ID)
metadata$Sample_ID <- sub("_merged$", "", metadata$Sample_ID)

# 设定 Baseline (非常重要，影响 coef 的方向)
metadata$Arg_Status <- ifelse(grepl("plus", metadata$True_Group), "ArgPlus", "ArgMinus")
metadata$Arg_Status <- factor(metadata$Arg_Status, levels = c("ArgMinus", "ArgPlus")) 

metadata$Virus_Status <- gsub("Arg_plus_|Arg_minus_", "", metadata$True_Group)
metadata$Virus_Status <- factor(metadata$Virus_Status, levels = c("mock", "LGTV", "SFTSV")) 

metadata$Group_Disp <- gsub("_plus_", "+_", as.character(metadata$True_Group))
metadata$Group_Disp <- gsub("_minus_", "-_", metadata$Group_Disp)
rownames(metadata) <- paste0(metadata$Group_Disp, ":", metadata$Sample_ID)
counts_raw <- counts_raw[, rownames(metadata)]

# =====================================================================
# 3. 运行交互作用 LRT 模型
# =====================================================================
message("🚀 正在运行 Expression ~ Arg + Virus + Arg:Virus 模型...")

dds <- DESeqDataSetFromMatrix(
  countData = counts_raw, 
  colData = metadata, 
  design = ~ Arg_Status + Virus_Status + Arg_Status:Virus_Status
)

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ Arg_Status + Virus_Status)

# 打印所有可用的系数名称，验证你的推断
message("🔍 当前模型包含的系数 (resultsNames):")
print(resultsNames(dds_lrt))

# =====================================================================
# 4. 提取整体显著性与特异性效应量 (核心高阶操作)
# =====================================================================
message("📊 正在提取整体交互 LRT 显著性 & 特异性交互 FoldChange...")

# A. 提取整体 LRT 的显著性 (这决定了基因是否具有总体交互作用)
res_lrt_global <- as.data.frame(results(dds_lrt))
res_lrt_global$Gene_ID <- rownames(res_lrt_global)
# 只保留我们关心的 pvalue 和 padj，避免后面的 log2FoldChange 混淆
res_lrt_global <- res_lrt_global[, c("Gene_ID", "pvalue", "padj")] 
colnames(res_lrt_global) <- c("Gene_ID", "LRT_pvalue", "LRT_padj")

# B. 提取 LGTV 特异性的交互作用系数 (Arg 对 LGTV 响应的改变程度)
res_lgtv_int <- as.data.frame(results(dds_lrt, name = "Arg_StatusArgPlus.Virus_StatusLGTV"))
lgtv_fc <- data.frame(
  Gene_ID = rownames(res_lgtv_int),
  LGTV_Interaction_log2FC = res_lgtv_int$log2FoldChange
)

# C. 提取 SFTSV 特异性的交互作用系数 (Arg 对 SFTSV 响应的改变程度)
res_sftsv_int <- as.data.frame(results(dds_lrt, name = "Arg_StatusArgPlus.Virus_StatusSFTSV"))
sftsv_fc <- data.frame(
  Gene_ID = rownames(res_sftsv_int),
  SFTSV_Interaction_log2FC = res_sftsv_int$log2FoldChange
)

# 组装超级数据框
master_res <- merge(gene_anno, res_lrt_global, by = "Gene_ID", all.y = TRUE)
master_res <- merge(master_res, lgtv_fc, by = "Gene_ID", all.x = TRUE)
master_res <- merge(master_res, sftsv_fc, by = "Gene_ID", all.x = TRUE)

# 按照 LRT padj 排序并保存全集
master_res <- master_res[order(master_res$LRT_padj), ]
write.csv(master_res, file.path(out_dir, "01_LRT_Interaction_Master_Table.csv"), row.names = FALSE)

# =====================================================================
# 5. 执行双重严苛筛选 (统计显著 + 效应量巨大)
# =====================================================================
message("🎯 正在执行终极双重筛选: LRT padj < 0.05 且 |Interaction_FC| > 1 ...")

sig_genes <- master_res %>%
  filter(
    !is.na(LRT_padj) & LRT_padj < 0.05 & 
    (abs(LGTV_Interaction_log2FC) > 1 | abs(SFTSV_Interaction_log2FC) > 1)
  )

write.csv(sig_genes, file.path(out_dir, "02_Arg_Dependent_Antiviral_Genes_Strict.csv"), row.names = FALSE)

message("✅ 分析完成！")
message("💡 结论：共挖掘出 ", nrow(sig_genes), " 个具有核心精氨酸依赖性的抗病毒/促病毒靶点基因。")
message("📁 结果已保存至：", out_dir)
```

## Trend Clustering & Visualization (趋势聚类与可视化)
# 15_Cluster_and_Plot.R
```
#!/usr/bin/env Rscript
# =====================================================================
# 核心靶点基因趋势分族图 (K-means 聚类 + 趋势箱型图，一点一基因)
# =====================================================================

if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2", repos = "https://cloud.r-project.org")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
if (!requireNamespace("ggsci", quietly = TRUE)) install.packages("ggsci", repos = "https://cloud.r-project.org")

library(DESeq2)
library(edgeR)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggsci)

# 1. 路径设置
work_dir <- "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
setwd(work_dir)

count_file <- "04_counts/annotated_comprehensive_counts.txt"
meta_file  <- "04_counts/sample_info.csv"
lrt_file   <- "11_LRT_analysis/02_Arg_Dependent_Antiviral_Genes_Strict.csv"
out_dir    <- "12_Trend_Clustering"
dir.create(out_dir, showWarnings = FALSE)

# =====================================================================
# 2. 读取底层矩阵与 Metadata (计算 log2CPM 用于聚类)
# =====================================================================
message("📦 正在加载底层矩阵与 174 个目标基因...")

# 加载数据
data <- read.table(count_file, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)
counts_raw <- data[, 8:ncol(data)]
rownames(counts_raw) <- data$Gene_ID

metadata <- read.csv(meta_file, stringsAsFactors = FALSE)
metadata$Sample_ID <- sub("_L1_UDI.*", "", metadata$Sample_ID)
metadata$Sample_ID <- sub("_merged$", "", metadata$Sample_ID)
metadata$Group_Disp <- gsub("_plus_", "+_", as.character(metadata$True_Group))
metadata$Group_Disp <- gsub("_minus_", "-_", metadata$Group_Disp)
rownames(metadata) <- paste0(metadata$Group_Disp, ":", metadata$Sample_ID)
counts_raw <- counts_raw[, rownames(metadata)]

# 计算 log2(CPM + 1)，对数转化能让聚类趋势更加平滑稳健
log_cpm_matrix <- cpm(counts_raw, log = TRUE, prior.count = 1)

# 读取之前筛选出的 174 个核心基因
sig_genes_df <- read.csv(lrt_file, stringsAsFactors = FALSE)
target_genes <- sig_genes_df$Gene_ID

if(length(target_genes) == 0) stop("❌ 错误：未读取到显著基因！")

# =====================================================================
# 3. 按组别计算基因表达均值 (实现“一点一基因”)
# =====================================================================
message("🔄 正在聚合生物学重复，提取基因组水平趋势...")

# 获取独立的 6 个实验组
unique_groups <- unique(metadata$True_Group)

# 创建一个均值矩阵：行是 174 个基因，列是 6 个实验组
mean_expr_matrix <- matrix(NA, nrow = length(target_genes), ncol = length(unique_groups))
rownames(mean_expr_matrix) <- target_genes
colnames(mean_expr_matrix) <- unique_groups

for (g in unique_groups) {
  samples_in_group <- rownames(metadata)[metadata$True_Group == g]
  # 取均值
  mean_expr_matrix[, g] <- rowMeans(log_cpm_matrix[target_genes, samples_in_group, drop = FALSE])
}

# =====================================================================
# 4. Z-score 标准化与 K-means 聚类分族
# =====================================================================
message("🧠 正在执行 Z-score 标准化与 K-means 趋势分族...")

# 对每个基因(按行)进行 Z-score 标准化
# 原理：scale() 默认按列标准化，所以我们要先转置 t()，标准化后再转置回来 t()
z_score_matrix <- t(scale(t(mean_expr_matrix)))

# 设置随机种子保证每次跑出来的颜色和聚类编号一致
set.seed(42) 

# K-means 聚类：这里默认分成 4 个族群 (如果跑出来发现可以分得更细，可以把 centers 改为 5 或 6)
k_num <- 4
km_res <- kmeans(z_score_matrix, centers = k_num, nstart = 25)

# 把聚类信息打在基因上
cluster_info <- data.frame(
  Gene_ID = target_genes,
  Trend_Cluster = paste0("Cluster ", km_res$cluster)
)

# 保存聚类结果表 (这表很重要，你要看里面是哪些基因)
final_cluster_table <- merge(sig_genes_df, cluster_info, by = "Gene_ID")
final_cluster_table <- final_cluster_table[order(final_cluster_table$Trend_Cluster), ]
write.csv(final_cluster_table, file.path(out_dir, "Gene_Trend_Clusters.csv"), row.names = FALSE)

# =====================================================================
# 5. 重构画图数据：宽表转长表
# =====================================================================
message("📊 正在绘制 SCI 级趋势分族箱型图...")

z_score_df <- as.data.frame(z_score_matrix)
z_score_df$Gene_ID <- rownames(z_score_df)
z_score_df <- merge(z_score_df, cluster_info, by = "Gene_ID")

plot_data <- melt(z_score_df, id.vars = c("Gene_ID", "Trend_Cluster"), variable.name = "True_Group", value.name = "Z_Score")

# 解析分组信息，方便用 ggplot 画图
plot_data$Arg_Status <- ifelse(grepl("plus", plot_data$True_Group), "Arg+", "Arg-")
plot_data$Arg_Status <- factor(plot_data$Arg_Status, levels = c("Arg-", "Arg+"))

plot_data$Virus_Status <- gsub("Arg_plus_|Arg_minus_", "", plot_data$True_Group)
plot_data$Virus_Status <- factor(plot_data$Virus_Status, levels = c("mock", "LGTV", "SFTSV"), labels = c("Mock", "LGTV", "SFTSV"))

# 计算各个 Cluster 的基因数量，加在图片标题上
cluster_counts <- table(plot_data$Trend_Cluster) / 6 # 除以 6 个组别才是真实基因数
plot_data$Cluster_Label <- sapply(plot_data$Trend_Cluster, function(c) paste0(c, " (n = ", cluster_counts[c], ")"))

# =====================================================================
# 6. 终极绘图：带抖动散点的分组趋势箱型图
# =====================================================================
p_trend <- ggplot(plot_data, aes(x = Virus_Status, y = Z_Score, fill = Arg_Status)) +
  # 绘制箱型图本体 (稍微调低透明度，让散点更明显)
  geom_boxplot(outlier.shape = NA, alpha = 0.6, lwd = 0.6, color = "black", 
               position = position_dodge(0.8), width = 0.6) +
  
  # 叠加散点 (这就是你要的“一个点代表一个基因”，使用了抖动避免点完全重合)
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.5, color = "black", alpha = 0.7, shape = 21, stroke = 0.5, 
             aes(fill = Arg_Status)) + # 给散点也填充一致的颜色
  
  # 叠加一个菱形表示均值
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "white", color = "black", stroke = 1,
               position = position_dodge(0.8)) +
  
  # 按照 Cluster 自动分面展示 4 张图
  facet_wrap(~ Cluster_Label, scales = "free_y", ncol = 2) +
  
  # 用柳叶刀配色
  scale_fill_lancet() +
  
  theme_bw(base_size = 15) +
  labs(x = "Virus Infection", y = "Relative Expression Level (Z-score)", 
       fill = "Arginine Status", title = "Expression Trends of Arg-Dependent Genes") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#e5e5e5", color = "black", linewidth = 1),
    axis.title.y = element_text(face = "bold", margin = margin(r = 12)),
    axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
    axis.text = element_text(color = "black", size = 12),
    panel.border = element_rect(color = "black", linewidth = 1),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  )

# 保存图片
pdf_file <- file.path(out_dir, "Gene_Trends_Boxplot_Clusters.pdf")
ggsave(pdf_file, p_trend, width = 12, height = 10, dpi = 300)

message("✅ 聚类与绘图完美收官！")
message("📁 图表与基因分类结果表保存在：", out_dir)
```
