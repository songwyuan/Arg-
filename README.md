### Analysis Pipeline (Workflow)
## Phase 1: Data Pre-processing & Quality Control (数据质控与清洗)
## qc.sh & qc2.sh: 自动化提取所有 .fastq.gz 文件执行 FastQC，并通过 MultiQC 汇总测序质量。区分了主测与补测批次。
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

## Phase 3: Pairwise Differential Expression (基础差异分析)
# 3_DE_analysis.R
```
#!/usr/bin/env Rscript
# 加载必要的包
library(DESeq2)
library(ggplot2)

# 1. 设置工作目录和文件路径 (请根据你的实际路径修改)
work_dir <- "/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis"
setwd(work_dir)

count_file <- "04_counts/final_counts.txt"
meta_file <- "04_counts/sample_info.csv"

# 2. 读取数据
counts <- read.table(count_file, header=TRUE, row.names=1, check.names=FALSE)
metadata <- read.csv(meta_file, row.names=1)

# 确保 count 矩阵的列名和 metadata 的行名顺序完全一致
counts <- counts[, rownames(metadata)]
stopifnot(all(colnames(counts) == rownames(metadata)))

# 3. 构建 DESeq2 对象
# 我们使用 True_Group 作为分组进行比较，这样最直接
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ True_Group)

# 过滤低表达基因 (至少在3个样本中 counts >= 10)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# 4. 运行差异表达分析
dds <- DESeq(dds)

# 5. 提取你关心的对比组结果 (示例)
# 提取结果的函数，可以按需调用
get_de_results <- function(dds, group1, group2, prefix) {
  # 注意：对比顺序是 group1 vs group2 (group1 相对于 group2 的变化)
  res <- results(dds, contrast=c("True_Group", group1, group2))
  res_ordered <- res[order(res$padj),]

  # 导出 CSV
  write.csv(as.data.frame(res_ordered),
            file=paste0("04_counts/DE_", prefix, "_", group1, "_vs_", group2, ".csv"))

  # 简单统计显著差异基因 (padj < 0.05 & |log2FC| > 1)
  sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0.58, na.rm=TRUE)
  sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < -0.58 , na.rm=TRUE)
  cat(sprintf("Comparison: %s vs %s -> Up: %d, Down: %d\n", group1, group2, sig_up, sig_down))
}

# --- 在这里定义你需要比较的组合 ---

# 示例 1：在 Arg+ 背景下，LGTV 感染 vs mock
get_de_results(dds, "Arg_plus_LGTV", "Arg_plus_mock", "ArgPlus")

# 示例 2：在 Arg- 背景下，LGTV 感染 vs mock
get_de_results(dds, "Arg_minus_LGTV", "Arg_minus_mock", "ArgMinus")

# 示例 3：LGTV 感染下，Arg+ vs Arg- (看突变效应)
get_de_results(dds, "Arg_plus_LGTV", "Arg_minus_LGTV", "LGTV_infection")

# 对 SFTSV 的比较
get_de_results(dds, "Arg_plus_SFTSV", "Arg_plus_mock", "ArgPlus")
get_de_results(dds, "Arg_minus_SFTSV", "Arg_minus_mock", "ArgMinus")

# 6. PCA 图评估样本聚类 (检查样本是否聚类合理，以及 1号和9号修正后是否和其他重复聚在一起)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("True_Group", "Genotype", "Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("04_counts/PCA_plot.pdf", width=8, height=6)
ggplot(pcaData, aes(PC1, PC2, color=True_Group, shape=Genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot (Corrected Labels)")
dev.off()

cat("🎉 Differential Expression Analysis Setup Complete. Check 04_counts directory for CSVs and PCA plot.\n")
```

## Phase 4: Advanced Likelihood Ratio Test (高级交互建模 LRT) ⭐核心⭐
# 14_Ultimate_LRT_Engine.R
```
#!/usr/bin/env Rscript
library(DESeq2)
library(dplyr)

setwd("/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis/Strategy4_Mine")

cat("🔄 1. 正在读取表达矩阵并注入 Human Symbol...\n")
counts <- read.table("../04_counts/final_counts.txt", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read.csv("../04_counts/sample_info.csv", row.names=1)
counts <- counts[, rownames(metadata)]

universal_mapping <- read.table("../03_counts/tick_to_human_mapping.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

tick_ids <- rownames(counts)
mapped_symbols <- universal_mapping$Human_Symbol[match(tick_ids, universal_mapping$Tick_ID)]
final_names <- ifelse(is.na(mapped_symbols) | mapped_symbols == "", tick_ids, mapped_symbols)

dup_counts <- table(final_names)
dup_names <- names(dup_counts[dup_counts > 1])
idx_dup <- final_names %in% dup_names
final_names[idx_dup] <- paste0(final_names[idx_dup], "_", tick_ids[idx_dup])
rownames(counts) <- final_names

cat("⚙️ 2. 正在构建 DESeq2 LRT 交互模型与分组组合模型...\n")
metadata$Genotype <- gsub("\\+", "_plus", metadata$Genotype)
metadata$Genotype <- gsub("-", "_minus", metadata$Genotype)
metadata$Genotype <- factor(metadata$Genotype, levels = c("Arg_plus", "Arg_minus"))
metadata$Treatment <- factor(metadata$Treatment, levels = c("mock", "LGTV", "SFTSV"))
metadata$Group <- factor(paste(metadata$Genotype, metadata$Treatment, sep="_"))

# LRT 模型 (算交互 P-value)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Genotype + Treatment + Genotype:Treatment)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ Genotype + Treatment)

# Group 模型 (提取特定的 FC)
dds_group <- DESeqDataSetFromMatrix(countData = counts[keep,], colData = metadata, design = ~ Group)
dds_group <- DESeq(dds_group)

cat("🛡️ 3. 正在执行【地狱级双重过滤】：病毒FC > 1.5 AND 交互Delta FC > 1.5 ...\n")
res_lrt <- results(dds_lrt)
sig_lrt_genes <- rownames(subset(res_lrt, padj < 0.05)) # 基础 LRT 显著

# 提取核心应答 FC (FC 1.5 倍的 log2 阈值是 0.585)
fc_ap_lgtv <- results(dds_group, contrast=c("Group", "Arg_plus_LGTV", "Arg_plus_mock"))$log2FoldChange
fc_am_lgtv <- results(dds_group, contrast=c("Group", "Arg_minus_LGTV", "Arg_minus_mock"))$log2FoldChange
fc_ap_sfts <- results(dds_group, contrast=c("Group", "Arg_plus_SFTSV", "Arg_plus_mock"))$log2FoldChange
fc_am_sfts <- results(dds_group, contrast=c("Group", "Arg_minus_SFTSV", "Arg_minus_mock"))$log2FoldChange

fc_matrix <- data.frame(
  row.names = rownames(dds_group),
  AP_LGTV = fc_ap_lgtv, AM_LGTV = fc_am_lgtv, AP_SFTS = fc_ap_sfts, AM_SFTS = fc_am_sfts
)

# ----------------- 核心过滤逻辑构建 -----------------
FC_THRESHOLD <- 0.585  # Log2(1.5) = 0.585

# 卡口 A: 病毒响应大小 (至少在缺 Arg 或有 Arg 的情况下，病毒能引起 1.5 倍变化)
resp_lgtv <- abs(fc_matrix$AP_LGTV) > FC_THRESHOLD | abs(fc_matrix$AM_LGTV) > FC_THRESHOLD
resp_sfts <- abs(fc_matrix$AP_SFTS) > FC_THRESHOLD | abs(fc_matrix$AM_SFTS) > FC_THRESHOLD

# 卡口 B: 交互效应大小 (缺 Arg 导致的病毒响应差距，必须大于 1.5 倍)
delta_lgtv <- abs(fc_matrix$AM_LGTV - fc_matrix$AP_LGTV)
delta_sfts <- abs(fc_matrix$AM_SFTS - fc_matrix$AP_SFTS)

# 开始提纯！
# 条件 1: LRT p-value 必须显著 (前面已取 sig_lrt_genes)
# 条件 2: LGTV 且 SFTSV 都会引起响应 (AND 逻辑)
# 条件 3: Arg 缺失对两种病毒的干预效应都大于 1.5 倍 (AND 逻辑)
strict_logic <- (resp_lgtv & resp_sfts) & (delta_lgtv > FC_THRESHOLD & delta_sfts > FC_THRESHOLD)

# 将逻辑掩码应用到全部基因
valid_genes_overall <- rownames(fc_matrix)[strict_logic]

# 最终取交集：即满足统计学 LRT，又满足绝对 FC
ultimate_genes <- intersect(sig_lrt_genes, valid_genes_overall)

# ----------------------------------------------------

# 如果要求太高导致基因全军覆没，给出退坡建议
if(length(ultimate_genes) == 0) {
  cat("\n⚠️ 警告：地狱级【AND + AND】过滤太过严苛，符合条件的基因数量为 0！\n")
  cat("建议退坡方案：将病毒要求改为 OR，或交互Delta改为 OR。请直接修改脚本第 55 行的 & 为 | 符号。\n")
} else {
  final_res <- res_lrt[ultimate_genes, ]
  final_res <- final_res[order(final_res$padj), ]
  write.csv(as.data.frame(final_res), "Strict_Interaction_Targets.csv")
  saveRDS(dds_lrt, "dds_lrt_model_annotated.rds")
  cat(sprintf("\n🎉 成功提纯！在【最严苛的 AND 逻辑】下，找到了 %d 个绝世罕见的超级核心靶标！\n", length(ultimate_genes)))
}
```

## Trend Clustering & Visualization (趋势聚类与可视化)
# 15_Cluster_and_Plot.R
```
#!/usr/bin/env Rscript
library(DESeq2)
library(tidyverse)
library(ggprism)

setwd("/mnt/alamo01/users/yuansongwei7/IDE8-ARG--LGTV_SFTSV/FJLBFC20260490-01/FJLBFC20260490-01/analysis/Strategy4_Mine")

cat("🔄 1. 加载 LRT 模型并进行 VST 聚类...\n")
dds_lrt <- readRDS("dds_lrt_model_annotated.rds")
# 正确读取 14 号脚本生成的差异基因列表！
lrt_res <- read.csv("Strict_Interaction_Targets.csv", row.names=1)
target_genes <- rownames(lrt_res)

# 执行 VST 转换并聚类 (k=4)
vsd <- vst(dds_lrt, blind=FALSE)
vst_mat <- assay(vsd)[target_genes, ]
set.seed(42)
km <- kmeans(vst_mat - rowMeans(vst_mat), centers = 4)
cluster_info <- data.frame(Gene = names(km$cluster), Cluster = paste0("Cluster_", km$cluster))

cat("🧮 2. 计算每个基因在各组的均值...\n")
norm_counts <- counts(dds_lrt, normalized=TRUE)[target_genes, ]
meta <- as.data.frame(colData(dds_lrt))
meta$Group <- paste(meta$Genotype, meta$Treatment, sep="_")

group_means <- t(apply(norm_counts, 1, function(row) tapply(row, meta$Group, mean)))

cat("📏 3. 计算相对于 Arg+ Mock 的 Log2 Fold Change...\n")
log2_means <- log2(group_means + 1)
baseline <- log2_means[, "Arg_plus_mock"]
log2fc_mat <- log2_means - baseline

cat("🎨 4. 正在绘制族群集体箱线图...\n")
plot_df <- as.data.frame(log2fc_mat) %>%
  rownames_to_column("Gene") %>% left_join(cluster_info, by="Gene") %>%
  pivot_longer(cols = -c(Gene, Cluster), names_to = "Group", values_to = "Log2FC") %>%
  mutate(Genotype = ifelse(grepl("Arg_plus", Group), "Arg+", "Arg-"),
         Treatment = gsub("Arg_plus_|Arg_minus_", "", Group))

plot_df$Treatment <- factor(plot_df$Treatment, levels=c("mock", "LGTV", "SFTSV"))
plot_df$Genotype <- factor(plot_df$Genotype, levels=c("Arg+", "Arg-"))

dodge_width <- 0.7
p <- ggplot(plot_df, aes(x = Treatment, y = Log2FC, color = Genotype, fill = Genotype)) +
  facet_wrap(~ Cluster, scales = "free_y", ncol = 2) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray50", linewidth=0.8) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, position = position_dodge(dodge_width), width=0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width),
              size = 1.8, alpha = 0.8, shape = 21, color = "white", stroke = 0.5) +
  stat_summary(fun = mean, geom = "line", aes(group = Genotype, linetype = Genotype),
               position = position_dodge(dodge_width), linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(dodge_width), size = 3) +
  scale_color_manual(values = c("Arg+" = "#1F78B4", "Arg-" = "#E31A1C")) +
  scale_fill_manual(values = c("Arg+" = "#1F78B4", "Arg-" = "#E31A1C")) +
  theme_prism(base_size = 14) +
  labs(title = "Collective Gene Patterns by Cluster (LRT Targets)",
       subtitle = "1 Point = 1 Gene. Y-axis: Log2FC relative to Arg+ Mock (Baseline = 0)",
       y = "Log2 Fold Change", x = "Viral Infection") +
  theme(legend.position = "bottom", strip.text = element_text(size=13, face="bold"))

pdf("Final_Collective_Cluster_Log2FC.pdf", width=11, height=8)
print(p)
dev.off()
cat("✅ 完美收官！请查看 Final_Collective_Cluster_Log2FC.pdf\n")
```
