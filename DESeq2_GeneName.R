# ==============================================================================
# 1. 라이브러리 로드
# ==============================================================================
library(DESeq2)
library(SummarizedExperiment)
library(tibble)

# ==============================================================================
# 2. 데이터 불러오기 (.rds 파일)
# ==============================================================================
se <- readRDS("null.merged.gene.SummarizedExperiment.rds")

# ==============================================================================
# 3. 데이터 정제
# ==============================================================================

# 카운트 데이터를 Matrix로 변환하고 Integer로 반올림
# nf-core 결과는 소수점이 포함되어 있으나 DESeq2는 정수 입력만 받음.
counts_data <- assay(se, "salmon.merged.gene_counts_length_scaled")
counts_data <- as.matrix(counts_data)
counts_data <- round(counts_data)

# 메타데이터 생성
# 현재 colData에는 condition 정보가 없으므로 샘플 이름에서 추출
# "Control"이라는 글자가 들어가면 Control, 아니면 Treated로 분류
col_data <- colData(se)
col_data$condition <- ifelse(grepl("Control", col_data$sample), "Control", "Treated")

# 기준 레벨(Reference Level) 설정
# Control을 분석의 기준으로 음수 -> Control 보다 발현량 적음
col_data$condition <- factor(col_data$condition)
col_data$condition <- relevel(col_data$condition, ref = "Control")


# ==============================================================================
# 4. DESeq2 객체 생성 및 분석 실행
# ==============================================================================
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ condition)

# se에 있던 gene_name 정보를 dds에 추가
rowData(dds)$gene_name <- rowData(se)$gene_name

# ==========================================================
# Pre-filtering
# ==========================================================
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)

# ==============================================================================
# 5. 결과 확인 및 저장
# ==============================================================================
res <- results(dds)
res <- as.data.frame(res)

# Ensembl ID는 rownames에 있으므로 칼럼으로 빼기
res$EnsemblID <- rownames(res)

# gene_name 추가
res$gene_name <- rowData(dds)$gene_name

# gene_name을 앞쪽으로 정리
res <- res[, c("gene_name", "EnsemblID", setdiff(colnames(res), c("gene_name", "EnsemblID")))]

# padj 순 정렬
res <- res[order(res$padj), ]

# 저장
write.csv(res, file = "DESeq2_results_with_gene_name.csv", row.names = FALSE)