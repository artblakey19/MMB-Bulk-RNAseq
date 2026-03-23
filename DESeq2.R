# 1. 라이브러리 로드
library(DESeq2)
library(tibble)

# 2. TSV 불러오기
tsv <- read.table("Cebu_salmon.merged.gene_counts_length_scaled.tsv",
                  header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE)
# gene_name 컬럼 제거 후 카운트 행렬만 추출
tsv <- tsv[, colnames(tsv) != "gene_name"]

# 3. 데이터 전처리
# 카운트 데이터를 Matrix로 변환하고 Integer로 반올림
# Salmon 결과는 소수점이 포함되어 있으나 DESeq2는 정수 입력만 받음.
counts_data <- as.matrix(tsv)
counts_data <- round(counts_data)

# 메타데이터 생성
# 샘플 컬럼 이름으로 실험군/대조군 확인
# "Control"이라는 글자가 들어가면 Control, 아니면 Treated로 분류
col_data <- data.frame(
  sample    = colnames(counts_data),
  row.names = colnames(counts_data)
)
col_data$condition <- ifelse(grepl("Control", col_data$sample), "Control", "Treated")

# 기준 레벨(Reference Level) 설정
# Control을 분석의 기준으로 음수 -> Control 보다 발현량 적음
col_data$condition <- factor(col_data$condition)
col_data$condition <- relevel(col_data$condition, ref = "Control")

# 4. DESeq2 객체 생성 및 분석 실행
# 반올림 과정에서 자료형을 Matrix로 바꾸었기 때문에
# DESeqDataSet이 아닌 DESeqDataSetFromMatrix사용
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ condition)
# Pre-filtering
# 최소 3개(가장 작은 그룹 크기) 이상의 샘플에서 10 이상의 count를 가진 유전자만 유지
# 해당 기준은 DESeq2 매뉴얼 문서를 근거로 함
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# DESeq(통계 분석) 실행
dds <- DESeq(dds)

# 5. 결과 확인 및 저장
# 결과 추출
res <- results(dds)
# Gene ID 칼럼 이름 붙이기
res <- rownames_to_column(as.data.frame(res), var = "EnsemblID")

# CSV 파일로 저장
# padj 순으로 정렬
res <- res[order(res$padj), ]
# row.names = FALSE를 하지 않으면 불필요한 row name이 형성됨
write.csv(res, file = "DESeq2_results.csv", row.names = FALSE)