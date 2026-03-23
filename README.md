# RNA-seq Differential Expression & GSEA Pipeline

nf-core/rnaseq (Salmon) 결과물을 입력으로 받아 **DESeq2 차등 발현 분석**과 **fgsea 기반 GSEA**를 수행하는 R 스크립트 모음입니다.

---

## 목차

- [입력 데이터](#입력-데이터)
- [스크립트 설명](#스크립트-설명)
- [실행 순서](#실행-순서)
- [환경 재현 (renv)](#환경-재현-renv)

---

## 입력 데이터

| 파일 | 설명 |
|------|------|
| `salmon.merged.gene_counts_length_scaled.tsv` | nf-core/rnaseq 파이프라인이 생성한 Salmon gene-level 정량 결과 (TSV). 라이브러리 크기와 평균 전사체 길이를 모두 보정한 count matrix |

> **참고:** TSV 파일의 첫 번째 열(`gene_id`)은 Ensembl ID, 두 번째 열(`gene_name`)은 HGNC Gene symbol이며, 나머지 열이 샘플별 count 값입니다.

---

## 스크립트 설명

### 1. `DESeq2.R` — 기본 차등 발현 분석

Salmon 정량 결과(TSV)를 로드하여 DESeq2 차등 발현 분석을 수행합니다.

- **조건 분류:** 샘플명에 `Control` 포함 여부로 Control / Treated 자동 분류
- **Pre-filtering:** 최소 3개 샘플에서 count ≥ 10인 유전자만 유지 ([DESeq2 매뉴얼 권장 기준](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html))
- **출력:** `DESeq2_results.csv` (padj 오름차순 정렬, EnsemblID 기반)

### 2. `DESeq2_GeneName.R` — Gene Name 포함 차등 발현 분석

`DESeq2.R`과 동일한 분석을 수행하되, **결과에 gene_name(유전자 심볼) 칼럼을 추가**합니다.

- TSV 파일의 `gene_name` 칼럼을 별도 보관한 뒤, Ensembl ID를 key로 DESeq2 객체에 매핑
- **출력:** `DESeq2_results_with_gene_name.csv` (gene_name, EnsemblID 순으로 칼럼 배치)
- `gene_name`은 HGNC Gene symbol로 여러 EnsemblID에 같은 symbol이 매핑될 수 있어 사용 시 중복 처리 필요
- EnsemblID에 대응하는 Gene symbol이 없을 시 EnsemblID가 gene_name으로 사용되니 유의

### 3. `GSEA_Adjusted_Wald.R` — GSEA (Wald 통계량 기반 Ranking)

DESeq2 결과를 입력으로 **fgsea 기반 Gene Set Enrichment Analysis**를 수행합니다.

- **Ranking metric:** DESeq2의 `stat` (Wald statistic) 칼럼
- **유전자 식별자:** Ensembl ID
- **유전자 세트:** MSigDB (msigdbr 패키지) — Hallmark, C2(CGP, CP 하위 컬렉션), C5(GO:BP)
- **출력:** `gsea_out_wald/` 디렉토리에 컬렉션별 CSV 파일

### 4. `GSEA_Adjusted_SignedP.R` — GSEA (Signed p-value 기반 Ranking)

`GSEA_Adjusted_Wald.R`과 동일한 분석이지만, **다른 Ranking metric**을 사용합니다.

- **Ranking metric:** `sign(log2FoldChange) × −log10(pvalue)` (Signed p-value)
  - p-value가 0인 경우 `1e-300`으로 대체하여 −log10 변환 시 오류 방지
- **출력:** `gsea_out_SignedP/` 디렉토리에 컬렉션별 CSV 파일

> **Wald vs Signed p-value:**  
> Wald 통계량은 fold change와 분산을 동시에 반영하며, Signed p-value는 통계적 유의성에 더 가중치를 둡니다.  
> 두 방법의 결과를 비교하여 보다 견고한 해석을 할 수 있습니다.

---

## 실행 순서

```text
1. DESeq2.R 또는 DESeq2_GeneName.R   →   DESeq2 결과 CSV 생성
2. GSEA_Adjusted_Wald.R              →   Wald 기반 GSEA 결과 생성
3. GSEA_Adjusted_SignedP.R           →   Signed p-value 기반 GSEA 결과 생성
```

> GSEA 스크립트는 `DESeq2_results.csv`를 입력으로 사용하므로, 반드시 DESeq2 스크립트를 먼저 실행해야 합니다.

---

## 환경 재현 (renv)

이 프로젝트는 [renv](https://rstudio.github.io/renv/)를 사용하여 R 패키지 버전을 관리합니다.

### 처음 클론했을 때 (환경 복원)

```r
# 1. renv 패키지 설치 (아직 없는 경우)
install.packages("renv")

# 2. renv.lock에 기록된 패키지 버전 그대로 복원
renv::restore()
```

`renv::restore()`를 실행하면 `renv.lock` 파일에 명시된 패키지들이 프로젝트 로컬 라이브러리(`renv/library/`)에 설치됩니다.

### 프로젝트 관리자용 (renv 초기 설정)

아직 renv이 초기화되지 않은 경우, 아래 순서로 설정합니다.

```r
# 1. renv 초기화 — renv/ 폴더와 .Rprofile이 자동 생성됨
renv::init()

# 2. Bioconductor 패키지 사용을 위한 설정
#    DESeq2, fgsea 등은 Bioconductor에서 설치
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "fgsea"))

# 3. CRAN 패키지 설치
install.packages(c("tibble", "readr", "dplyr", "purrr", "msigdbr"))

# 4. 현재 설치된 패키지 상태를 renv.lock에 스냅샷
renv::snapshot()
```

### 패키지 추가 / 업데이트 후

```r
# 새 패키지 설치 후 lock 파일 갱신
renv::snapshot()

# 변경 사항을 GitHub에 push
# renv.lock 파일은 반드시 커밋에 포함해야 합니다
```

### .gitignore 권장 설정

```gitignore
# renv
renv/library/
renv/staging/
renv/cellar/
renv/lock/
renv/python/
renv/sandbox/

# R
.Rhistory
.RData
.Rproj.user/

# 데이터 파일 (용량이 큰 경우)
*.rds
```

> **핵심 파일:** `renv.lock`과 `renv/activate.R`은 반드시 Git에 포함되어야 합니다.  
> 이 두 파일만 있으면 누구든 `renv::restore()`로 동일한 환경을 구성할 수 있습니다.

---

## 필요 패키지

| 패키지 | 출처 | 용도 |
|--------|------|------|
| `DESeq2` | Bioconductor | 차등 발현 분석 |
| `fgsea` | Bioconductor | GSEA 실행 |
| `tibble` | CRAN | 데이터프레임 유틸리티 |
| `readr` | CRAN | CSV 읽기/쓰기 |
| `dplyr` | CRAN | 데이터 조작 |
| `purrr` | CRAN | 함수형 매핑 |
| `msigdbr` | CRAN | MSigDB 유전자 세트 |
