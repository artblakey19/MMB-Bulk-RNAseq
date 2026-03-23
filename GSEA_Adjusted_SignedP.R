# ============================================================
# GSEA (fgsea)
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(fgsea)
  library(msigdbr)
})

# 1. csv 로드
INPUT_CSV <- "DESeq2_results.csv"
deg_results <- read_csv(INPUT_CSV, show_col_types = FALSE)

# 2. 컬럼명 지정
GENE_COL   <- "EnsemblID"

# 3. 랭킹 벡터 생성 (Signed p-value)
# metric = sign(log2FoldChange) * -log10(pvalue)
# pvalue가 0인 경우에는 1e-300으로 처리
ranks <- deg_results %>%
  filter(!is.na(.data[[GENE_COL]])) %>%
  mutate(
    metric = sign(log2FoldChange) * (-log10(pmax(pvalue, 1e-300)))
  ) %>%
  #is.finites는 NA가 포함된 데이터를 제거하는 역할을 함
  filter(is.finite(metric)) %>%
  arrange(desc(metric)) %>%
  { setNames(.$metric, as.character(.[[GENE_COL]])) }

# 4. MSigDB 데이터 준비
# Subcollection이 없는 Collection은 NA 사용
collections <- tibble::tribble(
  ~label,                 ~collection, ~subcollection,
  "Hallmark",             "H",       NA,
  "C2-CGP",               "C2",      "CGP",
  "C2-CP-Biocarta",       "C2",      "CP:BIOCARTA",
  "C2-CP-KEGG_MEDICUS",   "C2",      "CP:KEGG_MEDICUS",
  "C2-CP-PID",            "C2",      "CP:PID",
  "C2-CP-Reactome",       "C2",      "CP:REACTOME",
  "C2-CP-KEGG_LEGACY",    "C2",      "CP:KEGG_LEGACY",
  "C5-GO-BP",             "C5",      "GO:BP"
)

get_pathways <- function(collection, subcollection) {
  # subcollection이 NA면 인자를 넘기지 않기
  if (is.null(subcollection) || is.na(subcollection) || !nzchar(subcollection)) {
    m <- msigdbr(species = "Homo sapiens", collection = collection)
  } else {
    m <- msigdbr(
      species = "Homo sapiens",
      collection = collection,
      subcollection = subcollection
    )
  }

  split(m$ensembl_gene, m$gs_name)
}

# 5. fgsea 실행 (각 컬렉션별)
set.seed(42)
fgsea_list <- collections %>%
  mutate(
    pathways = purrr::map2(collection, subcollection, get_pathways),
    fgsea_res = purrr::map(
      pathways,
      ~ fgsea(
        pathways = .x,
        stats = ranks,
        minSize = 15,
        maxSize = 500
      )
    )
  )

fgsea_res <- fgsea_list %>%
  transmute(
    res = purrr::map2(
      fgsea_res, label,
      ~ as_tibble(.x) %>% mutate(collection = .y)
    )
  ) %>%
  pull(res) %>%
  bind_rows() %>%
  arrange(padj)

# 6. 결과 저장 (각 컬렉션 별로)
if (!dir.exists("gsea_out_SignedP")) dir.create("gsea_out_SignedP")

fgsea_split <- fgsea_list %>%
  transmute(
    label = label,
    res = purrr::map2(
      fgsea_res, label,
      ~ as_tibble(.x) %>% mutate(collection = .y)
    )
  )

purrr::pwalk(
  fgsea_split,
  function(label, res) {
    out_path <- file.path("gsea_out_SignedP", paste0("gsea_", label, ".csv"))
    write_csv(res %>% select(-leadingEdge), out_path)
  }
)