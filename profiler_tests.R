## TEST ##

# --- dependencies -------------------------------------------------------------

source("./auxiliary_functions.R")

library(dplyr)
# library(DESeq2)
# library(fgsea)
# library(rjson)

local_path <- "//wsl.localhost/Manjaro/home/FeAR/PROJECTS/transportome_profiler"
 

# --- metasample ---------------------------------------------------------------

# Beforehand:
#   run `kerblam run gen_test_data`
#
# `metasample` is used by the `gen_test_data` pipeline for the generation of the
# reduced dataset that can be used, in turn, for testing purposes. `metasample`
# subsamples the original dataset by rows and columns, keeping the original
# dataset structure, i.e., the relative proportions between the amount of
# samples belonging to the different experimental groups, identified on the
# basis of the `--metavars` parameter (default value: `_primary_site`).
# This test verifies that within-group subsampling actually occurs according to
# the chosen proportion.

file.path(local_path, "src/workflows/gen_test_data.sh") |> readLines() |>
  getCaptured("--metavars (\\w+) ") |> {\(x) paste0("X",x)}() -> metavars
cat("\n\"", metavars, "\" will be used to subsample the large dataset.\n\n",
    sep = "")

file.path(local_path, "src/workflows/gen_test_data.sh") |> readLines() |>
  getCaptured("\"(\\d+)%\"") |> {\(x) as.numeric(x)/100}() -> ss_ratio
cat("\n\"", ss_ratio, "\" detected subsampling factor.\n\n",
    sep = "")

read.csv(file.path(local_path, "data/test_expression_matrix.csv"),
         header = FALSE,
         nrows = 1,
         row.names = 1) |> t() -> sampleset
read.csv(file.path(local_path,
                   "data/expression_matrix_metadata.csv")) -> metadata

# Inner join
metaset <- merge(sampleset, metadata,
                 by.x = "sample", by.y = "sample")

# If 75% reduction => Each group in the reduced dataset should
#                     be (about) 0.25 times the original one 
check(nrow(sampleset)/nrow(metadata), ss_ratio, 1e-3) # Global check
var_levels <- unique(metadata[,metavars])
report_1 <- data.frame(row.names = var_levels)
for (value in var_levels) {
  x <- sum(metaset[,metavars] == value)
  y <- sum(metadata[,metavars] == value)
  report_1[value,1:4] <- c(x, y, round(x/y, digits = 4),
                           ifelse(y > 50, check(x/y, ss_ratio, 1e-2), NA))
}
colnames(report_1) <- c("reduced", "full", "ratio", "check")
print(report_1)



# --- metasplit ----------------------------------------------------------------

# Beforehand:
#   comment the following two lines (86-87) in `select_and_run.py`
#       os.remove(output_dir / f"{key}_case")
#       os.remove(output_dir / f"{key}_control")
#   run `kerblam run heatmaps -l --profile test`
#
# `metasplit` is the program used by `select_and_run.py` to parse the JSON query
# file (`dea_queries.json`) and extract within-group case and control
# submatrices from the global expression matrix for subsequent differential
# expression statistical analysis or ranking metric calculation.
# For each experimental group, this test verifies that both case and control
# submatrices actually contain all and only the expected samples, through a
# hardcoded reimplementation of the queries.

file.path(local_path, "data/in/config/DEA_queries/dea_queries.json") |>
    readLines() -> file_lines

file_lines |> getCaptured("&_primary_site=\\[?(.+?)\\]?&") -> sites
file_lines |> getCaptured("\"(.*?)\": {") -> types

# Show tumor 'type' definitions
names(sites) <- types
print(as.data.frame(sites))

report_2 <- data.frame(row.names = types)
for (cancer in types) {
  
  read.csv(file.path(local_path, "data/deas",
                     paste0(cancer, "_case")), row.names = 1) |>
    colnames() |> gsub("\\.", "-", x=_) -> case_samples
  
  metaset |> filter(
      X_study == "TCGA" &
          X_primary_site %in% strsplit(sites[cancer], ",")[[1]] &
          !(X_sample_type %in% c("Solid Tissue Normal", "Control Analyte"))) |>
      pull(1) -> case_expected
  
  report_2[cancer,1:4] <- c(length(case_samples), length(case_expected),
                            check(length(case_samples), length(case_expected)),
                            all(sort(case_samples) == sort(case_expected)))
  
  read.csv(file.path(local_path, "data/deas",
                     paste0(cancer, "_control")), row.names = 1) |>
    colnames() |> gsub("\\.", "-", x=_) -> ctrl_samples
  
  metaset |> filter(
      X_study %in% c("TCGA", "GTEX") &
          X_primary_site %in% strsplit(sites[cancer], ",")[[1]] &
          X_sample_type %in% c("Normal Tissue", "Solid Tissue Normal")) |>
      pull(1) -> ctrl_expected
  
  report_2[cancer,5:8] <- c(length(ctrl_samples), length(ctrl_expected),
                            check(length(ctrl_samples), length(ctrl_expected)),
                            all(sort(ctrl_samples) == sort(ctrl_expected)))
}
colnames(report_2) <- c("n_case", "expected", "check", "full_match",
                        "n_control", "expected", "check", "full_match")
print(report_2)



# --- generanker ---------------------------------------------------------------

# `generanker` is the program used by `select_and_run.py` to rank genes based on
# their case-control differential expression. For each cancer type, it applies
# one of the implemented metrics to the two submatrices returned by `metasplit`.
# This test performs a spot check to verify that the computed ranking statistic
# actually matches the expected ranking statistic. Note that, depending on the
# chosen metric, the sample tested will be different:
# - if the ranking metric calculation can be applied to the single row of the
#   raw-count expression matrix, then the check is done on 1 random gene from
#   each tumor type;
# - if the ranking metric calculation requires the (pre)processing of the entire
#   expression matrix, the check is done on the first 50 genes from 3 random
#   tumor types.

file.path(local_path, "src/workflows/gen_test_data.sh") |>
    readLines() |> getCaptured(" (\\d+)$") |> as.numeric() -> ngene
cat("\n\"", ngene, "\" rows (genes) detected.\n\n", sep = "")

file.path(local_path, "data/in/config/heatmaps_runtime_options.json") |>
    readLines() |> getCaptured("\"rank_method\": \"(\\w+)\",$") -> metric
cat("\n\"", metric, "\" metric detected.\n\n", sep = "")

# Dynamically call a function based on a string representing its name
get_metric <- get(metric |> sub("norm_", "", x=_))

# If the ranking metric calculation can be applied to the single row of the
# raw-count expression matrix.
if (metric %in% c("fold_change", "s2n_ratio", "cohen_d", "bws_test")) {
  
  report_3 <- data.frame(row.names = types)
  for (cancer in types) {
    
    gene_index <- sample(1:ngene, 1)
    
    read.csv(file.path(local_path, "data/deas",
                       paste0(cancer, "_case")))[gene_index,] |>
      {\(x)x[,-1] |> t() |> as.data.frame() |> `colnames<-`(x[,1])}() -> case
    
    read.csv(file.path(local_path, "data/deas",
                       paste0(cancer, "_control")))[gene_index,] |>
      {\(x)x[,-1] |> t() |> as.data.frame() |> `colnames<-`(x[,1])}() -> control
    
    # Compute the expected statistic
    stat <- get_metric(case, control)
    
    # Extract the actual ranking value
    read.csv(file.path(local_path, "data/deas",
                       paste0(cancer, "_deseq.csv")),
             row.names = 1)[names(case),] -> ranking
    
    report_3[cancer, 1:6] <- c(names(case), names(control),
                               check(names(case), names(control)),
                               ranking, stat, check(ranking, stat, 1e-9))
  }
  colnames(report_3) <- c("case_gene", "control_gene", "check",
                          "rank_stat", "expected", "check")
  print(report_3)
}

# If the ranking metric calculation requires the (pre)processing of the entire
# expression matrix
if (metric %in% c("norm_fold_change", "deseq_shrinkage",
                  "norm_s2n_ratio", "norm_cohen_d", "norm_bws_test")) {
  
  for (cancer in sample(types, 3)) {
  
    case_mat <- read.csv(file.path(local_path, "data/deas",
                                   paste0(cancer, "_case")), row.names = 1)
    ctrl_mat <- read.csv(file.path(local_path, "data/deas",
                                   paste0(cancer, "_control")), row.names = 1)
    dataset <- merge(case_mat, ctrl_mat,
                     by.x = "row.names", by.y = "row.names")
    rownames(dataset) <- dataset$Row.names
    dataset <- dataset[,-which(colnames(dataset) == "Row.names")]
    
    if(!(nrow(case_mat) == nrow(ctrl_mat) && nrow(case_mat) == nrow(dataset))) {
      stop("ERROR: row entry mismatch (check genes)")
    }
    if(!ncol(dataset) == (ncol(ctrl_mat) + ncol(case_mat))) {
      stop("ERROR: column entry mismatch (check samples)")
    }
    
    # Create a DESeq2Dataset object
    condition <- factor(c(rep("case", ncol(case_mat)),
                          rep("control", ncol(ctrl_mat))),
                        levels = c("control", "case"))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = unlog(dataset),
                                          colData = data.frame(condition),
                                          design = ~ condition)
    
    # Median of Ratios (MOR) normalization method
    DESeq2::estimateSizeFactors(dds) |>
      DESeq2::counts(normalized = TRUE) |> relog() -> norm_counts
    
    # View total counts before and after normalization
    data.frame(tot_raw = colSums(unlog(dataset)),
               tot_norm = colSums(unlog(norm_counts)),
               tot_tpm = colSums(DESeq2::fpm(dds, robust=FALSE))) -> counts
    print(cancer)
    #print(counts)
    counts |> filter(stringr::str_detect(rownames(counts), "^TCGA")) |>
        colMeans(na.rm = T) |> `names<-`(paste0("mean_", colnames(counts))) |>
        print()
    counts |> filter(stringr::str_detect(rownames(counts), "^GTEX")) |>
        colMeans(na.rm = T) |> `names<-`(paste0("mean_", colnames(counts))) |>
        print()
    
    if (metric != "deseq_shrinkage") {
      
      # Compute the expected statistic
      stat <- apply(norm_counts[1:50,], 1,
                    \(x)get_metric(x[which(condition == "case")],
                                   x[which(condition == "control")]))
      
    } else if (metric == "deseq_shrinkage") {
        
      # DESeq2 standard analysis
      dds2 <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)
      # res <- DESeq2::results(dds2,
      #                        name = "condition_case_vs_control")   # Standard FCs
      res <- DESeq2::lfcShrink(dds2,
                               coef = "condition_case_vs_control",
                               type = "apeglm",
                               apeAdapt = FALSE)                     # Shrunken FCs
      
      # Extract the expected statistic...
      #stat <- res$log2FoldChange[1:50]
      # ...for the most significant DEGs
      res |> as.data.frame() |> filter(padj < 1e-25) |>
          {\(x) x[1:50, "log2FoldChange", drop = FALSE]}() |> t() |>
          {\(x) x |> as.vector() |> `names<-`(colnames(x))}() -> stat
    }
    
    # Extract the actual ranking value
    read.csv(file.path(local_path, "data/deas",
                       paste0(cancer, "_deseq.csv")),
             row.names = 1)[names(stat), 1, drop = FALSE] -> ranking
    colnames(ranking) <- NULL
    
    # Always compute the norm_fold_change metric as a reference
    ref_logFC <- apply(norm_counts[1:50,], 1,
                       \(x)fold_change(x[which(condition == "case")],
                                       x[which(condition == "control")]))
    
    report_3 <- data.frame(rank_stat = ranking, expected = stat,
                           check = check(ranking[,1], stat, 1e-9),
                           delta = (ranking[,1] - stat),
                           delta_perc = 1e2*(ranking[,1] - stat)/(stat),
                           ref_logFC = ref_logFC)
    print(report_3)
  }
}



# --- run_gsea -----------------------------------------------------------------



file.path(local_path, "data/genesets.json") |>
  rjson::fromJSON(file=_) -> gene_sets



result <- fgsea::fgsea(
  pathways = genesets,
  stats = ranks,
  gseaParam = 1
)




data <- rjson::fromJSON(readr::read_file(file))

# We need a list of id: data
data <- lapply(data, \(x) {x[["data"]]})




