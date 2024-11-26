## TEST ##

# --- dependencies -------------------------------------------------------------

source("./auxiliary_functions.R")

library(dplyr)
# library(DESeq2)
# library(fgsea)
# library(rjson)
library(DBI)
library(RSQLite)


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

# Each group in the reduced dataset should be (about) 'ss_ratio' times the
# original one. Skip check for groups smaller than 50.
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



# --- MTP-DB -------------------------------------------------------------------

mtpdb <- file.path(local_path, "data/MTPDB.sqlite")
if (Sys.info()["sysname"] == "Windows") {
    # Can't access the DB directly on WSL...
    mtpdb <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "MTPDB.sqlite")
    file.copy(from = file.path(local_path, "data/MTPDB.sqlite"), to = mtpdb)
}

connection <- dbConnect(SQLite(), dbname = mtpdb)

dbListTables(connection) |>
    sapply(\(x) connection |> dbReadTable(x) |> duplicated() |> sum()) |>
    as.data.frame() |> `colnames<-`("duplicated") |>
    mutate(check = ifelse(duplicated == 0, "OK", "WARNING!!"))

dbDisconnect(connection)



# --- make_genesets  -----------------------------------------------------------

# Beforehand:
#   generate the maximum amount of gene sets (still larger than 10) by setting:
#    - '--no_prune' option in `make_genesets.py` (in the 'heatmaps.makefile')
#    - 'min_pop_score=0' in the 'generate_gene_list_trees' function
#    - 'min_recurse_set_size=0' in the 'generate_gene_list_trees' function
#
#   run `kerblam run heatmaps -l --profile test`

file.path(local_path, "data/genesets.json") |>
    rjson::fromJSON(file=_) -> gene_sets

gene_sets |> names() |> sapply(get_full_name, gene_sets) -> full_names

connection <- dbConnect(SQLite(), dbname = mtpdb)

# Get channel table
connection |> dbReadTable("channels") |>
    select(ensg, gating_mechanism, carried_solute) -> channel_table
connection |> dbReadTable("structure") |>
    select(ensg, pore_loops) |>
    mutate(pore_forming = ifelse(pore_loops > 0, 1, 0)) |>
    select(-pore_loops) -> pore_table
merge(channel_table, pore_table, by = "ensg", all.x = TRUE) -> channel_table
# Group and summarize
channel_table |> get_genes_by(gating_mechanism) -> channels_by_G
channel_table |> get_genes_by(carried_solute) -> channels_by_S
channel_table |> get_genes_by(pore_forming) -> channels_by_P
channel_table |> get_genes_by(gating_mechanism, carried_solute) -> channels_by_GS
channel_table |> get_genes_by(gating_mechanism, pore_forming) -> channels_by_GP
channel_table |> get_genes_by(carried_solute, pore_forming) -> channels_by_SP
channel_table |> get_genes_by(gating_mechanism, carried_solute, pore_forming) -> channels_by_GSP

base_cat <- "whole_transportome///pores///channels"

# Check ancestor (level 0)
setequal(get_gene_set(paste0("^", base_cat ,"$"), full_names, gene_sets),
         channel_table$ensg)

# Check level 1
check_channels <- make_checker(base_cat)

channels_by_G$gating_mechanism |>
    sapply(\(x) check_channels("gating_mechanism", x, channels_by_G))

channels_by_S$carried_solute |>
    sapply(\(x) check_channels("carried_solute", x, channels_by_S))

channels_by_P$pore_forming |> as.character() |>
    sapply(\(x) check_channels("pore_forming", x, channels_by_P))

# Check level 2
combo <- get_random_pairs(channels_by_GS)
setequal(
    paste0("^", base_cat,
           "///gating_mechanism::", combo$gating_mechanism,
           "///carried_solute::", esc_value(combo$carried_solute),
           "$") |>
        get_gene_set(full_names, gene_sets),
    channels_by_GS |>
        filter(gating_mechanism == combo$gating_mechanism &
                   carried_solute == combo$carried_solute) |>
    select(elements) |> unlist())

combo <- get_random_pairs(channels_by_GP)
setequal(
    paste0("^", base_cat,
           "///gating_mechanism::", combo$gating_mechanism,
           "///pore_forming::", esc_value(combo$pore_forming),
           "$") |>
        get_gene_set(full_names, gene_sets),
    channels_by_GP |>
        filter(gating_mechanism == combo$gating_mechanism &
                   pore_forming == combo$pore_forming) |>
        select(elements) |> unlist())

combo <- get_random_pairs(channels_by_SP)
setequal(
    paste0("^", base_cat,
           "///carried_solute::", esc_value(combo$carried_solute),
           "///pore_forming::", esc_value(combo$pore_forming),
           "$") |>
        get_gene_set(full_names, gene_sets),
    channels_by_SP |>
        filter(carried_solute == combo$carried_solute &
                   pore_forming == combo$pore_forming) |>
        select(elements) |> unlist())

# Check level 3 #### This always fails.... #### This always fails.... #### This always fails....
combo <- get_random_pairs(channels_by_GSP)
setequal(
    paste0("^", base_cat,
           "///gating_mechanism::", combo$gating_mechanism,
           "///carried_solute::", esc_value(combo$carried_solute),
           "///pore_forming::", esc_value(combo$pore_forming),
           "$") |>
        get_gene_set(full_names, gene_sets),
    channels_by_GSP |>
        filter(gating_mechanism == combo$gating_mechanism &
                   carried_solute == combo$carried_solute &
                   pore_forming == combo$pore_forming) |>
        select(elements) |> unlist())



# Get AQP table --- OK
connection |> dbReadTable("aquaporins") -> AQP_table
base_cat <- "whole_transportome///pores///aquaporins"
setequal(get_gene_set(paste0("^", base_cat ,"$"), full_names, gene_sets),
         AQP_table$ensg)



# Get SLC table
connection |> dbReadTable("solute_carriers") |>
    select(ensg, carried_solute, port_type, direction) -> SLC_table
# Group and summarize
SLC_table |> get_genes_by(carried_solute) -> SLCs_by_S
SLC_table |> get_genes_by(port_type) -> SLCs_by_T
SLC_table |> get_genes_by(direction) -> SLCs_by_D
SLC_table |> get_genes_by(carried_solute, port_type) -> SLCs_by_ST
SLC_table |> get_genes_by(carried_solute, direction) -> SLCs_by_SD
SLC_table |> get_genes_by(port_type, direction) -> SLCs_by_TD
SLC_table |> get_genes_by(carried_solute, port_type, direction) -> SLCs_by_STD

base_cat <- "whole_transportome///transporters///solute_carriers"

# Check ancestor (level 0)
setequal(get_gene_set(paste0("^", base_cat ,"$"), full_names, gene_sets),
         SLC_table$ensg)

# Check level 1
check_SLCs <- make_checker(base_cat)

SLCs_by_S$carried_solute |>
    sapply(\(x) check_SLCs("carried_solute", x, SLCs_by_S))

SLCs_by_T$port_type |>
    sapply(\(x) check_SLCs("port_type", x, SLCs_by_T))

SLCs_by_D$direction |>
    sapply(\(x) check_SLCs("direction", x, SLCs_by_D))

# Check level 2
combo <- get_random_pairs(SLCs_by_ST)
setequal(
    paste0("^", base_cat,
           "///carried_solute::", esc_value(combo$carried_solute),
           "///port_type::", combo$port_type,
           "$") |>
        get_gene_set(full_names, gene_sets),
    SLCs_by_ST |>
        filter(carried_solute == combo$carried_solute &
                   port_type == combo$port_type) |>
        select(elements) |> unlist())

combo <- get_random_pairs(SLCs_by_SD)
setequal(
    paste0("^", base_cat,
           "///carried_solute::", esc_value(combo$carried_solute),
           "///direction::", combo$direction,
           "$") |>
        get_gene_set(full_names, gene_sets),
    SLCs_by_SD |>
        filter(carried_solute == combo$carried_solute &
                   direction == combo$direction) |>
        select(elements) |> unlist())

combo <- get_random_pairs(SLCs_by_TD)
setequal(
    paste0("^", base_cat,
           "///port_type::", combo$port_type,
           "///direction::", combo$direction,
           "$") |>
        get_gene_set(full_names, gene_sets),
    SLCs_by_TD |>
        filter(port_type == combo$port_type &
                   direction == combo$direction) |>
        select(elements) |> unlist())

# Check level 3 #### This always fails.... #### This always fails.... #### This always fails....
combo <- get_random_pairs(SLCs_by_STD)
setequal(
    paste0("^", base_cat,
           "///carried_solute::", esc_value(combo$carried_solute),
           "///port_type::", combo$port_type,
           "///direction::", combo$direction,
           "$") |>
        get_gene_set(full_names, gene_sets),
    SLCs_by_STD |>
        filter(carried_solute == combo$carried_solute &
                   port_type == combo$port_type &
                   direction == combo$direction) |>
        select(elements) |> unlist())



# Get pumps table --- OK
connection |> dbReadTable("pumps") |>
    select(ensg, carried_solute, direction) -> pump_table
# Group and summarize
pump_table |> get_genes_by(carried_solute) -> pumps_by_S
pump_table |> get_genes_by(direction) -> pumps_by_D
pump_table |> get_genes_by(carried_solute, direction) -> pumps_by_SD

base_cat <- "whole_transportome///transporters///atp_driven///pumps"

# Check ancestor (level 0)
setequal(get_gene_set(paste0("^", base_cat ,"$"), full_names, gene_sets),
         pump_table$ensg)

# Check level 1
check_pumps <- make_checker(base_cat)

pumps_by_S$carried_solute |>
    sapply(\(x) check_pumps("carried_solute", x, pumps_by_S))

pumps_by_D$direction |>
    sapply(\(x) check_pumps("direction", x, pumps_by_D))

# Check level 2
combo <- get_random_pairs(pumps_by_SD)
setequal(
    paste0("^", base_cat,
           "///carried_solute::", esc_value(combo$carried_solute),
           "///direction::", combo$direction,
           "$") |>
        get_gene_set(full_names, gene_sets),
    pumps_by_SD |>
        filter(carried_solute == combo$carried_solute &
                   direction == combo$direction) |>
        select(elements) |> unlist())



# Get ABC table --- OK
connection |> dbReadTable("ABC_transporters") |>
    select(ensg, carried_solute, direction) -> ABC_table
# Group and summarize
ABC_table |> get_genes_by(carried_solute) -> ABCs_by_S
ABC_table |> get_genes_by(direction) -> ABCs_by_D
ABC_table |> get_genes_by(carried_solute, direction) -> ABCs_by_SD # Empty

base_cat <- "whole_transportome///transporters///atp_driven///ABC"

# Check ancestor (level 0)
setequal(get_gene_set(paste0("^", base_cat ,"$"), full_names, gene_sets),
         ABC_table$ensg)

# Check level 1
check_ABCs <- make_checker(base_cat)

ABCs_by_S$carried_solute |>
    sapply(\(x) check_ABCs("carried_solute", x, ABCs_by_S))

ABCs_by_D$direction |>
    sapply(\(x) check_ABCs("direction", x, ABCs_by_D))







dbDisconnect(connection)
















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




