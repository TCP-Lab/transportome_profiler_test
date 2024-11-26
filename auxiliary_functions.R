
# --- Auxiliary Functions ------------------------------------------------------

# Custom log and log^(-1) functions
unlog <- function(x){round((2^x)-1)}
relog <- function(x){log2(x+1)}

# Check correspondence between actual and expected outcomes
check <- function(x, y, err = 0) {
  if (is.numeric(x) && is.numeric(y)) {
    outcome <- ifelse(abs(x-y) <= err, "OK", "WARNING!!")
  } else if (is.character(x) && is.character(y)) {
    outcome <- ifelse(x == y, "OK", "WARNING!!")
  } else {
    stop("Cannot compare different data type!")
  }
  return(outcome)
}

# Here `pattern` should be a regex with a (lazy) Capturing group
getCaptured <- function(text_lines, pattern) {
  regexec(pattern, text_lines, perl = TRUE) |>
    regmatches(text_lines, m=_) |>
    sapply(\(x) x[2]) |>
    na.omit() |>
    unique() -> values
  return(values)
}

# Make gene sets from MTP-DB tables.
# - Omits NAs to keep only complete cases
# - Remove duplicates
# - Discard gene sets smaller than 'min_set_size' (10)
get_genes_by <- function(db_table, ...) {
    
    min_set_size <- 10
    #min_recurse_set_size <- 40
    
    db_table |> group_by(...) |>
        summarize(elements_dups = list(ensg), .groups = "drop") |>
        as.data.frame() |> na.omit() |>
        mutate(elements = elements_dups |> sapply(unique), .keep = "unused") |>
        filter(elements |> sapply(length) >= min_set_size)
}

# Get the full name of the gene set by magical recursion
get_full_name <- function(gene_set_id, catalogue) {
    
    gene_set_name <- catalogue[[gene_set_id]]$name
    parent_id <- catalogue[[gene_set_id]]$parent
    
    if (is.null(parent_id)) {
        return("whole_transportome")
    } else {
        return(paste0(get_parent(parent_id, catalogue), "///", gene_set_name))
    }
}

# Retrieve gene set from transportome_profiler 'genesets.json' file 
get_gene_set <- function(set_name, full_names_vector, catalogue) {
    full_names_vector[grep(set_name, full_names_vector)] |> names() -> set_id
    catalogue[[set_id]]$data
}


# Escape cation charges and adjust 'pore_forming' value to match names in
# 'genesets.json' file 
esc_value <- function(value) {
    value |>
        gsub("+", "\\+", x=_, fixed = TRUE) |>
        gsub("0", "0.0", x=_, fixed = TRUE) |>
        gsub("1", "1.0", x=_, fixed = TRUE)
}

# Function factory: returns functions to compare actual gene sets with expected
# ones
make_checker <- function(base_category) {
    function(feature, value, summary_table) {
        
        value <- as.character(value)
        
        paste0("^", base_category, "///", feature, "::", esc_value(value), "$") |>
            get_gene_set(full_names, gene_sets) -> x
        
        summary_table[summary_table[,feature] == value, "elements"] |>
            unlist() -> y
        
        return(setequal(x,y))
    }
}

# Get a random feature-value pair tuple from a given table
get_random_pairs <- function(summary_table) {
    index <- sample(1:nrow(summary_table),1)
    random_pairs <- summary_table[index, 1:ncol(summary_table)-1]
    print(random_pairs)
    return(random_pairs)
}








# All the following functions take in two numeric vectors or single-column
# matrices or data frames.

# log_2(FC)
fold_change <- function(case, ctrl) {
  FC <- mean(unlist(case)) - mean(unlist(ctrl))
  return(FC)
}

# S2N ratio
s2n_ratio <- function(case, ctrl) {
  signal <- fold_change(case, ctrl)
  noise <- sqrt(var(unlist(case))) + sqrt(var(unlist(ctrl)))
  s2n <- signal/noise
  if ((abs(s2n) == Inf) || is.nan(s2n)) {
    s2n <- 0
  }
  return(s2n)
}

# Pooled Standard Deviation for two independent samples (Cohen definition)
pooled_SD <- function(case, ctrl) {
  case |> unlist() |> length() -> n_case
  ctrl |> unlist() |> length() -> n_ctrl
  
  (((n_case - 1)*var(unlist(case)) + (n_ctrl - 1)*var(unlist(ctrl)))/
    (n_case + n_ctrl - 2)) |> sqrt() |> as.numeric() -> pSD
  return(pSD)
}

# Cohen's d
cohen_d <- function(case, ctrl) {
  d <- fold_change(case, ctrl)/pooled_SD(case, ctrl)
  if ((abs(d) == Inf) || is.nan(d)) {
    d <- 0
  }
  return(d)
}

# This function implements the computation of the BWS statistic for two
# independent samples as defined in
#
#   W. Baumgartner, P. Weiß and H. Schindler 'A Nonparametric Test for the
#   General Two-Sample Problem'. Biometrics, Vol. 54, No. 3 (Sep., 1998),
#   pp. 1129-1135.
#
# along with its 'one-sided' modification introduced by Neuhauser in 2001:
#
#   Neuhauser M. 'One-sided two-sample and trend tests based on a modified
#   Baumgartner-Weig-Schindler statistic'. J. Nonparam. Statist. 13,729-739.
#
# Notably, this function uses the option 'ties.method = "max"' for ranking to
# be as consistent as possible with the definition of B given by the authors in
# original the paper.
#
#   <<  ...the rank Gi (Hj) of each element Xi (Yj) is defined as the number of
#       data in both sets smaller or equal to Xi (Yj). >>
#
# So, for 'alternative = "two-sided"' this function returns the exact same
# values returned by the following function from the BWStest package:
#
#   BWStest::bws_stat(unlist(case), unlist(ctrl))
#
# On the contrary, SciPy's 'bws_test' function output can be reproduced by
# changing the setting to 'ties.method = "average"'.
#
bws_test <- function(case, ctrl, alternative = "one-sided") {
    
    ctrl |> unlist() |> length() -> n
    case |> unlist() |> length() -> m
    
    # Compute the ranks of the combined samples...
    combined_rank <- rank(c(unlist(ctrl), unlist(case)), ties.method = "max")
    # ...and reassign to case and control
    combined_rank[1:n] |> sort() -> Ri
    combined_rank[n+1:m] |> sort() -> Hj
    
    i <- seq(1, n)
    j <- seq(1, m)
    
    Bx_num <- Ri - ((m + n) / n) * i
    By_num <- Hj - ((m + n) / m) * j
    
    if (alternative == "two-sided") {
        Bx_num <- Bx_num^2
        By_num <- By_num^2
    } else if (alternative == "one-sided") {
        Bx_num <- Bx_num * abs(Bx_num)
        By_num <- By_num * abs(By_num)
    }
    
    Bx_den <- (i / (n + 1)) * (1 - i / (n + 1)) * m * ((m + n) / n)
    By_den <- (j / (m + 1)) * (1 - j / (m + 1)) * n * ((m + n) / m)
    
    Bx <- (1/n) * sum(Bx_num / Bx_den, na.rm = TRUE)
    By <- (1/m) * sum(By_num / By_den, na.rm = TRUE)
    
    if (alternative == "two-sided") {
        B <- (Bx + By) / 2
    } else if (alternative == "one-sided") {
        B <- (By - Bx) / 2
    }
    
    return(B)
}




