
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

# The BWS test statistic
bws_test <- function(case, ctrl) {
  B <- BWStest::bws_stat(unlist(case), unlist(ctrl))
  return(B)
}

# Dummy function
deseq_shrinkage <- function() {
  return(0)
}




