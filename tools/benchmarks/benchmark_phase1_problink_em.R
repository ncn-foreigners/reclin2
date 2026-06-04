#!/usr/bin/env Rscript

script_path <- function() {
  args <- commandArgs(FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[1])))
  normalizePath(sys.frames()[[1]]$ofile)
}

git_out <- function(repo, args) {
  system2("git", c("-C", repo, args), stdout = TRUE, stderr = TRUE)
}

repo_root <- function() {
  normalizePath(git_out(getwd(), c("rev-parse", "--show-toplevel"))[1])
}

source_reclin2_r <- function(repo) {
  files <- c(
    "R/extend_to.R",
    "R/set_seed.R",
    "R/problink_em.R",
    "R/problink_em_gamma.R",
    "R/problink_em_mixed.R"
  )
  for (f in files) source(file.path(repo, f), chdir = TRUE)
}

mk_binary_patterns <- function() {
  K <- 12L
  lev <- rep(list(c(FALSE, TRUE)), K)
  patterns <- do.call(expand.grid, c(lev, KEEP.OUT.ATTRS = FALSE))
  names(patterns) <- paste0("b", seq_len(K))
  set.seed(1001)
  signal <- rowSums(patterns)
  patterns$n <- as.integer(1 + rpois(nrow(patterns), lambda = 5 + signal^2))
  form <- as.formula(paste("~", paste(names(patterns)[seq_len(K)], collapse = "+")))
  list(formula = form, patterns = patterns)
}

mk_gamma_matrix <- function(N = 5000L, K = 4L, hurdle = TRUE) {
  set.seed(if (hurdle) 1002 else 1003)
  match <- rbinom(N, 1, 0.05) == 1
  X <- matrix(NA_real_, N, K)
  for (k in seq_len(K)) {
    if (hurdle) {
      X[match, k] <- ifelse(rbinom(sum(match), 1, 0.8) == 1, 0,
                            rgamma(sum(match), shape = 2, scale = 0.25))
      X[!match, k] <- ifelse(rbinom(sum(!match), 1, 0.002) == 1, 0,
                             rgamma(sum(!match), shape = 2, scale = 3))
    } else {
      X[match, k] <- rgamma(sum(match), shape = 2, scale = 0.35) + 1e-8
      X[!match, k] <- rgamma(sum(!match), shape = 2, scale = 3) + 1e-8
    }
  }
  colnames(X) <- paste0("g", seq_len(K))
  X
}

mk_mixed_matrix <- function(N = 5000L) {
  set.seed(1004)
  match <- rbinom(N, 1, 0.05) == 1
  X <- matrix(0, N, 6L)
  colnames(X) <- c("b1", "b2", "c1", "c2", "x1", "x2")
  for (k in 1:2) {
    X[match, k] <- rbinom(sum(match), 1, 0.92)
    X[!match, k] <- rbinom(sum(!match), 1, 0.08)
  }
  for (k in 3:4) {
    X[match, k] <- sample(0:4, sum(match), TRUE, prob = c(0.84, rep(0.04, 4)))
    X[!match, k] <- sample(0:4, sum(!match), TRUE, prob = c(0.08, rep(0.23, 4)))
  }
  for (k in 5:6) {
    X[match, k] <- ifelse(rbinom(sum(match), 1, 0.75) == 1, 0,
                          rgamma(sum(match), shape = 2, scale = 0.25))
    X[!match, k] <- rgamma(sum(!match), shape = 2, scale = 3)
  }
  X
}

time_case <- function(name, fun, reps, warmup) {
  for (i in seq_len(warmup)) invisible(suppressWarnings(suppressMessages(fun())))
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    gc()
    times[i] <- system.time(invisible(suppressWarnings(suppressMessages(fun()))))[["elapsed"]]
  }
  data.frame(
    case = name,
    median_elapsed = median(times),
    min_elapsed = min(times),
    max_elapsed = max(times),
    reps = reps,
    stringsAsFactors = FALSE
  )
}

run_benchmarks <- function(repo, reps, warmup) {
  source_reclin2_r(repo)

  binary <- mk_binary_patterns()
  gamma_hurdle <- mk_gamma_matrix(hurdle = TRUE)
  gamma_pure <- mk_gamma_matrix(hurdle = FALSE)
  mixed <- mk_mixed_matrix()
  all_cont <- gamma_hurdle
  mixed_types <- list(
    b1 = "binary", b2 = "binary",
    c1 = "categorical", c2 = "categorical",
    x1 = "continuous", x2 = "continuous"
  )

  cases <- list(
    time_case(
      "problink_em_binary_patterns",
      function() problink_em(binary$formula, patterns = binary$patterns,
                             tol = 1e-6),
      reps, warmup
    ),
    time_case(
      "problink_em_gamma_hurdle",
      function() problink_em_gamma(comparison_matrix = gamma_hurdle,
                                   p0 = 0.05, p0M = 0.8,
                                   tol = 1e-5, maxits = 80),
      reps, warmup
    ),
    time_case(
      "problink_em_gamma_pure",
      function() problink_em_gamma(comparison_matrix = gamma_pure,
                                   p0 = 0.05, tol = 1e-5, maxits = 80),
      reps, warmup
    ),
    time_case(
      "problink_em_mixed",
      function() problink_em_mixed(comparison_matrix = mixed, types = mixed_types,
                                   p0 = 0.05, p0M = 0.8,
                                   tol = 1e-5, maxits = 80),
      reps, warmup
    ),
    time_case(
      "problink_em_mixed_all_continuous",
      function() problink_em_mixed(comparison_matrix = all_cont,
                                   types = rep("continuous", ncol(all_cont)),
                                   p0 = 0.05, p0M = 0.8,
                                   tol = 1e-5, maxits = 80),
      reps, warmup
    )
  )
  do.call(rbind, cases)
}

run_child <- function(label, repo, reps, warmup) {
  out <- tempfile(paste0("reclin2-bench-", label, "-"), fileext = ".csv")
  env <- c(
    "RECLIN2_BENCH_CHILD=1",
    paste0("RECLIN2_BENCH_REPO=", repo),
    paste0("RECLIN2_BENCH_OUT=", out),
    paste0("RECLIN2_BENCH_REPS=", reps),
    paste0("RECLIN2_BENCH_WARMUP=", warmup)
  )
  status <- system2(file.path(R.home("bin"), "Rscript"), script_path(), env = env)
  if (!identical(status, 0L)) stop("benchmark child failed for ", label)
  res <- read.csv(out, stringsAsFactors = FALSE)
  res$version <- label
  res
}

if (identical(Sys.getenv("RECLIN2_BENCH_CHILD"), "1")) {
  repo <- normalizePath(Sys.getenv("RECLIN2_BENCH_REPO"))
  reps <- as.integer(Sys.getenv("RECLIN2_BENCH_REPS", "3"))
  warmup <- as.integer(Sys.getenv("RECLIN2_BENCH_WARMUP", "1"))
  res <- run_benchmarks(repo, reps = reps, warmup = warmup)
  write.csv(res, Sys.getenv("RECLIN2_BENCH_OUT"), row.names = FALSE)
  quit(save = "no", status = 0)
}

repo <- repo_root()
reps <- as.integer(Sys.getenv("RECLIN2_BENCH_REPS", "3"))
warmup <- as.integer(Sys.getenv("RECLIN2_BENCH_WARMUP", "1"))
base_ref <- Sys.getenv("RECLIN2_BENCH_BASE_REF")
if (!nzchar(base_ref)) {
  dirty <- length(git_out(repo, c("status", "--porcelain"))) > 0L
  base_ref <- if (dirty) "HEAD" else "HEAD^"
}

base_wt <- tempfile("reclin2-phase1-base-", tmpdir = tempdir())
cat("Benchmarking current worktree against base ref ", base_ref, "\n", sep = "")
status <- system2("git", c("-C", repo, "worktree", "add", "--detach",
                           base_wt, base_ref))
if (!identical(status, 0L)) stop("failed to create base worktree")
on.exit({
  system2("git", c("-C", repo, "worktree", "remove", "--force", base_wt),
          stdout = FALSE, stderr = FALSE)
}, add = TRUE)

current <- run_child("current", repo, reps, warmup)
base <- run_child("base", base_wt, reps, warmup)

tab <- merge(base, current, by = "case", suffixes = c("_base", "_current"))
tab$speedup <- tab$median_elapsed_base / tab$median_elapsed_current
tab <- tab[order(tab$case), c("case", "median_elapsed_base",
                             "median_elapsed_current", "speedup",
                             "min_elapsed_base", "min_elapsed_current",
                             "max_elapsed_base", "max_elapsed_current",
                             "reps_base")]
names(tab)[names(tab) == "reps_base"] <- "reps"
num <- vapply(tab, is.numeric, logical(1))
tab[num] <- lapply(tab[num], function(x) round(x, 4))
print(tab, row.names = FALSE)
