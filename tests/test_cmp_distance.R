
library(reclin2)
source("helpers.R")

## --- cmp_distance: two-argument form ---
cmp <- cmp_distance()
d <- cmp(c(1.0, 2.5, 3.0, NA), c(1.0, 3.0, 1.0, 2.0))
expect_equal(d[1], 0)
expect_equal(d[2], 0.5)
expect_equal(d[3], 2.0)
stopifnot(is.na(d[4]))

## --- cmp_distance: one-argument (binary) form ---
b <- cmp(d)
expect_equal(b[1], TRUE)
expect_equal(b[2], FALSE)
expect_equal(b[3], FALSE)
expect_equal(b[4], FALSE)

## --- cmp_distance: scale parameter ---
cmp2 <- cmp_distance(scale = 10)
d2 <- cmp2(c(0, 5, 10), c(0, 0, 0))
expect_equal(d2, c(0, 0.5, 1.0))

## --- cmp_distance: transform parameter ---
cmp3 <- cmp_distance(transform = log1p)
d3 <- cmp3(c(0, 1, 3), c(0, 0, 0))
expect_equal(d3, log1p(c(0, 1, 3)))

## --- cmp_absdist is shorthand ---
cmp4 <- cmp_absdist()
d4 <- cmp4(c(3, 7), c(1, 7))
expect_equal(d4, c(2, 0))

## --- cmp_distance: scale must be a positive finite scalar (FIX 4) ---
expect_error(cmp_distance(scale = 0))
expect_error(cmp_distance(scale = -1))
expect_error(cmp_distance(scale = c(1, 2)))
expect_error(cmp_distance(scale = NA_real_))
expect_error(cmp_distance(scale = Inf))

## --- cmp_distance: a transform producing negatives warns (FIX 4) ---
## log on (scaled) distances < 1 yields negative values; the hurdle-gamma model
## needs non-negative distances, so this must warn.
cmp_log <- cmp_distance(transform = log)
expect_warning(cmp_log(c(0.5, 2), c(0, 0)))   # |0.5-0| = 0.5 -> log < 0
cat("  cmp_distance scale/transform validation (FIX 4): PASS\n")

cat("PASS: test_cmp_distance.R\n")
