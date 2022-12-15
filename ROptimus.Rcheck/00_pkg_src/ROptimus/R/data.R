################################################################################
#' Tutorial 5 Genomic Contact Data
#'
#' A dataset containing 734 genomic contact pairs for reproducing Tutorial 5.
#' Here, we take a set of 734 pairs of i and j 40-kilobase DNA regions that are
#' known to be in contact inside a cell. Binning the DNA (e.g., a single chromosome)
#' into 40-kb regions, each region is represented as a single integer that is equal
#' to its end position divided by the length of the region, which is 40 kb. For
#' instance, the 1st region, with start and end positions at the 1st and 40000th
#' nucleotides, respectively, is denoted as 1 (40000th base / 40000 bases = 1).
#' This simplifies the notation for a contact between two regions to a pair of
#' positive integers as in used in this dataset.
#'
#' @format A data frame with 734 rows and 2 variables:
#' \describe{
#'   \item{i}{Integer representing 40-kb contact region partner of j}
#'   \item{j}{Integer representing 40-kb contact region partner of i}
#' }
"IJ_ORIG"

#' Tutorial 5 m() function
#'
#' The m() function supplied to Optimus() for reproducing Tutorial 5. This is
#' for the small executable example.
#'
#' @format Function
"ex.m.fun"

#' Tutorial 5 u() function
#'
#' The u() function supplied to Optimus() for reproducing Tutorial 5. This is
#' for the small executable example.
#'
#' @format Function
"ex.u.fun"

#' Tutorial 5 r() function
#'
#' The r() function supplied to Optimus() for reproducing Tutorial 5. This is
#' for the small executable example.
#'
#' @format Function
"ex.r.fun"
