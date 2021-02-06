
#' Field codes for BLAST results returned in tabular format
#' 
#' A dataset containing the field codes for indicating the desired format of tabular BLAST results,
#' as well as the description of each field.
#' 
#' @format A data frame with 50 rows and 2 variables
#' \describe{
#'   \item{field}{output format (`outfmt`) field code}
#'   \item{description}{brief description of information specified by the field}
#' }
#' @source `system("blastn -help")`
"fmtspec"
#> [1] "fmtspec"

#' List of datasets available for query through BLAST+/blastr
#' 
#' @format A character vector with 34 database names
#' @source `system("update_blastdb.pl --showall", intern=TRUE)`
"blast_dbs"
#> [1] "blast_dbs"
