#' Search for clonotypes in external databases
#'
#' @param .data Input dataset.
#' @param .db Either data frame, data table or any other object with columns, or a path to the database file if ".db.name" is provided.
#' @param .seq.col Name of the CDR3 sequence column in the database.
#' @param .v.col Provide FALSE to match only by ".seq.col" if ".db.name" is provided.
#' @param .db.name Either NA or the name of the external database. Provide one of the supported database names to
#' skip specifying column names. Supported databases:
#'
#' - "mcpas" - McPAS-TCR database, http://friedmanlab.weizmann.ac.il/McPAS-TCR/
#'
#' - "vdjdb" - VDJDB database, https://vdjdb.cdr3.net/search
#'
#' @export lookup_clonotypes
lookup_clonotypes <- function (.data, .db, .seq.col = NA, .v.col = NA, .db.name = NA) {
  # .db.name
  # - "mcpas"
  # - "vdjdb"

}
