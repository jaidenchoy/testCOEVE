#' Explore basics of repertoire
#'
#' @aliases repExplore
#'
#' @description The \code{repExplore} function calculates the basic statistics of
#' repertoire: the number of unique immune receptor clonotypes, their relative abundances,
#' and sequence length distribution across the input dataset.
#'
#' @param .data The data to be processed. Can be \link{data.frame},
#' \link{data.table}, or a list of these objects.
#'
#' Every object must have columns in the immunarch compatible format.
#' \link{immunarch_data_format}
#'
#' Competent users may provide advanced data representations:
#' DBI database connections, Apache Spark DataFrame from \link{copy_to} or a list
#' of these objects. They are supported with the same limitations as basic objects.
#'
#' Note: each connection must represent a separate repertoire.
#'
#' @param .method A string that specifies the method of analysis. It can be
#' either "volume", "count" or "len".
#'
#' When .method is set to "volume" the repExplore calculates the number of unique
#' clonotypes in the input data.
#'
#' When .method is set to "count" the repExplore calculates relative abundance of
#' different receptors in the input data.
#'
#' When .method is set to "len" the repExplore calculates the distribution of
#' CDR3 sequence lengths.
#'
#' @param .col A string that specifies the column to be processed. Pass "nt" for
#' nucleotide sequence or "aa" for amino acid sequence.
#'
#' @param .coding If \code{TRUE}, then only coding sequences will be analysed.
#'
#' @seealso \link{vis.immunr_exp_vol}
#'
#' @examples
#' \dontrun{
#' # Calculate statistics and generate a visual output with vis()
#' vis(repExplore(your_data$data, .method = "volume"))
#'
#' vis(repExplore(your_data$data, .method = "count"))
#'
#' vis(repExplore(your_data$data, .method = "len"))
#' }
#'
#' @export repExplore
repExplore <- function (.data, .method = c("volume", "count", "len"), .col = c("nt", "aa"), .coding=T) {
  if (!has_class(.data, "list")) {
    .data = list(.data)
  }

  if (.coding) {
    .data = coding(.data)
  }

  if (.method[1] == "volume") {
    res = add_class(data.frame(Sample = names(.data),
                               Volume = sapply(.data, function (df) {
                                           nrow(collect(select(df, IMMCOL$count), n = Inf))
                                         }),
                               stringsAsFactors = F), "immunr_exp_vol")
  } else if (.method[1] == "count") {
    res = lapply(1:length(.data), function (i) {
      count_table = table(collect(select(.data[[i]], IMMCOL$count), n = Inf))
      data.frame(Sample = names(.data)[i],
                 Clone.num = as.numeric(names(count_table)),
                 Clonotypes = as.vector(count_table),
                 stringsAsFactors = F)
      })

    res = do.call(rbind, res)
    res = add_class(res, "immunr_exp_count")
  } else if (.method[1] == "len") {
    seq_col = switch(.col[1], nt = IMMCOL$cdr3nt, aa = IMMCOL$cdr3aa, stop("Unknown sequence column"))
    res = lapply(.data, function (df) { table(nchar(collect(select(df, seq_col), n = Inf)[[1]])) })
    res = lapply(1:length(.data), function (i) {
            data.frame(Sample = names(.data)[i],
                       Length = as.numeric(names(res[[i]])),
                       Count = as.vector(res[[i]]))
      } )
    res = do.call(rbind, res)
    res = add_class(res, "immunr_exp_len")
  } else {
    stop("Unknown method")
    return(NULL)
  }

  res
}

rep.ex <- repExplore
