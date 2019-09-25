#' Clonotype tracking across different samples.
#'
#' @aliases trackClonotypes tr.clo poisson_model norm_model
#'
#' @description
#' Track the dynamics of clonotypes across different repertoires.
#'
#' @usage
#' trackClonotypes(.data, .which = c(1, 15), .model = c("none", "poisson", "norm"),
#' .col = c("nuc", "aa"), .q = c(.025, .975), .sd = 1.05, .dup = c("merge", "remove"))
#'
#' tr.clo(...)
#'
#' poisson_model(.data, .q = c(.025, .975))
#'
#' norm_model(.data, .q = c(.025, .975), .sd = 1.05)
#'
#' @param .data List of repertoires or a single repertoire.
#' @param .which Either 1) a vector of length two `c(X, Y)`, where `Y` is the number of the most abundant clonotypes to take from `X`'th data frame in the input list; 2) or a character vector of sequences to take from all data frames; 3) or a data frame with two / three columns - first for sequences, and second/third for gene segments.
#' @param .model Which model to use for estimation of confidence intervals for clonotypes: "none" for no model, "poisson" for poisson model, "norm" for normal model.
#' @param .col Column with sequences, either "nuc" for nucleotide sequences or "aa" for amino acid sequences.
#' @param .q Quantiles for models.
#' @param .sd Standart deviation for normal model.
#' @param .dup What to do with duplicated sequences - either merge and sum up their counts ("merge") or remove sequences with lower counts ("remove").
#'
#' @return Data frame with input sequences and counts / estimations for each of the input repertoires.
#'
#' @examples
#' \dontrun{
#' immtr = trackClonotypes(immdata$data, "norm", c(3,100));
#' tidyr::drop_na(dcast(immtr, Sequence ~ Sample, value.var = c("Data")))
#' }
trackClonotypes <- function (.data, .which = c(1, 15), .model = c("none", "poisson", "norm"), .col = c("nuc", "aa"),
                             .q = c(.025, .975), .sd = 1.05, .dup = c("merge", "remove")) {
  .model = .model[1]
  .proc_duplicates <- function (.df, .dup, .col) {
    if (.dup == "remove") {
      .df
    } else {
      .df %>% group_by_(.col) %>% summarise(Count = sum(IMMCOL$count))
    }
  }

  if (!has_class(.data, "list")) {
    .data = list(Data = .data)
  }

  fun = switch(.model[1],
               none = function (x, ...) { cbind(Data = as.matrix(collect(select(x, Count))), Lo = 0, Hi = as.matrix(collect(select(x, Count)))) },
               poisson = function (x, .q, ...) poisson_model(x, .q),
               norm = function (x, .q, .sd, ...) norm_model(x, .q, .sd),
               stop("Unknown model's name"))

  .col = .col[1]
  .dup = .dup[1]
  seq_vec = c()
  data_index = .which[1]
  for (i in c(data_index, c(1:length(.data))[-data_index])) {
    # process duplicates
    seq_data = .proc_duplicates(.data[[i]], .dup[1], "IMMCOL$cdr3aa")

    # run the model
    model_counts = fun(seq_data, .q = .q, .sd = .sd)

    # get the sequences
    if (has_class(.which, "character")) {
      #
      # grep them real good
      #
      match
    } else {
      top_seq = .which[2]
      if (i == data_index) {
        seq_vec = head(collect(select(.data[[i]], IMMCOL$cdr3aa)), top_seq)[[1]]
        res = data.frame(Sequence = seq_vec, model_counts[1:top_seq,], Sample = names(.data)[i], stringsAsFactors = F)
      } else {
        seq_data = collect(select(.data[[i]], IMMCOL$cdr3aa))[[1]]
        inds = match(seq_vec, seq_data, nomatch = -1)
        model_counts = model_counts[inds[inds != -1],]
        if (is.null(dim(model_counts))) {
          model_counts = matrix(model_counts, nrow=1)
          colnames(model_counts) = c("Data", "Lo", "Hi")
        }
        res = rbind(res, data.frame(Sequence = seq_vec[inds != -1], model_counts, Sample = names(.data)[i], stringsAsFactors = F))
      }
    }
  }

  add_class(res, "immunr_dynamics")
}

tr.clo <- trackClonotypes

poisson_model <- function (.data, .q = c(.025, .975)) {
  col = collect(select(.data, Count))[[1]]
  lo = qpois(.q[1], col)
  hi = qpois(.q[2], col)
  res = cbind(Data = col, Lo = lo, Hi = hi)
  add_class(res, "immunr_dynamics")
}

norm_model <- function (.data, .q = c(.025, .975), .sd = 1.05) {
  col = collect(select(.data, Count))[[1]]
  lo = qnorm(.q[1], col, .sd)
  hi = qnorm(.q[2], col, .sd)
  res = cbind(Data = col, Lo = lo, Hi = hi)
  add_class(res, "immunr_dynamics")
}
