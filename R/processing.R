error_correction <- function () {

}

decontamination <- function () {

}


#' Resample data frame using values from the column with number of clonesets.
#'
#' @aliases resample_col resample downsample_col downsample prop_sample_col prop_sample
#'
#' @description
#' Resample the data using the number of clones. Function \code{*_col} ???
#'
#' @usage
#' resample(.data, .n = -1)
#'
#' downsample(.data, .n)
#'
#' prop_sample(.data, .perc = 50)
#'
#' @param .data Data frame with the column \code{.col} or list of such data frames.
#' @param .n Number of values / reads / UMIs to choose.
#' @param .perc Percentage (0 - 100). See "Details" for more info.
#'
#' @return Subsampled data frame.
#'
#' @details
#' \code{resample}. Using multinomial distribution, compute the number of occurences for each cloneset, than remove zero-number clonotypes and
#' return resulting data frame. Probabilities for \code{rmultinom} for each cloneset is a percentage of this cloneset in
#' the \code{.col} column. It's a some sort of simulation of how clonotypes are chosen from the organisms. For now it's not working
#' very well, so use \code{downsample} instead.
#'
#' \code{downsample}. Choose \code{.n} clones (not clonotypes!) from the input repertoires without any probabilistic simulation, but
#' exactly computing each choosed clones. Its output is same as for \code{resample} (repertoires), but is more consistent and
#' biologically pleasant.
#'
#' \code{prop_sample}. Choose the first N clonotypes which occupies \code{.perc} percents of overall UMIs / reads.
#'
#' @seealso \link{rmultinom}, \link{clonal_proportion}
#'
#' @examples
#' \dontrun{
#' # Get 100K reads (not clones!).
#' immdata.1.100k <- resample(immdata[[1]], 100000, .col = "read.count")
#' }
#'
#' @export downsample
resample_col <- function (.data, .n = -1) {
  col_vec <- collect(select(.data, IMMCOL$prop), n = Inf)[[1]]
  if (.n[1] == -1) { .n <- sum(col_vec)}
  rmultinom(1, .n, col_vec)
}

resample <- function (.data, .n = -1) {
  if (has_class(.data, 'list')) {
    if (length(.n) != length(.data)) {
      .n <- c(.n, rep.int(-1, length(.data) - length(.n)))
    }
    return(lapply(.data, resample, .n = .n))
  }
  new_col <- resample_col(.data, .n)

  .data = collect(.data, n = Inf)[new_col != 0, ]
  new_col = new_col[new_col != 0]
  .data[[IMMCOL$count]] = new_col
  .data[[IMMCOL$prop]] = .data$Count / sum(.data[[IMMCOL$count]])
  .data[order(.data[[IMMCOL$prop]], decreasing = T),]
}

downsample_col <- function (.data, .n) {
  read_vec <- collect(select(.data, IMMCOL$count), n = Inf)[[1]]
  read_indices <- rep(0, sum(read_vec))
  Rcpp::cppFunction(
    '
    NumericVector fill_vec(NumericVector read_vec, NumericVector read_indices) {
      int dummy = 0;
      for (int i = 0; i < read_vec.size(); i++) {
      for (int j = dummy; j < (read_vec[i] + dummy); j++) {
      read_indices[j] = i;
      }
      dummy = dummy + read_vec[i];
      }
      return read_indices;
    }
    '
  )
  read_indices <- fill_vec(read_vec, read_indices)
  new_counts <- sample(read_indices, .n)
  new_reads <- rep(0, length(read_vec))
  Rcpp::cppFunction(
    '
    NumericVector fill_reads(NumericVector new_reads, NumericVector new_counts) {
      for (int i = 0; i < new_counts.size(); i++) {
      new_reads[new_counts[i]] = new_reads[new_counts[i]] + 1;
      }
      return new_reads;
    }
    '
  )
  fill_reads(new_reads, new_counts)
}

downsample <- function (.data, .n) {
  if (has_class(.data, 'list')) {
    return(lapply(.data, downsample, .n = .n))
  }

  new_col = downsample_col(.data, .n)

  .data = collect(.data, n = Inf)[new_col > 0, ]
  new_col = new_col[new_col > 0]
  .data[[IMMCOL$count]] = new_col
  .data[[IMMCOL$prop]] = new_col / sum(new_col)
  .data[order(.data[[IMMCOL$prop]], decreasing = T),]
}

prop_sample_col <- function (.data, .perc = 50) {
  res = clonal_proportion(.data, .perc)
  if (has_class(.data, "list")) {
    res[,1]
  } else {
    res[1]
  }
}

prop_sample <- function (.data, .perc = 50) {
  cols = prop_sample_col(.data, .perc)

  if (length(cols) == 1) {
    collect(.data, n = cols[1])
  } else {
    lapply(1:length(.data), function (i) {
      collect(.data[[i]], n = cols[i])
    })
  }
}


#' Get ranks / indices for counts of clonotypes
#'
#' @aliases set.rank set.index
#'
#' @description
#' Set new columns "Rank" and "Index":
#'
#' set.rank <==> .data$Rank = rank(.data[, .col], ties.method = 'average')
#'
#' set.index <==> .data$Index = 1:nrow(.data) in a sorted data frame by \code{.col}
#'
#' @param .data Data frame or list with data frames.
#' @param .col Character vector with name of the column to use for ranking or indexing.
#'
#' @return Data frame with new column "Rank" or "Index" or list with such data frames.
get_rank <- function (.data) {
  if (has.class(.data, 'list')) {
    lapply(.data, get_rank)
  } else {
    y <- 1 / collect(select(.data, IMMCOL$count), n = Inf)
    .data$Rank <- rank(y, ties.method = 'average')
    .data
  }
}

get_index <- function (.data) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, get_index))
  }
  .data = collect(select(.data, IMMCOL$count), n = Inf)
  .data <- .data[order(.data, decreasing = T), ]
  .data$Index <- 1:nrow(.data)
  .data
}
