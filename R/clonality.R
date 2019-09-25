#' Assess clonal proportions of repertoire
#'
#' @aliases clonality clonal.prop clonal_proportion top_proportion tail_proportion clonal_space_homeostasis
#' @concept diversity
#'
#' @description \code{repClonality} function encompasses several methods to measure
#' clonal proportions in a given repertoire.
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
#' @param .method A String with one of the following options: \code{"clonal.prop"},
#' \code{"homeo"}, \code{"top"} or \code{"tail"}.
#'
#' Set \code{"clonal.prop"} to compute clonal proportions or in other words
#' percentage of clonotypes required to occupy specified by \code{.perc} percent
#' of the total immune repertoire.
#'
#' Set \code{"homeo"} to analyse relative abundance (also known as clonal space homeostasis), which is defined as the
#' proportion of repertoire occupied by clonal groups.
#'
#' Set \code{"top"} to estimate relative abundance for the groups of top clonotypes in
#' repertoire. Use \code{".head"} to define index intervals.
#'
#' Set \code{"tail"} to estimate relative abundance for the groups of rare clonotypes
#' with low counts. Use \code{".bound"} to define the boundaries of clonotype groups.
#'
#' @param .perc A single numerical value ranging from 0 to 100.
#' @param .clone.types A named numerical vector with the boundaries of the half-closed
#' intervals that mark off clonal groups.
#' @param .head A numerical vector with ranges of the top clonotypes.
#' @param .bound A numerical vector with ranges of abundance for the rare clonotypes in
#' the dataset.
#'
#' @details Clonal proportion assessment is a different approach to estimate
#' repertoire diversity. When visualised, it allows for thorough examination of
#' immune repertoire structure and composition.
#'
#' In its core this type of analysis is similar to the relative species abundance
#' concept in ecology. Relative abundance is the percent composition of an organism
#' of a particular kind relative to the total number of organisms in the area.
#'
#' A stacked barplot of relative clonotype abundances can be therefore viewed as
#' a non-parametric approach to comparing their underlying distributions.
#'
#' @seealso \link{diversity}
#'
#' @examples
#' \dontrun{}
#'
#' @export repClonality rep.clo
repClonality <- function (.data, .method = c("clonal.prop", "homeo", "top", "tail"),
                          .perc = 10, .clone.types = c(Rare = .00001, Small = .0001,
                                                       Medium = .001, Large = .01, Hyperexpanded = 1),
                          .head = c(10, 100, 1000, 3000, 10000, 30000, 100000),
                          .bound = c(1, 3, 10, 30, 100)) {
  res = switch(.method[1],
               clonal.prop = clonal_proportion(.data, .perc),
               homeo = clonal_space_homeostasis(.data, .clone.types),
               top = top_proportion(.data, .head),
               tail = tail_proportion(.data, .bound),
               stop("Wrong method"))
  res
}

rep.clo <- repClonality

clonal_proportion <- function (.data, .perc = 10) {
  if (has_class(.data, 'list')) {
    res = t(sapply(.data, clonal_proportion, .perc = .perc))
  } else {
    .perc = min(.perc, 100)
    prop <- 0
    n <- 0
    col <- collect(select(.data, IMMCOL$count), n = Inf)[[1]]
    col.sum <- sum(col)
    while (prop < col.sum * (.perc / 100)) {
      n <- n + 1
      prop <- prop + col[n]
    }
    res = c(Clones = n, Percentage = 100 * signif((prop / col.sum), 3), Clonal.count.prop = n / length(col))
  }

  add_class(res, "immunr_clonal_prop")
}

clonal_space_homeostasis <- function (.data, .clone.types = c(Rare = .00001, Small = .0001,
                                                              Medium = .001, Large = .01, Hyperexpanded = 1)) {
  .clone.types <- c(None = 0, .clone.types)

  if (!has_class(.data, "list")) {
    .data <- list(Sample = .data)
  }

  mat <- matrix(0, length(.data), length(.clone.types) - 1, dimnames = list(names(.data), names(.clone.types)[-1]))
  .data <- lapply(.data, function (d) collect(select(d, IMMCOL$prop), n = Inf)[[1]])
  for (i in 2:length(.clone.types)) {
    mat[,i-1] <- sapply(.data, function (x) sum(x[x > .clone.types[i-1] & x <= .clone.types[i]]))
    colnames(mat)[i-1] <- paste0(names(.clone.types[i]), ' (', .clone.types[i-1], ' < X <= ', .clone.types[i], ')')
  }

  add_class(mat, "immunr_homeo")
}

top_proportion <- function (.data, .head = c(10, 100, 1000, 3000, 10000, 30000, 100000)) {
  if (has_class(.data, 'list')) {
    res = sapply(.head, function (.h) sapply(.data, top_proportion, .head = .h))
    colnames(res) = .head
    row.names(res) = names(.data)
  } else {
    .data = collect(select(.data, IMMCOL$prop), n = Inf)[[1]]
    res = sapply(.head, function (.h) sum(head(.data, .h)) / sum(.data) )
    names(res) = .head
  }
  add_class(res, "immunr_top_prop")
}

tail_proportion <- function (.data, .bound = c(1, 3, 10, 30, 100)) {
  if (has_class(.data, 'list')) {
    res = sapply(.bound, function (.b) sapply(.data, tail_proportion, .bound = .b))
    colnames(res) = .bound
    row.names(res) = names(.data)
  } else {
    .data = collect(select(.data, IMMCOL$count), n = Inf)[[1]]
    res = sapply(.bound, function (.b) { sum(.data[.data <= .bound] ) / sum(.data) } )
    names(res) = .bound
  }
  add_class(res, "immunr_tail_prop")
}
