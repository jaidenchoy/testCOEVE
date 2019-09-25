#' Calculate the overlap between repertoires
#'
#' @importFrom data.table data.table
#'
#' @description The \code{repOverlap} function is designed to analyse the overlap between
#' two or more repertoires. It contains a number of methods to compare immune receptor
#' sequences that are shared between individuals.
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
#' @param .method A string that specifies the method of analysis or a combination of
#' methods. The \code{repOverlap} function supports following basic methods:
#' "public", "overlap", "jaccard", "tversky", "cosine", "morisita".
#'
#' @param .col A string that specifies the column to be processed. Pass "nuc" for
#' nucleotide sequence or "aa" for amino acid sequence.
#' @param .quant Select the column with data to evaluate
#' @param .a,.b Alpha and beta parameters for Tversky Index. Default values give
#' the Jaccard index measure.
#' @param .verbose if T then output the progress.
#' @param .dup Defines the duplicates' behaviour. Pass "merge" or "remove".
#' @param .step An integer that defines the step of incremetal overlap (Note! Currently,
#' top+overlap, or "top+shared" and "top+morisita".
#'
#' @details "public" and "shared" are synonyms that exist
#' for the convenience of researchers.
#'
#' The "overlap" coefficient is a similarity measure that measures the overlap between two finite sets.
#'
#' The "jaccard" index is conceptually a percentage of how many objects two sets
#' have in common out of how many objects they have total.
#'
#' The "tversky" index is an asymmetric similarity measure on sets that
#' compares a variant to a prototype.
#'
#' The "cosine" index is a measure of similarity between two non-zero vectors
#' of an inner product space that measures the cosine of the angle between them.
#'
#' The "morisita" index measures how many times it is more likely to randomly
#' select two sampled points from the same quadrat (the dataset is covered by a
#' regular grid of changing size) than it would be in the case of a random
#' distribution generated from a Poisson process.
#'
#' @seealso \link{jaccard_index}
#'
#' @export repOverlap rep.ov
repOverlap <- function (.data,
                        .method = c("public", "overlap", "jaccard", "tversky", "cosine", "morisita", "top+shared", "top+morisita"),
                        .col = "nuc",
                        .quant = c("count", "prop"),
                        .a = .5,
                        .b = .5,
                        .verbose = T,
                        .dup = c("merge", "remove")
                        ) {
  assertthat::assert_that(has_class(.data, "list"))

  .col = unlist(strsplit(.col, split="\\+"))
  .col = sapply(.col, switch_type, USE.NAMES=FALSE)
  .dup = .dup[1]
  .method = .method[1]
  .quant = .quant_column_choice(.quant[1])

  if (.method == "cosine" || .method == "morisita") {
    .col = c(.col, IMMCOL$count)
  }

  for (i in 1:length(.data)) {
    .data[[i]] = collect(select_(.data[[i]], .dots = .col), n = Inf)
  }

  res = switch(.method,
               shared = num_shared_clonotypes(.data),
               public = num_shared_clonotypes(.data),
               overlap = apply_symm(.data, overlap_coef, .diag = NA, .verbose = .verbose),
               jaccard = apply_symm(.data, jaccard_index, .diag = NA, .verbose = .verbose),
               tversky = apply_symm(.data, tversky_index, .a = .a, .b = .b, .diag = NA, .verbose = .verbose),
               cosine = apply_symm(.data, cosine_sim, .quant=.quant, .diag = NA, .verbose = .verbose),
               morisita = apply_symm(.data, morisita_index, .quant=.quant, .dup = .dup, .diag = NA, .verbose = .verbose),
               stop("You entered the wrong method! Please, try again."))
  if (length(res) == 4) {
    res = res[1,2]
  } else {
    res = res
  }

  add_class(res, "immunr_ov_matrix")
}

rep.ov <- repOverlap


#' Number of shared clonotypes.
#'
#' @param .data List of repertoires of length 2 or more - either data frames, data tables or MonetDB tables.
#' @param .columns Character vector, which columns to use.
#' @param .norm Logical, if T than normalise the overlap value by multiplication of repertoires' sizes.
#' @param .head Overlap only .head clonotypes from each repertoire. -1 means overlap all clonotypes.
#' @param .verbose if T then output the progress.
#'
#' @return If length of .data is 2 than a single values, otherwise matrix with overlap values.
#'
#' @export
num_shared_clonotypes <- function (.data) {
  res = matrix(0, length(.data), length(.data))

  for (i in 1:(length(.data)-1)) {
    data_i = collect(.data[[i]], n = Inf)
    # res[i,i] = dplyr::count(data_i)[[1]]
    for (j in (i+1):length(.data)) {
      data_j = collect(.data[[j]], n = Inf)
      # print(dplyr::count(dplyr::intersect(data_i, data_j)))
      res[i,j] = dplyr::count(dplyr::intersect(data_i, data_j))[[1]]
      res[j,i] = res[i,j]
    }
  }

  diag(res) = NA

  row.names(res) = names(.data)
  colnames(res) = names(.data)

  if (length(.data) == 2) { res[1,2] }
  else { res }
}


#' Various measures for overlapping sets of objects.
#'
#' @aliases overlap_coef jaccard_index tversky_index morisita_index cosine_sim horn_index
#'
#' @description
#' Functions for computing similarity between two vectors or sets.
#'
#' - Overlap cofficient is a similarity measure related to the Jaccard index that measures the overlap between two sets, and is defined as the size of the intersection divided by the smaller of the size of the two sets.
#'
#' - Jaccard index is a statistic used for comparing the similarity and diversity of sample sets.
#'
#' - Tversky index is an asymmetric similarity measure on sets that compares a variant to a prototype and is related to Jaccard index.
#'
#' - Cosine similarity is a measure of similarity between two vectors of an inner product space that measures the cosine of the angle between them.
#'
#' - Morisita's overlap index is a statistical measure of dispersion of individuals in a population. It is used to compare overlap among samples (Morisita 1959). This formula is based on the assumption that increasing the size of the samples will increase the diversity because it will include different habitats (i.e. different faunas).
#'
#' - Horn's overlap index is based on Shannon's entropy.
#'
#' Use the \link{repOverlap} function for computing similarities of clonesets.
#'
#' @usage
#' overlap_coef(.x, .y)
#'
#' jaccard_index(.x, .y)
#'
#' tversky_index(.x, .y, .a = .5, .b = .5)
#'
#' cosine_sim(.x, .y, .do.norm = c(NA, T, F), .laplace = 0)
#'
#' morisita_index(.x, .y, .dup = c("merge", "remove"))
#'
#' horn_index(.x, .y, .dup = c("merge", "remove"))
#'
#' @param .x,.y Character vectors or data frames, sets of objects to intersect, see "Details"
#' for more information.
#' @param .a,.b Alpha and beta parameters for Tversky Index. Default values gives the Jaccard index measure.
#' @param .dup Character value, either "remove" to completely remove duplicated values from the first columns
#' of .x and .y, or "merge" (by default) to sum up values of duplicated elements.
#'
#' @details
#' For \code{cosine_sim} .x and .y are two numeric vectors.
#' For \code{overlap_coef}, \code{jaccard_index} and \code{tversky_index} .x and .y are character vectors.
#' For \code{morisita_index} and \code{horn_index} .x and .y are two data frames or data tables with character
#' vectors as first columns and numeric (i.e., number of such elements) as second.
#'
#' @seealso \link{repOverlap}
#'
#' @export overlap_coef jaccard_index tversky_index cosine_sim morisita_index horn_index
overlap_coef <- function (.x, .y) {
  UseMethod('overlap_coef')
}

overlap_coef.default <- function(.x, .y) {
  .x = collect(.x, n = Inf)
  .y = collect(.y, n = Inf)
  collect(count(dplyr::intersect(.x, .y)))[[1]]/min(nrow(.x), nrow(.y))
}

overlap_coef.character <- function(.x, .y) {
  length(dplyr::intersect(.x, .y)) / min(length(.x), length(.y))
}


jaccard_index <- function (.x, .y) {
  UseMethod('jaccard_index')
}

jaccard_index.default <- function(.x, .y) {
  .x = collect(.x, n = Inf)
  .y = collect(.y, n = Inf)
  intersection = collect(count(dplyr::intersect(.x, .y)))[[1]]
  intersection/(collect(count(.x))[[1]] + collect(count(.y))[[1]] + intersection)
}

jaccard_index.character <- function(.x, .y) {
  intersection = length(dplyr::intersect(.x, .y))
  intersection/(length(.x) + length(.y) + intersection)
}

tversky_index <- function (.x, .y, .a = .5, .b = .5) {
  UseMethod('tversky_index')
}

tversky_index.default <- function(.x, .y, .a = .5, .b = .5) {
  .x = collect(.x, n = Inf)
  .y = collect(.y, n = Inf)
  intersection = collect(count(dplyr::intersect(.x, .y)))[[1]]
  intersection/(.a * collect(count(dplyr::setdiff(.x, .y)))[[1]] + .b * collect(count(dplyr::setdiff(.y, .x)))[[1]] + intersection)
}

tversky_index.character <- function(.x, .y, .a = .5, .b = .5) {
  intersection = length(dplyr::intersect(.x, .y))
  intersection/(.a * length(dplyr::setdiff(.x, .y)) + .b * length(dplyr::setdiff(.y, .x)) + intersection)
}

cosine_sim <- function (.x, .y, .quant) {
  UseMethod('cosine_sim')
}

cosine_sim.default <- function(.x, .y, .quant) {
  .x = collect(.x, n = Inf)
  .y = collect(.y, n = Inf)
  col_name = colnames(.x)[ colnames(.x) != .quant]
  joined_set = full_join(.x, .y, by = col_name, suffix = c("_a", "_b"))
  alpha_col = paste0(.quant, "_a")
  beta_col = paste0(.quant, "_b")
  new_set = collect(select_(joined_set, .dots=c(alpha_col, beta_col)))
  first_col = new_set[,1]
  second_col = new_set[,2]
  first_col[is.na(first_col)] = 0
  second_col[is.na(second_col)] = 0
  sum(first_col*second_col)/(sqrt(sum(first_col*first_col)) * sqrt(sum(second_col*second_col)))
}

cosine_sim.numeric <- function(.x, .y, .quant) {
  df = rbind(.x, .y)
  sum(.x * .y) / (sqrt(rowSums(df^2))[1] * sqrt(rowSums(df^2))[2])[[1]]
}

morisita_index <- function (.x, .y, .quant=IMMCOL$count, .dup = c("merge", "remove")) {
  .x = collect(.x, n = Inf)
  .y = collect(.y, n = Inf)
  col_name = colnames(.x)[colnames(.x) != .quant]
  names(.x)[colnames(.x) == .quant] = "Count"
  names(.y)[colnames(.y) == .quant] = "Count"

  if (!has_class(.x, "tbl_sql")) {
    .x = tbl_dt(.x)
    .y = tbl_dt(.y)
  }

  if (.dup[1] == "remove")  {
    sorted_alpha = arrange(.x, desc(Count))
    sorted_beta = arrange(.y, desc(Count))
    alpha_prepared = distinct_(sorted_alpha, .dots = col_name, .keep_all = TRUE)
    beta_prepared = distinct_(sorted_beta, .dots = col_name, .keep_all = TRUE)
  }
  else {
    alpha_sorted = group_by_(.x, .dots = col_name)
    beta_sorted = group_by_(.y, .dots = col_name)
    # expr <- lazyeval::interp(~sum(.quant), .quant = as.name(.quant))
    alpha_prepared <- summarise(alpha_sorted, Count_a = sum(Count))
    beta_prepared <- summarise(beta_sorted, Count_b = sum(Count))
  }

  if (has_class(alpha_prepared, "data.table")) {
    alpha_len = nrow(alpha_prepared)
    beta_len = nrow(beta_prepared)
  } else {
    alpha_len = collect(count(alpha_prepared))[[1]]
    beta_len = collect(count(beta_prepared))[[1]]
  }

  joined_set = select_(full_join(alpha_prepared, beta_prepared, by = col_name), .dots = c(col_name, "Count_a", "Count_b"))
  # joined_set = select_(full_join(alpha_prepared, beta_prepared, by = col_name,  suffix = c("", "")), .dots = c(col_name, "Count_a", "Count_b"))
  alpha_col = "Count_a"
  beta_col= "Count_b"
  joined_set = ungroup(joined_set)
  new_set = collect(select_(joined_set, .dots=c(alpha_col, beta_col)))
  first_col = new_set[,1]
  second_col = new_set[,2]
  first_col[is.na(first_col)] = 0 # not is na - this would be a length
  second_col[is.na(second_col)] = 0
  2 * sum(first_col*second_col) / ((sum(first_col * first_col)/alpha_len^2 + sum(second_col * second_col)/beta_len^2) * (as.numeric(alpha_len) * as.numeric(beta_len)))
}

horn_index <- function (.x, .y, .dup = c("merge", "remove")) {
  .x[,1] <- as.character(.x[,1])
  .y[,1] <- as.character(.y[,1])
  colnames(.x) <- c('Species', 'Count')
  colnames(.y) <- c('Species', 'Count')
  if (.do.unique) {
    .x <- .x[!duplicated(.x[,1]),]
    .y <- .y[!duplicated(.y[,1]),]
  }
  .x[,2] <- as.numeric(.x[,2]) / sum(.x[,2])
  .y[,2] <- as.numeric(.y[,2]) / sum(.y[,2])
  merged <- merge(.x, .y, by = 'Species', all = T)
  merged[is.na(merged)] <- 0
  rel.12 <- merged[,2] / merged[,3]
  rel.12[merged[,3] == 0] <- 0
  rel.21 <- merged[,3] / merged[,2]
  rel.21[merged[,2] == 0] <- 0
  1 / log(2) * sum(merged[,2] / 2 * log(1 + rel.21) + merged[,3] / 2 * log(1 + rel.12))
}


#' Sequential overlap of top-by-count clonotypes.
#'
#' @param .data Either data frame, data table, MonetDB table or list of those, list of repertoires.
#'
#' @export
top_overlap <- function (.data, .fun, .step = 1000, ...) {
  step_vec = seq(.step, min(sapply(.data, function (d) nrow(collect(select(d, 1))))), .step)

  res = lapply(step_vec, function (i) {
    .fun(.data, ...)
  })
  names(res) = step_vec
  res
}
