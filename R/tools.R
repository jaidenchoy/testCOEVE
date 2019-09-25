#' Check that distribution is sums up to one.
#'
#' @description
#' Check if the given .data is a distribution and normalise it if necessary with an optional laplace correction.
#'
#' @param .data Numeric vector of values.
#' @param .do.norm One of the three values - NA, T or F. If NA than check for distrubution (sum(.data) == 1)
#' and normalise if needed with the given laplace correction value. if T then do normalisation and laplace
#' correction. If F than don't do normalisaton and laplace correction.
#' @param .laplace Value for the laplace correction.
#' @param .na.val Replace all NAs with this value.
#' @param .warn.zero if T then the function checks if in the resulted vector (after normalisation)
#' are any zeros, and prints a warning message if there are some.
#' @param .warn.sum if T then the function checks if the sum of resulted vector (after normalisation)
#' is equal to one, and prints a warning message if not.
#'
#' @return Numeric vector.
#'
#' @export check_distribution
check_distribution <- function (.data, .do.norm = NA, .laplace = 1, .na.val = 0, .warn.zero = F, .warn.sum = T) {
  if (sum(is.na(.data)) == length(.data)) {
    warning("Error! Input vector is completely filled with NAs. Check your input data to avoid this. Returning a vector with zeros.\n")
    return(rep.int(0, length(.data)))
  }

  if (is.na(.do.norm)) {
    .data[is.na(.data)] <- .na.val
    if (sum(.data) != 1) {
      .data <- .data + .laplace
      .data <- prop.table(.data + .laplace)
    }
  } else if (.do.norm) {
    .data[is.na(.data)] <- .na.val
    .data <- prop.table(.data + .laplace)
    .warn.sum = F
  }

  if (.warn.zero && (0 %in% .data)) {
    warningText <- paste("Warning! There are", sum(which(.data == 0)), "zeros in the input vector. Function may produce incorrect results.\n")
    if(.laplace != 1){
      warningText <- paste(warningText, "To fix this try to set .laplace = 1 or any other small number in the function's parameters\n")
    }else{
      warningText <- paste(warningText, "To fix this try to set .laplace to any other small number in the function's parameters\n")
    }
    warning(warningText)
  }

  if (.warn.sum && sum(.data) != 1) {
    if (abs(sum(.data) - 1) < 1e-14) {
      # message("Note: difference between the sum of the input vector and 1 is ", (sum(.data) - 1), ", which may be caused by internal R subroutines and may not affect the result at all.\n")
      # Just skip this case - it means all is OK
    } else {
      warningText <- "Warning! Sum of the input vector is NOT equal to 1. Function may produce incorrect results.\n"
      if(!isTRUE(.do.norm)) warningText <- paste(warningText, "To fix this try to set .do.norm = TRUE in the function's parameters.\n")
      warning(warningText)
    }
  }

  .data
}


#' Add a new class attribute to the input object.
#'
#' @export
add_class <- function (.obj, .class) {
  class(.obj) <- c(class(.obj), .class)
  .obj
}


#' Check class
#'
#' @export
has_class <- function (.data, .class)
{
  .class %in% class(.data)
}


#' Set and update progress bars
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @aliases set_pb add_pb
#'
#' @export set_pb add_pb
set_pb <- function (.max)
{
  txtProgressBar(min = 0, max = .max, style = 3)
}

add_pb <- function (.pb, .value = 1)
{
  setTxtProgressBar(.pb, .pb$getVal() + .value)
}


# allele2segment
# segment2allele
# allele2family
# segment2family
#
# matrixSubgroups
#'
#' @export .quant_column_choice
.quant_column_choice <- function (x, .verbose = T) {
  x <- switch(x[1],
              count = IMMCOL$count,
              Count = IMMCOL$count,
              prop = IMMCOL$prop,
              Prop = IMMCOL$prop,
              proportion = IMMCOL$prop,
              Proportion = IMMCOL$prop,
              freq = IMMCOL$prop,
              { cat("You have specified an invalid column identifier. Choosed column: Clones\n"); IMMCOL$count })
  x
}


#' @export matrixdiagcopy
matrixdiagcopy <- function (.mat) {
  .mat[lower.tri(.mat)] <- t(.mat)[lower.tri(.mat)]
  .mat
}


#' Translate
#'
#' @aliases bunch_translate translate_bunch
#'
#' @export
bunch_translate <- function (.seq, .two.way = T) {
  .seq = toupper(.seq)
  .seq[grepl("N", .seq)] <- NA

  sapply(.seq, function (y) {
    if (!is.na(y)) {
      ny <- nchar(y)
      ny3 <- ny %/% 3
      tmp <- ''
      if (.two.way) {
        if (ny %% 3 != 0) { tmp <- paste0(rep('N', times = 3), collapse = '') }
        y <- paste0(substr(y, 1, 3*((ny3 %/% 2) +  (ny %% 2))),
                    tmp,
                    substr(y, 3*((ny3 %/% 2) +  (ny3 %% 2)) + (ny %% 3) + 1, ny),
                    collapse = '')
      } else {
        y <- substring(y, seq(1, nchar(y) - 2, 3), seq(3, nchar(y), 3))
      }
      paste0(AA_TABLE[unlist(strsplit(gsub("(...)", "\\1_", y), "_"))],collapse="")
    } else {
      NA
    }
  }, USE.NAMES = F)
}

translate_bunch <- bunch_translate

#' Check if everything is ok with group names
check_group_names <- function (.meta, .by) {
  names_to_check = c()
  if (is.null(.by)) {
    names_to_check = .by
  } else {
    names_to_check = names(.by)
  }

  for (i in 1:length(names_to_check)) {
    if (!(names_to_check[i] %in% colnames(.meta))) {
      message("Check failed: '", names_to_check[i], "' not in the metadata table!")
      return(F)
    }
  }

  return(T)
}


#' Get a character vector of samples' groups from the input metadata file
#'
group_from_metadata <- function (.by, .metadata, .sep="; ") {
  if (length(.by) == 1) {
    collect(select(.metadata, .by))[[1]]
  } else {
    do.call(paste, c(list(sep = .sep), lapply(1:length(.by), function (i) { collect(select(.metadata, .by[i]))[[1]] } )))
  }
}


# .meta == NA => .by is a vector of values to group by
# .meta != NA => .by is a name of the column in the metadata file
# .meta == NA & .by == NA => just choose the default column .data.sample.col for grouping
process_metadata_arguments <- function (.data, .by, .meta = NA, .data.sample.col = "Sample",
                                        .meta.sample.col = "Sample") {
  .data[[.data.sample.col]] = as.character(.data[[.data.sample.col]])
  if (!is.na(.by)[1]) {
    if (!is.na(.meta)[1]) {
      data_groups = group_from_metadata(.by, .meta)
      group_name = .by
      is_grouped = T
      data_group_names =.meta[[.meta.sample.col]]
    } else {
      if (length(.by) == length(.data[[.data.sample.col]])) {
        data_groups = .by
        data_group_names = .data[[.data.sample.col]]
        group_name = "Group"
        is_grouped = T
      } else {
        stop("Error: length of the input vector '.by' isn't the same as the length of the input data. Please provide vector of the same length.")
      }
    }
  } else {
    data_groups = unique(.data[[.data.sample.col]])
    group_name = .data.sample.col
    is_grouped = F
    data_group_names = unique(.data[[.data.sample.col]])

    if (length(data_groups) != length(data_group_names)) {
      stop("Error: number of samples doesn't equal to the number of samples in the metadata")
    }
  }

  names(data_groups) = data_group_names
  group_column = data_groups[.data[[.data.sample.col]]]
  list(groups=data_groups, group_column=group_column, group_names=data_group_names, name=group_name, is_grouped=is_grouped)
}


rename_column <- function (.data, .old, .new) {
  colnames(.data)[match(.old, colnames(.data))] = .new
  .data
}


#' Apply function to each pair of data frames from a list.
#'
#' @aliases apply_symm apply_asymm
#'
#' @description
#' Apply the given function to every pair in the given datalist. Function either
#' symmetrical (i.e. fun(x,y) == fun(y,x)) or assymmetrical (i.e. fun(x,y) != fun(y,x)).
#'
#' @usage
#' apply.symm(.datalist, .fun, ..., .diag = NA, .verbose = T)
#'
#' apply.asymm(.datalist, .fun, ..., .diag = NA, .verbose = T)
#'
#' @param .datalist List with some data.frames.
#' @param .fun Function to apply, which return basic class value.
#' @param ... Arguments passsed to .fun.
#' @param .diag Either NA for NA or something else != NULL for .fun(x,x).
#' @param .verbose if T then output a progress bar.
#'
#' @return Matrix with values M[i,j] = fun(datalist[i], datalist[j])
#'
#' @export apply_symm apply_asymm
apply_symm <- function (.datalist, .fun, ..., .diag = NA, .verbose = T) {
  res <- matrix(0, length(.datalist), length(.datalist))
  if (.verbose) pb <- set_pb(length(.datalist)^2 / 2 + length(.datalist)/2)
  for (i in 1:length(.datalist))
    for (j in i:length(.datalist)) {
      if (i == j && is.na(.diag)) { res[i,j] <- NA }
      else { res[i,j] <- .fun(.datalist[[i]], .datalist[[j]], ...) }
      if (.verbose) add_pb(pb)
    }
  if (.verbose) close(pb)
  row.names(res) <- names(.datalist)
  colnames(res) <- names(.datalist)
  matrixdiagcopy(res)
}

apply_asymm <- function (.datalist, .fun, ..., .diag = NA, .verbose = T) {
  res <- matrix(0, length(.datalist), length(.datalist))
  if (.verbose) pb <- set_pb(length(.datalist)^2)
  for (i in 1:length(.datalist))
    for (j in 1:length(.datalist)) {
      if (i == j && is.na(.diag)) { res[i,j] <- NA }
      else { res[i,j] <- .fun(.datalist[[i]], .datalist[[j]], ...) }
      if (.verbose) add_pb(pb)
    }
  if (.verbose) close(pb)
  row.names(res) <- names(.datalist)
  colnames(res) <- names(.datalist)
  res
}


#' Work in progress
#'
#' @export top
top <- function (.data, .n = 10) {
  if (has_class(.data, "list")) {
    lapply(.data, top, .n = .n)
  } else {
    top_n(.data, .n, Clones) %>% collect(n = Inf)
  }
}


#' Work in progress
#'
#' @aliases coding noncoding
#'
#' @export coding noncoding inframes outofframes
coding <- function (.data) {
  if (has_class(.data, "list")) {
    lapply(.data, coding)
  } else {
    # immdata$data[[1]] %>% mutate_(Len = "nchar(CDR3.nt)", Nonc = "CDR3.nt %like% '[*, ~]'") %>% filter((Len %% 3 == 0) & (!Nonc))

    d = collect(.data, n = Inf)
    d[grep('[*, ~]', d[[IMMCOL$cdr3aa]], invert = T), ]
  }
}

noncoding <- function (.data) {
  if (has_class(.data, "list")) {
    lapply(.data, noncoding)
  } else {
    d = collect(.data, n = Inf)
    d[grep('[*, ~]', d[[IMMCOL$cdr3aa]], invert = F), ]
  }
}

inframes <- function (.data) {
  if (has_class(.data, "list")) {
    lapply(.data, inframes)
  } else {
    d = collect(.data, n = Inf)
    subset(.data, nchar(.data[[IMMCOL$cdr3nt]]) %% 3 == 0)
  }
}

outofframes <- function (.data) {
  if (has_class(.data, "list")) {
    lapply(.data, outofframes)
  } else {
    d = collect(.data, n = Inf)
    subset(.data, nchar(.data[[IMMCOL$cdr3nt]]) %% 3 != 0)
  }
}


#' WIP
#' @export switch_type
switch_type <- function(type) {
  switch(tolower(type),
         nuc = IMMCOL$cdr3nt,
         nt = IMMCOL$cdr3nt,
         v = IMMCOL$v,
         j = IMMCOL$j,
         aa = IMMCOL$cdr3aa
  )
}

return_segments <- function (.gene) {
  stringr::str_replace_all(.gene, "\\*[[:digit:]]+", "")
}

return_families <- function (.gene) {
  stringr::str_replace_all(return_segments(.gene), "\\-[[:digit:]]+", "")
}
