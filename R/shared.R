#' @importFrom data.table setcolorder
#'
#' @export pubRep publicRepertoire
pubRep <- function (.data, .col = "aa+v", .quant = c("count", "prop"), .coding=T, .min.samples = 1, .db = NA, .monet.dir = NA, .use.dt = T, .verbose=T) {
  assertthat::assert_that(has_class(.data, "list"))

  .col = sapply(unlist(strsplit(.col, split="\\+")), switch_type, USE.NAMES=FALSE)

  .quant = .quant_column_choice(.quant[1])

  if (.use.dt) {
    data_list = lapply(.data, function (df) data.table(collect(select(df, c(.col, .quant)), n = Inf), key = .col))

    if (.verbose) { pb = set_pb(length(.data)) }
  } else {
    data_list = lapply(.data, select, c(.col, .quant))
  }

  res = data_list[[1]]
  # res = res %>% group_by_(.dots = .col) %>% summarise(Count = paste0()sum(.quant))
  res = res %>% group_by_(.dots = .col) %>% summarise_(Count = paste0("sum(", .quant, ")"))
  if (.use.dt) {
    colnames(res)[ncol(res)] = names(.data)[1]

    if (.verbose) { add_pb(pb) }
  }

  for (i in 2:length(.data)) {
    res = full_join(res, data_list[[i]] %>% group_by_(.dots = .col) %>% summarise_(Count = paste0("sum(", .quant, ")")),
                    by = .col, suffix = c(as.character(i-1), as.character(i)))
    if (.use.dt) {
      colnames(res)[ncol(res)] = names(.data)[i]

      if (.verbose) { add_pb(pb) }
    }
  }
  res = collect(res, n = Inf)
  if (.use.dt) { if (.verbose) { add_pb(pb) } }

  res = res %>% ungroup()

  # Add #samples column after .col
  res[["Samples"]] = rowSums(!is.na(as.matrix(res[, (length(.col) +1) : ncol(res), with = F])))
  res = res %>% dplyr::filter(Samples >= .min.samples)
  if (.verbose) { add_pb(pb) }

  col_samples = colnames(res)[(length(.col)+1):(ncol(res)-1)]

  if (.use.dt) {
    res = setcolorder(res, c(.col, "Samples", col_samples))
  } else {
    res = res[c(.col, "Samples", col_samples)]
  }

  res = res[order(res$Samples, decreasing = T),]

	if (!is.na(.monet.dir)) {
    con = DBI::dbConnect(MonetDBLite::MonetDBLite(), embedded = .monet.dir)
    ms = MonetDBLite::src_monetdblite(dbdir = .monet.dir)
    .db = list()
    .db$name = "public"
    .db$ms = ms
    .db$con = con
  }

  if (!is.na(.db)) {
    DBI::dbWriteTable(.db$con, .db$name, res, overwrite=TRUE)
    res = dplyr::tbl(.db$ms, .db$name)
  }

  add_class(res, "immunr_public_repertoire")
}


publicRepertoire <- pubRep


#' @export
public_matrix <- function (.data) {
  sample_i = match("Samples", colnames(.data)) + 1
  max_col = dim(.data)[2]
  as.matrix(collect(.data %>% dplyr::select(sample_i:max_col), n = Inf))
}


num2bin <- function(number, n_bits) {
  vec = as.numeric(intToBits(number))
  if (missing(n_bits)) {
    vec
  } else {
    vec[1:n_bits]
  }
}


get_public_repertoire_names <- function (.pr) {
  sample_i = match("Samples", names(.pr)) + 1
  names(.pr)[sample_i:ncol(.pr)]
}


#' @export pubRepFilter publicRepertoireFilter
pubRepFilter <- function (.pr, .meta, .by, .min.samples = 1) {
  assertthat::assert_that(has_class(.pr, "immunr_public_repertoire"))
  assertthat::assert_that(.min.samples > 0)

  if (!check_group_names(.meta, .by)) { return(NA) }

  data_groups = lapply(1:length(.by), function (i) { .meta[["Sample"]][.meta[[names(.by)[i]]] == .by[i]] })

  samples_of_interest = data_groups[[1]]

  if (length(.by) > 1) {
    for (i in 2:length(data_groups)) {
      samples_of_interest = intersect(data_groups[[i]], samples_of_interest)
    }
  }
  samples_of_interest = intersect(samples_of_interest, get_public_repertoire_names(.pr))

  if (length(samples_of_interest) == 0) {
    message("Warning: no samples found, check group values in the .by argument!")
    return(NA)
  }

  sample_i = match("Samples", colnames(.pr))
  indices = c(1:(match("Samples", colnames(.pr))), match(samples_of_interest, colnames(.pr)))
  new.pr = .pr %>% dplyr::select(indices)

  new.pr[["Samples"]] = rowSums(!is.na(as.matrix(new.pr[, (sample_i+1) : ncol(new.pr), with = F])))
  new.pr = new.pr %>% dplyr::filter(Samples >= .min.samples)

  add_class(new.pr, "immunr_public_repertoire")
}

publicRepertoireFilter <- pubRepFilter


#' @export pubRepApply publicRepertoireApply
pubRepApply <- function (.pr1, .pr2, .fun = function (x) log10(x[1]) / log10(x[2])) {
  col_before_samples = names(.pr1)[1:(match("Samples", colnames(.pr1))-1)]

  # tmp = apply(public_matrix(.pr1), 1, .inner.fun)
  tmp = rowMeans(public_matrix(.pr1), na.rm = T)
  .pr1[, (match("Samples", colnames(.pr1))+1):ncol(.pr1)] = NULL
  .pr1[["Quant"]] = tmp

  # tmp = apply(public_matrix(.pr2), 1, .inner.fun)
  tmp = rowMeans(public_matrix(.pr2), na.rm = T)
  .pr2[, (match("Samples", colnames(.pr2))+1):ncol(.pr2)] = NULL
  .pr2[["Quant"]] = tmp

  pr.res = dplyr::inner_join(.pr1, .pr2, by = col_before_samples)
  pr.res[["Samples.x"]] = pr.res[["Samples.x"]] + pr.res[["Samples.y"]]
  pr.res[, Samples.y := NULL]
  names(pr.res)[match("Samples.x", colnames(pr.res))] = "Samples"

  pr.res[["Result"]] = apply(pr.res[, c("Quant.x", "Quant.y")], 1, .fun)

  add_class(pr.res, "immunr_public_repertoire_apply")
}

publicRepertoireApply <- pubRepApply


# pubRepAnalysis <- function (.pr, .method, .by = NA, .meta = NA, .min.samples = 1) {
#
# }

