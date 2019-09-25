# make Monet DB from the list of files

#' Various functions for testing purposes
#'
load_test_data <- function () {
  # data(ms_small, package = "immunarch")
  data(immdata)
  if (require("MonetDBLite")) {
    immdata = test_make_db(immdata$data, immdata$meta)
  }
  # rm(ms_data)
  immdata
}

test_make_db <- function (.data, .meta = NA) {
  assertthat::assert_that(has_class(.data, "list"))

  dbdir = tempdir()
  con = DBI::dbConnect(MonetDBLite::MonetDBLite(), embedded = dbdir)

  for (i in 1:length(.data)) {
    DBI::dbWriteTable(con, names(.data)[i], .data[[i]], overwrite=TRUE)
  }

  ms = MonetDBLite::src_monetdblite(dbdir = dbdir)
  res_db = list()
  for (i in 1:length(.data)) {
    res_db[[names(.data)[i]]] = dplyr::tbl(ms, names(.data)[i])
  }

  if (is.na(.meta)) {
    res_db
  } else {
    list(data = res_db, meta = .meta)
  }
}

test_make_dt <- function (.data, .meta = NA) {
  require(dtplyr)
}

test_make_spark <- function (.data, .meta = NA) {

}

test_make_df <- function (.data, .meta = NA) {
  require(sparklyr)
}


# ms_data = list()
# for (file in list.files("~/Downloads/vdjtools-examples-master/ms_immunr/data/")) {
#   name = stringr::str_sub(file, 1, stringr::str_locate(file, ".txt")[1]-1)
#   print(name)
#   ms_data[[name]] = head(parse_vdjtools(stringr::str_c("~/Downloads/vdjtools-examples-master/ms_immunr/data/", file)), round(runif(1, 70000, 150000)))
# }
# ms_meta = readr::read_tsv("~/Downloads/vdjtools-examples-master/ms_immunr/metadata.txt")
# ms_data = list(data = ms_data, meta = ms_meta)
# save(ms_data, file="data/ms.rda", compress = "bzip2", compression_level = 9)


# tmp = list()
# tmp[[1]] = parse_vdjtools("~/Downloads/same_person/F87-1.txt")
# tmp[[2]] = parse_vdjtools("~/Downloads/same_person/M25-1.txt")
# tmp[[3]] = parse_vdjtools("~/Downloads/same_person/NA93-1.txt")
# sapply(tmp, nrow)
#
# split_data <- function (.df, .num) {
#   .df = .df[sample(1:nrow(.df), size = nrow(.df), ), ]
#   .df = .df[order(.df$Count, decreasing = T), ]
#   res = list()
#   for (n in 1:.num) {
#     res[[n]] = .df[sample(1:nrow(.df), size = nrow(.df) %/% .num, replace = F), ]
#   }
#
#   print(sapply(res, nrow))
#   res
# }
#
# TEST_DATA = list()
# tos = c(5,7,3)
# k = 1
# for (i in 1:3) {
#   tmp2 = split_data(tmp[[1]], tos[i])
#   for (j in 1:length(tmp2)) {
#     TEST_DATA[[k]] = tmp2[[j]][order(tmp2[[j]]$Count, decreasing = T),]
#     k = k + 1
#   }
# }
# names(TEST_DATA) = paste0("S", 1:15)
