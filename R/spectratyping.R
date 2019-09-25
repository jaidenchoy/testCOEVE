#' Analyse spectratype of immune repertoire
#'
#' @description
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
#' @param .quant Select the column with data to evaluate
#' @param .col A string that specifies the column to be processed. Pass "nuc" for
#' nucleotide sequence or "aa" for amino acid sequence.
#' @param .gene A string with either "v", "j" or NA.
#'
#' @export spectratype vis.immunr_spectr vis.immunr_spectr_nogene
spectratype <- function (.data, .quant = c("id", "count"), .col = c("nuc", "aa"), .gene = "v") {
  assertthat::assert_that(!has_class(.data, "list"))

  which_cols = switch(.col[1],
                      nuc = IMMCOL$cdr3nt,
                      aa = IMMCOL$cdr3aa,
                      stop("Wrong column name, choose either 'nuc' or 'aa'"))

  which_cols = c(which_cols, IMMCOL$count)

  if (!is.na(.gene)) {
    .gene = paste0(toupper(.gene), ".name")
    which_cols = c(which_cols, .gene)
  }

  .data = collect(select_(.data, .dots = which_cols), n = Inf)
  .data$Length = nchar(.data[[which_cols[1]]])
  if (.quant[1] == "id") {
    .data[[IMMCOL$count]] = 1
  }

  if (is.na(.gene)) {
    df = summarise(group_by(do.call(select_, list(.data, "Length", IMMCOL$count)), Length), Val = sum(Count))
    add_class(as.data.frame(df), "immunr_spectr_nogene")
  } else {
    .gene = switch(tolower(substr(.gene, 1, 1)), v = "V.name", j = "J.name")
    df = do.call(select_, list(.data, "Length", IMMCOL$count, .gene))
    df = group_by_(df, "Length", .gene)
    df = summarise_(df, Val = paste0("sum(", IMMCOL$count, ")"))
    df$Gene = df[[.gene]]
    df[[.gene]] = NULL
    df = df[order(df$Val, decreasing = T),]

    add_class(as.data.frame(df), "immunr_spectr")
  }
}

vis.immunr_spectr <- function (.data, .main = "Spectratype", .legend = "Gene segment", .labs = c("CDR3 length", NA)) {
  dup = which(cumsum(!duplicated(.data$Gene)) == 12)[1]
  if (length(dup)) {
    uniq = unique(.data$Gene[1:(dup - 1)])
    .data$Gene[!(.data$Gene %in% uniq)] = "Z"  # <- dirty hack to avoid factors
  }

  p = ggplot() + geom_bar(aes(x = Length, y = Val, fill = Gene), data = .data, stat = "identity") +
    scale_fill_manual(name = .legend, breaks = c(sort(uniq), "Z"),
                      labels=c(sort(uniq), "Other"),
                      values = c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "grey75")) +
    ggtitle(.main)

    if (!is.na(.labs[2])) {
      p = p + ylab(.labs[2])
    } else {
      p = p + ylab("Count")
    }
    p = p + xlab(.labs[1])
    p + theme_linedraw()
}

vis.immunr_spectr_nogene <- function (.data, .main = "Spectratype", .legend = "Gene segment", .labs = c("CDR3 length", NA)) {
  p = ggplot() + geom_bar(aes(x = Length, y = Val), data = .data, stat = "identity")

  if (!is.na(.labs[2])) {
    p = p + ylab(.labs[2])
  } else {
    p = p + ylab("Count")
  }
  p = p + xlab(.labs[1])
  p + theme_linedraw()
}
