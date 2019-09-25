##### Main IO functions #####


#' Load repertoire files to the workspace
#'
#' @importFrom readr read_delim read_tsv read_csv col_integer col_character col_double col_logical col_guess cols write_lines
#' @importFrom stringr str_split
#'
#' @description The \code{repLoad} function loads repertoire files
#' into R workspace in the immunarch format where you can immediately use them for
#' the analysis. \code{repLoad} automatically detects the right format for
#' your files, so all you need is simply provide the path to your files.
#'
#' See "Details" for more information on supported formats. See "Examples" for
#' diving right into it.
#'
#' @param .path A character string specifying the path to the input data.
#' Input data can be one of the following:
#'
#' - a single repertoire file.
#' In this case \code{repLoad} returns an R \link{data.frame};
#'
#' - a vector of paths to repertoire files.
#' In this case \code{repLoad} returns a list of R \link{data.frame}s;
#'
#' - a path to the folder with repertoire files and, if available, metadata file.
#' If the metadata file is not presented, then the \code{repLoad} returns a list of R \link{data.frame}s.
#' If the metadata file if presented, then the \code{repLoad} returns a list with two elements "data" and "meta".
#' "data" is an another list with repertoire R \link{data.frame}s. "meta" is a data frame with the metadata.
#'
#' @param .format A character string specifying what format to use. See "Details" for more information on supported formats.
#'
#' Leave NA (which is default) if you want `immunarch` to detect formats automatically.
#'
#' @details
#' You can load the data either from a single file or from a folder with repertoire files.
#' In the second case you may want to provide a metadata file. In this case you should locate it in the folder.
#' It is necessary to name it exactly "metadata.txt".
#' If metadata is provided, \code{repLoad} will create a list in with 2 elements your environment, namely "data" and "meta".
#' Otherwise \code{repLoad} will create a list with the number of elements matching the number of your files.
#'
#' The metadata has to be a tab delimited file with first column named "Sample".
#' It can have any number of additional columns with arbitrary names.
#' The first column should contain base names of files without extensions in
#' your folder. Example:
#' \tabular{llll}{
#'  Sample \tab Sex \tab Age \tab Status\cr
#'  immunoseq_1 \tab M \tab 1 \tab C\cr
#'  immunoseq_2 \tab M \tab 2 \tab C\cr
#'  immunoseq_3 \tab F \tab 3 \tab A
#' }
#'
#' \code{repLoad} has the ".format" argument that sets the format for input repertoire files.
#' Do not pass it if you want immunarch to detect the formats and parse files automatically!
#' In case you want to force the package to parse the data in a specific format,
#' you can choose one of the several options for the argument:
#'
#' - "immunoseq" - ImmunoSEQ of any version. http://www.adaptivebiotech.com/immunoseq
#'
#' - "mitcr" - MiTCR. https://github.com/milaboratory/mitcr
#'
#' - "mixcr" - MiXCR of any version. https://github.com/milaboratory/mixcr
#'
#' - "migec" - MiGEC. http://migec.readthedocs.io/en/latest/
#'
#' - "migmap" - For parsing IgBLAST results postprocessed with MigMap. https://github.com/mikessh/migmap
#'
#' - "tcr" - tcR, our previous package. https://imminfo.github.io/tcr/
#'
#' - "vdjtools" - VDJtools of any version. http://vdjtools-doc.readthedocs.io/en/latest/
#'
#' - "imgt" - IMGT HighV-QUEST. http://www.imgt.org/HighV-QUEST/
#'
#' - "airr" - adaptive immune receptor repertoire (AIRR) data format. http://docs.airr-community.org/en/latest/datarep/overview.html
#'
#' - "10x" - 10XGenomics clonotype annotations tables. https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/output/annotation
#'
#' @seealso \link{immunr_data_format} for immunarch data format; \link{repSave} for file saving;
#' \link{repOverlap}, \link{geneUsage} and \link{repDiversity} for starting with immune repertoires basic statistics.
#'
#' @examples
#' \dontrun{
#' # To load the data from a single file without specifying the data format:
#' immdata <- repLoad("path/to/your/folder/immunoseq_1.txt")
#' # To load the data from a single ImmunoSEQ file, go with:
#' immdata <- repLoad("path/to/your/folder/immunoseq_1.txt", .format = "immunoseq")
#'
#' # Suppose you have a following structure in your folder:
#' # >_ ls
#' # immunoseq1.txt
#' # immunoseq2.txt
#' # immunoseq3.txt
#' # metadata.txt
#'
#' # To load the whole folder with every file in it type:
#' immdata <- repLoad("path/to/your/folder/")
#'
#' # We recommend creating a metadata file named exactly "metadata.txt" in the folder.
#'
#' # In that case, when you load your data you will see:
#' # > immdata <- repLoad("path/to/your/folder/")
#' # > names(immdata)
#' # [1] "data" "meta"
#'
#' # If you do not have "metadata.txt", you will see:
#' # > immdata <- repLoad("path/to/your/folder/")
#' # > names(immdata)
#' # [1] "immunoseq_1" "immunoseq_2" "immunoseq_3"
#' }
#'
#' @export repLoad
repLoad <- function (.path, .format = NA) {
  require(stringr)

  res = list()
  metadata = c()

  for (i in 1:length(.path)) {
    if (dir.exists(.path[i])) {
      cat("Parsing", .path[i], "...\n")
      file_list = list.files(.path[i], full.names = T)
      res <- c(res, repLoad(file_list, .format[1]))
    } else if (file.exists(.path[i])) {
      cat("Parsing", .path[i], "-- ")

      # Parsing metadata
      if (stringr::str_detect(.path[i], "metadata")) {
        cat("metadata\n")

        metadata = readr::read_tsv(.path[i], trim_ws = T, col_types = cols())

        to_drop = c()
        for (col_name in names(metadata)) {
          if (sum(is.na(metadata[[col_name]])) == nrow(metadata)) {
            to_drop = c(to_drop, col_name)
          }
        }

        if (length(to_drop)) {
          metadata = metadata[-match(to_drop, names(metadata))]
          cat("Dropping ", length(to_drop), "column(s) from metadata.txt. Do you have spaces or tabs after the name of the last column? You shouldn't have them.\n")
        }
      }
      # Parsing repertoire files
      else {
        # Detect format
        if (is.na(.format)) {
          cur_format = .detect_format(.path[i])
        } else {
          cur_format = .format[1]
        }

        if (!is.na(cur_format)) {
          cat(paste0(cur_format, "\n"))

          parse_fun = switch(cur_format,
                             mitcr = parse_mitcr,
                             migec = parse_migec,
                             migmap = parse_migmap,
                             mixcr = parse_mixcr,
                             tcr = parse_tcr,
                             vdjtools = parse_vdjtools,
                             imgt = parse_imgt,
                             immunoseq = parse_immunoseq,
                             airr = parse_airr,
                             `10x` = parse_10x,
                             immunarch = parse_immunarch,
                             NA)

          if (!suppressWarnings(is.na(parse_fun))) {
            parse_res = parse_fun(.path[i])
            if (has_class(parse_res, "list")) {
              res = c(res, parse_res)
            } else {
              res = c(res, list(parse_res))
              names(res)[length(res)] = .remove.ext(.path[i])
            }
          } else {
            cat("unsupported format, skipping\n")
          }
        } else {
          cat("unsupported format, skipping\n")
        }
      }
    } else {
      cat('Can\'t find folder or file:\t"', .path[i], '"', sep = '', end = '\n')
    }
  }

  if (length(metadata)) {
    list(data = res, meta = metadata)
  } else {
    if (length(res) == 1) {
      res[[1]]
    } else {
      res
    }
  }
}


#' Save repertoire files to the disk.
#'
#' @description The \code{repSave} function is deigned to save your data to the disk
#' in desirable format. Currently supports "immunarch" and "vdjtools" file formats.
#'
#' @param .data An R dataframe, a list of R dataframes or a list with \code{data} and
#' \code{meta} where first element is a list of dataframes and the latter is a dataframe
#' with metadata.
#' @param .path A string with the path to the output directory. It should include file
#' name if a single dataframe is provided to .data argument.
#' @param .format A string with desirable format specification. Current options are
#' "immunarch" and "vdjtools".
#' @param .compress A boolean value. Defines whether the output will be compressed or not.
#'
#' @details It is not necessary to create directories beforehand. If the provided directory
#' does not exist it will be created automatically.
#'
#' @export repSave
repSave <- function (.data, .path, .format = c("immunarch", "vdjtools"),
                     .compress = T) {
  .format = .format[1]
  if ( !.format %in% c("immunarch", "vdjtools") ) {
    stop("Unknown format. Please provide either 'immunarch' or 'vdjtools'")
  }

  if ( !has_class(.data, "list") ) {
    success = switch(.format,
                     immunarch = save_immunarch(.data, .path, .compress),
                     vdjtools = save_vdjtools(.data, .path, .compress),
                     NA)
    # readr::write_tsv(.data, path = .path)
  } else {
    if ( "meta" %in% names(.data) ) {
      ifelse(!dir.exists(.path), dir.create(.path), FALSE)
      readr::write_tsv(.data$meta, path = paste0(.path, "/metadata.txt"))
      for ( name in names(.data$data) ) {
        success = switch(.format,
                         immunarch = save_immunarch(.data$data[[name]],
                                                 paste0(.path, "/", name),
                                                 .compress),
                         vdjtools = save_vdjtools(.data$data[[name]],
                                                  paste0(.path, "/", name),
                                                  .compress),
                         NA)
        # readr::write_tsv(.data$data[[name]], path = paste0(.path, "/", name))
      }
    } else {
      ifelse(!dir.exists(.path), dir.create(.path), FALSE)
      for ( name in names(.data) ) {
        success = switch(.format,
                         immunarch = save_immunarch(.data[[name]],
                                                 paste0(.path, "/", name),
                                                 .compress),
                         vdjtools = save_vdjtools(.data[[name]],
                                                  paste0(.path, "/", name),
                                                  .compress),
                         NA)
        # readr::write_tsv(.data$data[[name]], path = paste0(.path, "/", name))
      }
    }
  }

}


##### IO utility functions #####

#' Split all-samples-in-one ImmunoSEQ files by sample
#'
#' @description Adaptive Biotechnologies ImmunoSEQ platform files are often downloaded as a big text file
#' with receptor data for all samples from the experiment.
#' Information about the source of each receptor is located in the "sample_name" column of the file.
#' Such format is impractical for using in R analysis pipelines and require preprocessing
#' to split the file into smaller files each corresponding to a specific sample.
#' The \code{split_immunoseq_files} function splits the file into smaller files and writes them in the same ImmunoSEQ
#' format to the output folder provided by the user.
#'
#' If you want to just load the big file as separated data frames into R to use with immunarch, just use the \link{repLoad} function
#' that automatically detectes the ImmunoSEQ format and splits the big file.
#'
#' @param .big.file Path to the input ImmunoSEQ file with the "sample_name" column with information about receptors' source sample.
#' @param .output.folder Path to the folder to write sample repertoires to.
#'
#' @seealso \link{repLoad}
#'
#' @examples
#' \dontrun{
#' # Split the file into smaller files
#' split_immunoseq_files("~/file/path/immunoseq_big.tsv", "~/file/path/experiment1")
#'
#' # If you just want to load all files as separate data frames to use with immunarch,
#' # just use the repLoad function
#' immdata = repLoad("~/file/path/immunoseq_big.tsv")
#' }
#' @export split_immunoseq_files
split_immunoseq_files <- function (.big.file, .output.folder) {
  cat("Reading", .big.file, "...\n")
  df = readr::read_tsv(.big.file)
  df_list = split(df, df$sample_name)
  cat("Writing files to", .output.folder, "...\n")
  for(name in names(df_list)) {
    file_name = stringr::str_c(.output.folder, name, ".txt", sep = "")
    cat("Output file:", file_name, "\n")
    readr::write_tsv(df_list[[name]], file_name)
  }
}

.remove.ext <- function (.str) {
  gsub(pattern = '.*/|[.].*$', replacement = '', x = .str)
}


.detect_format <- function (.filename) {
  res_format = NA

  f = file(.filename, "r")
  l = readLines(f, 1)
  close(f)

  if (any(str_detect(l, c("MiTCRFullExport", "mitcr")))) {
    res_format = "mitcr"
  } else if (str_detect(l, "CDR3 amino acid sequence") && str_detect(l, "V segment") && !str_detect(l, "Good events")) {
    res_format = "mitcr"
  } else if (str_detect(l, "CDR3 amino acid sequence") && str_detect(l, "V segment") && str_detect(l, "Good events")) {
    res_format = "migec"
  } else if (str_detect(l, "v.end.in.cdr3") && str_detect(l, "cdr3aa")) {
    res_format = "migmap"
  } else if (str_detect(l, "CDR3.amino.acid.sequence") && str_detect(l, "Umi.count")) {
    res_format = "tcr"
  } else if (str_detect(tolower(l), "cdr3nt") && str_detect(tolower(l), "vend") && str_detect(tolower(l), 'v')) {
    res_format = "vdjtools"
  } else if (str_detect(tolower(l), "count") && str_detect(tolower(l), "sequence") && str_detect(tolower(l), 'd segment')) {
    res_format = "vdjtools"
  } else if (str_detect(tolower(l), "junction start") && str_detect(tolower(l), "v-d-j-region end") && str_detect(tolower(l), 'v-region')) {
    res_format = "imgt"
  } else if (str_detect(tolower(l), "v_resolved") && str_detect(tolower(l), "amino_acid")) {
    res_format = "immunoseq"
  } else if (str_detect(tolower(l), "maxresolved")) {
      res_format = "immunoseq"
  } else if (str_detect(tolower(l), "v_gene") && str_detect(tolower(l), "templates") && str_detect(tolower(l), "amino_acid")) {
    res_format = "immunoseq"
  } else if (str_detect(tolower(l), "allvalignment") && str_detect(tolower(l), "vhit")) {
    res_format = "mixcr"
  } else if (str_detect(tolower(l), "clonal sequence")) {
    res_format = "mixcr"
  } else if (str_detect(tolower(l), "junction_aa") && str_detect(tolower(l), "cigar")) {
    res_format = "airr"
  } else if (str_detect(tolower(l), "clonotype_id") && str_detect(tolower(l), "v_gene")) {
    res_format = "10x"
  } else if (str_detect(tolower(l), "exported from immunarch")) {
    res_format = "immunarch"
  }

  res_format
}


#' Additional functions for help in parsing
#'
.make_names <- function (.char) {
  if (is.na(.char[1])) { NA }
  # else { tolower(make.names(.char)) }
  else { tolower(.char) }
}


.which_recomb_type <- function (.name) {

  recomb_type = NA

  i = 1

  while (is.na(recomb_type) && i < 100) {
    if (any(str_detect(.name[i], c("TCRA", "TRAV", "TCRG", "TRGV", "IGKV", "IGLV")))) {
      recomb_type = "VJ"
    } else if (any(str_detect(.name[i], c("TCRB", "TRBV", "TCRD", "TRDV", "IGHV")))) {
      recomb_type = "VDJ"
    }
    i = i + 1
  }
  if (is.na(recomb_type)) {
    warning("Can't determine the type of V(D)J recombination. No insertions will be presented in the resulting data table.")
  }

  recomb_type
}


.get_coltypes <- function (.filename, .nuc.seq, .aa.seq, .count,
                           .vgenes, .jgenes, .dgenes,
                           .vend, .jstart, .dstart, .dend,
                           .vd.insertions, .dj.insertions, .total.insertions,
                           .skip, .sep, .add = NA) {
  table.colnames <- colnames(readr::read_delim(.filename,
                                               col_types = cols(),
                                               delim = .sep,
                                               quote = "",
                                               escape_double = F,
                                               comment = "",
                                               n_max = 1,
                                               trim_ws = T,
                                               skip = .skip))

  swlist <- list(col_character(), col_character(),
                 col_integer(),
                 col_character(), col_character(), col_character(),
                 col_integer(), col_integer(), col_integer(), col_integer(),
                 col_integer(), col_integer(), col_integer())

  names(swlist) <- tolower(c(.nuc.seq, ifelse(is.na(.aa.seq), "NA", .aa.seq),
                             .count,
                             .vgenes, .jgenes, .dgenes,
                             .vend, .jstart, .dstart, .dend,
                             .vd.insertions, .dj.insertions, .total.insertions))
  if (!is.na(.add)) {
    swlist = c(swlist, rep(col_guess(), length(.add)))
    names(swlist)[tail(1:length(swlist), length(.add))] = .add
  }

  swlist <- c(swlist, '_')

  if (is.na(.aa.seq)) {
    swlist = swlist[-2]
  }

  col.classes <- list(sapply(tolower(table.colnames), function (x) {
    do.call(switch, c(x, swlist))
  }, USE.NAMES = F))[[1]]
  names(col.classes) = table.colnames

  col.classes
}

.remove.alleles <- function (.data) {
  if (has_class(.data, "list")) {
    lapply(.data, .remove.alleles)
  } else {
    .data[[IMMCOL$v]] = return_segments(.data[[IMMCOL$v]])
    .data[[IMMCOL$j]] = return_segments(.data[[IMMCOL$j]])
    .data
  }
}

.postprocess <- function (.data) {
  logic = is.na(.data[[IMMCOL$cdr3aa]]) & !is.na(.data[[IMMCOL$cdr3nt]])
  if (sum(logic)) {
    .data[[IMMCOL$cdr3aa]][logic] = bunch_translate(.data[[IMMCOL$cdr3nt]][logic])
  }

  for (colname in c(IMMCOL$ve, IMMCOL$ds, IMMCOL$de, IMMCOL$js, IMMCOL$vnj, IMMCOL$vnd, IMMCOL$dnj)) {
    if (colname %in% colnames(.data)) {
      logic = is.na(.data[[colname]])
      .data[[colname]][logic] = -1

      logic = .data[[colname]] < 0
      .data[[colname]][logic] = NA
    }
  }

  for (col_i in 1:length(IMMCOL$order)) {
    colname = IMMCOL$order[col_i]
    if (colname %in% colnames(.data)) {
      if (!has_class(.data[[colname]], IMMCOL$type[col_i])) {
        .data[[colname]] = as(.data[[colname]], IMMCOL$type[col_i])
      }
    }
  }

  .data
}


##### Parsers #####

parse_repertoire <- function (.filename, .nuc.seq, .aa.seq, .count,
                              .vgenes, .jgenes, .dgenes,
                              .vend, .jstart, .dstart, .dend,
                              .vd.insertions, .dj.insertions, .total.insertions,
                              .skip = 0, .sep = '\t', .add = NA) {
  .nuc.seq <- .make_names(.nuc.seq)
  .aa.seq <- .make_names(.aa.seq)
  .count <- .make_names(.count)
  .vgenes <- .make_names(.vgenes)
  .jgenes <- .make_names(.jgenes)
  .dgenes <- .make_names(.dgenes)
  .vend <- .make_names(.vend)
  .jstart <- .make_names(.jstart)
  .vd.insertions <- .make_names(.vd.insertions)
  .dj.insertions <- .make_names(.dj.insertions)
  .total.insertions <- .make_names(.total.insertions)
  .dstart <- .make_names(.dstart)
  .dend <- .make_names(.dend)
  .add = .make_names(.add)

  col.classes = .get_coltypes(.filename, .nuc.seq, .aa.seq, .count,
                              .vgenes, .jgenes, .dgenes,
                              .vend, .jstart, .dstart, .dend,
                              .vd.insertions, .dj.insertions, .total.insertions,
                              .skip = .skip, .sep = '\t', .add)

  suppressMessages(df <- readr::read_delim(.filename, col_names=T,
                                           col_types = col.classes, delim = .sep,
                                           quote = "", escape_double = F,
                                           comment = "", trim_ws = T,
                                           skip = .skip, na = c("", "NA", ".")))

  names(df) = tolower(names(df))
  recomb_type = .which_recomb_type(df[[.vgenes]])

  table.colnames = names(col.classes)

  df[[.nuc.seq]] = toupper(df[[.nuc.seq]])

  if(is.na(.aa.seq)) {
    df$CDR3.amino.acid.sequence = bunch_translate(df[[.nuc.seq]])
    .aa.seq = "CDR3.amino.acid.sequence"
  }

  if (is.na(.count)) {
    .count = "Count"
    df$Count = 1
  }
  df$Proportion = df[[.count]] / sum(df[[.count]])
  .prop = "Proportion"

  ins_ok = F
  if (is.na(.vd.insertions)) {
    .vd.insertions <- "VD.insertions"
    df$VD.insertions <- NA
  }

  if (!(.vd.insertions %in% table.colnames)) {
    .vd.insertions <- "VD.insertions"
    df$VD.insertions <- NA

    if (!is.na(.vend) && !is.na(.dstart)) {
      if (!is.na(recomb_type) && recomb_type == "VDJ") {
        df$VD.insertions <- df[[.dstart]] - df[[.vend]] - 1
        df$VD.insertions[is.na(df[[.dstart]])] <- NA
        df$VD.insertions[is.na(df[[.vend]])] <- NA

        ins_ok = T
      }
    }
  }

  if(!ins_ok) {
    df$V.end <- NA
    df$D.start <- NA
    df$D.end <- NA
    .vend <- "V.end"
    .dstart <- "D.start"
    .dend = "D.end"
  }

  ins_ok = F
  if (is.na(.dj.insertions)) {
    .dj.insertions <- "DJ.insertions"
    df$DJ.insertions <- NA
  }

  if (!(.dj.insertions %in% table.colnames)) {
    .dj.insertions <- "DJ.insertions"
    df$DJ.insertions <- NA

    if (!is.na(.jstart) && !is.na(.dend)) {
      if (!is.na(recomb_type) && recomb_type == "VDJ") {
        df$DJ.insertions <- df[[.jstart]] - df[[.dend]] - 1
        df$DJ.insertions[is.na(df[[.dend]])] <- NA
        df$DJ.insertions[is.na(df[[.jstart]])] <- NA

        ins_ok = T
      }
    }
  }
  if(!ins_ok) {
    df$J.start <- NA
    df$D.start <- NA
    df$D.end <- NA
    .jstart <- "J.start"
    .dstart <- "D.start"
    .dend = "D.end"
  }

  ins_ok = F
  if (is.na(.total.insertions)) {
    .total.insertions <- "Total.insertions"
    df$Total.insertions <- NA
  }

  if (!(.total.insertions %in% table.colnames)) {
    .total.insertions <- "Total.insertions"
    df$Total.insertions <- -1
    if (!is.na(recomb_type)) {
      if (recomb_type == "VJ") {
        df$Total.insertions <- df[[.jstart]] - df[[.vend]] - 1
        df$Total.insertions[df$Total.insertions < 0] <- 0
      } else if (recomb_type == "VDJ" ) {
        df$Total.insertions <- df[[.vd.insertions]] + df[[.dj.insertions]]
      }
    } else {
      df$Total.insertions = NA
    }
  }

  vec_names = c(.count, .prop, .nuc.seq, .aa.seq,
                .vgenes, .dgenes, .jgenes,
                .vend, .dstart, .dend, .jstart,
                .total.insertions, .vd.insertions, .dj.insertions)
  if (!is.na(.add)) { vec_names = c(vec_names, .add) }

  df <- df[, vec_names]

  colnames(df)[1] = IMMCOL$count
  colnames(df)[2] = IMMCOL$prop
  colnames(df)[3] = IMMCOL$cdr3nt
  colnames(df)[4] = IMMCOL$cdr3aa
  colnames(df)[5] = IMMCOL$v
  colnames(df)[6] = IMMCOL$d
  colnames(df)[7] = IMMCOL$j
  colnames(df)[8] = IMMCOL$ve
  colnames(df)[9] = IMMCOL$ds
  colnames(df)[10] = IMMCOL$de
  colnames(df)[11] = IMMCOL$js
  colnames(df)[12] = IMMCOL$vnj
  colnames(df)[13] = IMMCOL$vnd
  colnames(df)[14] = IMMCOL$dnj

  .postprocess(df)
}

parse_immunoseq <- function (.filename, .wash.alleles = T) {

  .fix.immunoseq.genes <- function (.col) {
    # fix ","
    .col <- gsub(",", ", ", .col, fixed = T, useBytes = T)
    # fix forward zeros
    .col <- gsub("-([0])([0-9])", "-\\2", .col, useBytes = T)
    .col <- gsub("([VDJ])([0])([0-9])", "\\1\\3", .col, useBytes = T)
    # fix gene names
    .col <- gsub("TCR", "TR", .col, fixed = T, useBytes = T)
    .col
  }

  filename <- .filename
  file_cols = list()
  file_cols[[IMMCOL$count]] = "templates"
  file_cols[[IMMCOL$cdr3nt]] = 'rearrangement'
  file_cols[[IMMCOL$cdr3aa]] = 'amino_acid'
  file_cols[[IMMCOL$v]] = 'v_resolved'
  file_cols[[IMMCOL$d]] = 'd_resolved'
  file_cols[[IMMCOL$j]] = 'j_resolved'
  file_cols[[IMMCOL$vnj]] = "n1_insertions"
  file_cols[[IMMCOL$vnd]] = "n1_insertions"
  file_cols[[IMMCOL$dnj]] = "n2_insertions"

  v_index_col_name = "v_index"
  d_index_col_name = "d_index"
  j_index_col_name = "j_index"
  n1_index_col_name = "n1_index"
  n2_index_col_name = "n2_index"

  #
  # Check for the version of ImmunoSEQ files
  #
  f = file(.filename, "r")
  l = readLines(f, 2)
  close(f)
  if (str_detect(l[[1]], "v_gene") && !str_detect(l[[1]], "v_resolved")) {
    file_cols[[IMMCOL$v]] = 'v_gene'
    file_cols[[IMMCOL$d]] = 'd_gene'
    file_cols[[IMMCOL$j]] = 'j_gene'
  } else if (str_detect(l[[1]], "MaxResolved")) {
    file_cols[[IMMCOL$v]] = 'vMaxResolved'
    file_cols[[IMMCOL$d]] = 'dMaxResolved'
    file_cols[[IMMCOL$j]] = 'jMaxResolved'

    file_cols[[IMMCOL$vnj]] = "n1insertion"
    file_cols[[IMMCOL$vnd]] = "n1insertion"
    file_cols[[IMMCOL$dnj]] = "n2insertion"
  }


  l_split = strsplit(l, "\t")
  if (str_detect(l[[1]], "templates")) {
    if (str_detect(l[[1]], "templates/reads")) {
      file_cols[[IMMCOL$count]] = "count (templates/reads)"
      file_cols[[IMMCOL$cdr3nt]] = 'nucleotide'
      file_cols[[IMMCOL$cdr3aa]] = 'aminoAcid'
    } else if (l_split[[2]][match("templates", l_split[[1]])] == "null") {
      file_cols[[IMMCOL$count]] = "reads"
    }
  }

  if (!str_detect(l[[1]], "v_index")) {
    v_index_col_name = "vindex"
    d_index_col_name = "dindex"
    j_index_col_name = "jindex"
    n1_index_col_name = "n1index"
    n2_index_col_name = "n2index"
  }

  for (col_name in names(file_cols)) {
    file_cols[[col_name]] = .make_names(file_cols[[col_name]])
  }

  file_cols[[IMMCOL$prop]] = IMMCOL$prop
  file_cols[[IMMCOL$ve]] = IMMCOL$ve
  file_cols[[IMMCOL$ds]] = IMMCOL$ds
  file_cols[[IMMCOL$de]] = IMMCOL$de
  file_cols[[IMMCOL$js]] = IMMCOL$js
  file_cols[[IMMCOL$seq]] = IMMCOL$seq

  suppressMessages(df <- readr::read_delim(.filename, col_names=T,
                    delim = "\t", quote = "",
                    escape_double = F, comment = "",
                    trim_ws = T, skip = 0))

  names(df) = tolower(names(df))

  df[[file_cols[[IMMCOL$prop]]]] = df[[file_cols[[IMMCOL$count]]]] / sum(df[[file_cols[[IMMCOL$count]]]])

  # Save full nuc sequences and cut them down to CDR3
  df[[IMMCOL$seq]] = df[[file_cols[[IMMCOL$cdr3nt]]]]

  # TODO: what if df[["v_index]] has "-1" or something like that?
  df[[file_cols[[IMMCOL$cdr3nt]]]] = stringr::str_sub(df[[IMMCOL$seq]], df[[v_index_col_name]]+1, nchar(df[[IMMCOL$seq]]))

  df[[file_cols[[IMMCOL$v]]]] = .fix.immunoseq.genes(df[[file_cols[[IMMCOL$v]]]])
  df[[file_cols[[IMMCOL$d]]]] = .fix.immunoseq.genes(df[[file_cols[[IMMCOL$d]]]])
  df[[file_cols[[IMMCOL$j]]]] = .fix.immunoseq.genes(df[[file_cols[[IMMCOL$j]]]])

  recomb_type = "VDJ"
  if (recomb_type == "VDJ") {
    df[[file_cols[[IMMCOL$ve]]]] = df[[n1_index_col_name]] - df[[v_index_col_name]]
    df[[file_cols[[IMMCOL$ds]]]] = df[[d_index_col_name]] - df[[v_index_col_name]]
    df[[file_cols[[IMMCOL$de]]]] = df[[n2_index_col_name]] - df[[v_index_col_name]]
    df[[file_cols[[IMMCOL$js]]]] = df[[j_index_col_name]] - df[[v_index_col_name]]
    file_cols[[IMMCOL$vnj]] = IMMCOL$vnj
    df[[IMMCOL$vnj]] = -1
  }

  sample_name_vec = NA
  if ("sample_name" %in% colnames(df)) {
    if (length(unique(df[["sample_name"]])) > 1) {
      sample_name_vec = df[["sample_name"]]
    }
  }

  df = df[unlist(file_cols[IMMCOL$order])]
  names(df) = IMMCOL$order

  if (.wash.alleles) {
    df = .remove.alleles(df)
    df[[IMMCOL$v]] = gsub("([VDJ][0-9]*)$", "\\1-1", df[[IMMCOL$v]], useBytes = T)
    df[[IMMCOL$j]] = gsub("([VDJ][0-9]*)$", "\\1-1", df[[IMMCOL$j]], useBytes = T)
  }

  if (nrow(df) > 0) {
    if (has_class(df[[IMMCOL$vnj]], "character")) {
      df[[IMMCOL$vnj]][df[[IMMCOL$vnj]] == "no data"] = NA
    }
    if (has_class(df[[IMMCOL$vnd]], "character")) {
      df[[IMMCOL$vnd]][df[[IMMCOL$vnd]] == "no data"] = NA
    }
    if (has_class(df[[IMMCOL$dnj]], "character")) {
      df[[IMMCOL$dnj]][df[[IMMCOL$dnj]] == "no data"] = NA
    }

    df[[IMMCOL$vnj]] = as.integer(df[[IMMCOL$vnj]])
    df[[IMMCOL$vnd]] = as.integer(df[[IMMCOL$vnd]])
    df[[IMMCOL$dnj]] = as.integer(df[[IMMCOL$dnj]])
  }

  if (!is.na(sample_name_vec[1])) {
    df = lapply(split(df, sample_name_vec), .postprocess)
  } else {
    .postprocess(df)
  }
}

parse_mitcr <- function (.filename) {
  .skip = 0
  f = file(.filename, "r")
  l = readLines(f, 1)
  # Check for different levels of the MiTCR output
  if (any(stringr::str_detect(l, c("MiTCRFullExport", "mitcr")))) { .skip = 1 }
  mitcr_format = 1
  if (stringr::str_detect(l, "MiTCRFullExport") || .skip == 0) { mitcr_format = 2 }
  close(f)

  if (mitcr_format == 1) {
    filename <- .filename
    .count = "count"
    nuc.seq <- 'cdr3nt'
    aa.seq <- "cdr3aa"
    vgenes <- 'v'
    jgenes <- 'j'
    dgenes <- 'd'
    vend <- 'VEnd'
    jstart <- 'JStart'
    dstart <- 'DStart'
    dend = 'DEnd'
    vd.insertions <- NA
    dj.insertions <- NA
    total.insertions <- NA
    .sep = '\t'
  } else {
    # Check if there are barcodes
    f = file(.filename, "r")
    l = readLines(f, 1 + .skip)[.skip + 1]
    barcodes = NA
    .count = 'Read count'
    if ("NNNs" %in% strsplit(l, "\t", T)[[1]]) {
      .count = "NNNs"
    }
    close(f)

    filename <- .filename
    nuc.seq <- 'CDR3 nucleotide sequence'
    aa.seq <- "CDR3 amino acid sequence"
    vgenes <- 'V segments'
    jgenes <- 'J segments'
    dgenes <- 'D segments'
    vend <- 'Last V nucleotide position'
    jstart <- 'First J nucleotide position'
    dstart <- 'First D nucleotide position'
    dend = 'Last D nucleotide position'
    vd.insertions <- 'VD insertions'
    dj.insertions <- 'DJ insertions'
    total.insertions <- 'Total insertions'
    .sep = '\t'
  }

  parse_repertoire(.filename = filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
                   .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
                   .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
                   .total.insertions = total.insertions, .skip = .skip, .sep = .sep)
}

parse_mixcr <- function (.filename) {
  fix.alleles <- function (.data) {
    .data[[IMMCOL$v]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$v]])
    .data[[IMMCOL$d]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$d]])
    .data[[IMMCOL$j]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$j]])
    .data
  }

  .filename <- .filename
  .nuc.seq <- 'nseqcdr3'
  .aa.seq <- 'aaseqcdr3'
  .count <- 'clonecount'
  .sep = '\t'
  .vend <- "allvalignments"
  .jstart <- "alljalignments"
  .dalignments <- "alldalignments"
  .vd.insertions <- "VD.insertions"
  .dj.insertions <- "DJ.insertions"
  .total.insertions <- "Total.insertions"

  table.colnames <- tolower(make.names(read.table(.filename, sep = .sep, skip = 0, nrows = 1, stringsAsFactors = F, strip.white = T, comment.char = "", quote = "")[1,]))
  table.colnames <- gsub(".", "", table.colnames, fixed = T)

  if (!("nseqcdr3" %in% table.colnames)) {
    if ('targetsequences' %in% table.colnames) {
      .nuc.seq <- 'targetsequences'
      .aa.seq <- 'targetsequences'
    } else {
      .nuc.seq <- "nseqcdr3"
      .aa.seq <- "aaseqcdr3"
    }
  }

  if (!("allvalignments" %in% table.colnames)) {
    .vend = "allvalignment"
  }
  if (!("alldalignments" %in% table.colnames)) {
    .dalignments = "alldalignment"
  }
  if (!("alljalignments" %in% table.colnames)) {
    .jstart = "alljalignment"
  }

  if ("bestvhit" %in% table.colnames) {
    .vgenes <- 'bestvhit'
  } else if ('allvhits' %in% table.colnames) {
    .vgenes <- 'allvhits'
  } else if ('vhits' %in% table.colnames) {
    .vgenes <- 'vhits'
  } else if ('allvhitswithscore' %in% table.colnames) {
    .vgenes <- 'allvhitswithscore'
  } else {
    cat("Error: can't find a column with V genes\n")
  }

  if ("bestjhit" %in% table.colnames) {
    .jgenes <- 'bestjhit'
  } else if ('alljhits' %in% table.colnames) {
    .jgenes <- 'alljhits'
  } else if ('jhits' %in% table.colnames) {
    .jgenes <- 'jhits'
  } else if ('alljhitswithscore' %in% table.colnames) {
    .jgenes <- 'alljhitswithscore'
  } else {
    cat("Error: can't find a column with J genes\n")
  }

  if ("bestdhit" %in% table.colnames) {
    .dgenes <- 'bestdhit'
  } else if ('alldhits' %in% table.colnames) {
    .dgenes <- 'alldhits'
  } else if ('dhits' %in% table.colnames) {
    .dgenes <- 'dhits'
  } else if ('alldhitswithscore' %in% table.colnames) {
    .dgenes <- 'alldhitswithscore'
  } else {
    cat("Error: can't find a column with D genes\n")
  }

  df <- read_delim(file = .filename, col_types = cols(),
                   delim = .sep, skip = 0, comment = "",
                   quote = "", escape_double = F, trim_ws = T)
  names(df) <- tolower(gsub(".", "", names(df), fixed = T))
  names(df) <- str_replace_all(names(df), " ", "")

  # check for VJ or VDJ recombination
  # VJ / VDJ / Undeterm
  recomb_type = "Undeterm"
  if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRA", "TRAV", "TRGV", "IGKV", "IGLV"))) {
    recomb_type = "VJ"
  } else if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRB", "TRBV", "TRDV", "IGHV"))) {
    recomb_type = "VDJ"
  }

  .vd.insertions <- "VD.insertions"
  df$VD.insertions <- -1
  if (recomb_type == "VJ") {
    df$VD.insertions <- -1
  } else if (recomb_type == "VDJ") {
    logic <- sapply(str_split(df[[.dalignments]], "|"), length) >= 4 &
      sapply(strsplit(df[[.vend]], "|", T, F, T), length) >= 5
    df$VD.insertions[logic] <-
      as.numeric(sapply(str_split(df[[.dalignments]][logic], "|"), "[[", 4)) -
      as.numeric(sapply(strsplit(df[[.vend]][logic], "|", T, F, T), "[[", 5)) - 1
  }

  .dj.insertions <- "DJ.insertions"
  df$DJ.insertions <- -1
  if (recomb_type == "VJ") {
    df$DJ.insertions <- -1
  } else if (recomb_type == "VDJ") {
    logic <- sapply(strsplit(df[[.jstart]], "|", T, F, T), length) >= 4 &
      sapply(str_split(df[[.dalignments]], "|"), length) >= 5
    df$DJ.insertions[logic] <-
      as.numeric(sapply(strsplit(df[[.jstart]][logic], "|", T, F, T), "[[", 4)) -
      as.numeric(sapply(str_split(df[[.dalignments]][logic], "|"), "[[", 5)) - 1
  }

  logic <- (sapply(strsplit(df[[.vend]], "|", T, F, T), length) > 4) & (sapply(strsplit(df[[.jstart]], "|", T, F, T), length) >= 4)
  .total.insertions <- "Total.insertions"
  if (recomb_type == "VJ") {
    df$Total.insertions <- NA
    if (length(which(logic)) > 0) {
      df$Total.insertions[logic] <-
        as.numeric(sapply(strsplit(df[[.jstart]][logic], "|", T, F, T), "[[", 4)) - as.numeric(sapply(strsplit(df[[.vend]][logic], "|", T, F, T), "[[", 5)) - 1
    }
  } else if (recomb_type == "VDJ") {
    df$Total.insertions <- df[[.vd.insertions]] + df[[.dj.insertions]]
  } else {
    df$Total.insertions <- NA
  }
  df$Total.insertions[df$Total.insertions < 0] <- -1

  df$V.end <- -1
  df$J.start <- -1
  df[[.vend]] = gsub(";", "", df[[.vend]], fixed = T)
  logic = sapply(strsplit(df[[.vend]], "|", T, F, T), length) >= 5
  df$V.end[logic] <- sapply(strsplit(df[[.vend]][logic], "|", T, F, T), "[[", 5)
  logic = sapply(strsplit(df[[.jstart]], "|", T, F, T), length) >= 4
  df$J.start[logic] <- sapply(strsplit(df[[.jstart]][logic], "|", T, F, T), "[[", 4)

  .vend <- "V.end"
  .jstart <- "J.start"

  logic <- sapply(str_split(df[[.dalignments]], "|"), length) >= 5
  df$D5.end <- -1
  df$D3.end <- -1
  df$D5.end[logic] <- sapply(str_split(df[[.dalignments]][logic], "|"), "[[", 4)
  df$D3.end[logic] <- sapply(str_split(df[[.dalignments]][logic], "|"), "[[", 5)
  .dalignments <- c('D5.end', 'D3.end')

  .freq = "Proportion"
  df$Proportion = df[[.count]] / sum(df[[.count]])

  if (.aa.seq == .nuc.seq) {
    .aa.seq = "CDRaaseq"
    df[[.aa.seq]] = bunch_translate(df[[.nuc.seq]])
  }

  df <- df[, make.names(c(.count, .freq,
                          .nuc.seq, .aa.seq,
                          .vgenes, .dgenes, .jgenes,
                          .vend, .dalignments, .jstart,
                          .total.insertions, .vd.insertions, .dj.insertions, .nuc.seq))]

  colnames(df) <- IMMCOL$order

  df[[IMMCOL$v]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.]*[[:digit:]]*[)])", "", df[[IMMCOL$v]])
  df[[IMMCOL$v]] <- gsub(",", ", ", df[[IMMCOL$v]])
  df[[IMMCOL$d]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.]*[[:digit:]]*[)])", "", df[[IMMCOL$d]])
  df[[IMMCOL$d]] <- gsub(",", ", ", df[[IMMCOL$d]])
  df[[IMMCOL$j]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.]*[[:digit:]]*[)])", "", df[[IMMCOL$j]])
  df[[IMMCOL$j]] <- gsub(",", ", ", df[[IMMCOL$j]])

  .postprocess(fix.alleles(df))
}

parse_migec <- function (.filename) {
  filename <- .filename
  nuc.seq <- 'CDR3 nucleotide sequence'
  aa.seq <- 'CDR3 amino acid sequence'
  .count <- 'Good events'
  vgenes <- 'V segments'
  jgenes <- 'J segments'
  dgenes <- 'D segments'
  vend <- 'Last V nucleotide position'
  jstart <- 'First J nucleotide position'
  dstart <- 'First D nucleotide position'
  dend = 'Last D nucleotide position'
  vd.insertions <- 'VD insertions'
  dj.insertions <- 'DJ insertions'
  total.insertions <- 'Total insertions'
  .skip = 0
  .sep = '\t'

  parse_repertoire(.filename = filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count,
                   .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
                   .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
                   .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
                   .total.insertions = total.insertions, .skip = .skip, .sep = .sep)
}

parse_migmap <- function (.filename) {
  filename <- .filename
  nuc.seq <- 'cdr3nt'
  aa.seq <- 'cdr3aa'
  .count <- 'count'
  vgenes <- 'v'
  jgenes <- 'j'
  dgenes <- 'd'
  vend <- 'v.end.in.cdr3'
  jstart <- 'j.start.in.cdr3'
  dstart <- 'd.start.in.cdr3'
  dend = 'd.end.in.cdr3'
  vd.insertions <- NA
  dj.insertions <- NA
  total.insertions <- NA
  .skip = 0
  .sep = '\t'

  parse_repertoire(.filename = filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
                   .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
                   .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
                   .total.insertions = total.insertions, .skip = .skip, .sep = .sep)
}

parse_tcr <- function (.filename) {
  f = file(.filename, "r")
  l = readLines(f, 2)[2]
  close(f)

  nuc.seq <- 'CDR3.nucleotide.sequence'
  aa.seq <- 'CDR3.amino.acid.sequence'
  .count <- 'Read.count'
  vgenes <- 'V.gene'
  jgenes <- 'J.gene'
  dgenes <- 'D.gene'
  vend <- 'V.end'
  jstart <- 'J.start'
  dstart <- 'D5.end'
  dend = 'D3.end'
  vd.insertions <- "VD.insertions"
  dj.insertions <- "DJ.insertions"
  total.insertions <- "Total.insertions"
  .skip = 0
  .sep = '\t'

  if (substr(l, 1, 2) != "NA") {
    .count <- 'Umi.count'
  }

  parse_repertoire(.filename = .filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
                   .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
                   .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
                   .total.insertions = total.insertions, .skip = .skip, .sep = .sep)
}

parse_vdjtools <- function (.filename) {
  skip = 0

  # Check for different VDJtools outputs
  f = file(.filename, "r")
  l = readLines(f, 1)
  close(f)

  .skip = 0
  .count = "count"
  filename <- .filename
  nuc.seq <- 'cdr3nt'
  aa.seq <- 'CDR3aa'
  vgenes <- 'V'
  jgenes <- 'J'
  dgenes <- 'D'
  vend <- 'Vend'
  jstart <- 'Jstart'
  dstart <- 'Dstart'
  dend <- 'Dend'
  vd.insertions <- NA
  dj.insertions <- NA
  total.insertions <- NA
  .sep = '\t'

  if (length(strsplit(l, "-", T)) > 0) {
    if (length(strsplit(l, "-", T)[[1]]) == 3) {
      if (strsplit(l, "-", T)[[1]][2] == "header") {
        .count <- "count"
        .skip <- 1
      }
    } else if (tolower(substr(l, 1, 2)) == "#s") {
      .count <- "#Seq. count"
      nuc.seq <- 'N Sequence'
      aa.seq <- 'AA Sequence'
      vgenes <- 'V segments'
      jgenes <- 'J segments'
      dgenes <- 'D segment'
      vend <- NA
      jstart <- NA
      dstart <- NA
      dend <- NA
    } else if (stringr::str_detect(l, "#")) {
      .count <- "X.count"
    } else {
      .count <- "count"
    }
  }

  parse_repertoire(.filename = filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
                   .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
                   .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
                   .total.insertions = total.insertions, .skip = .skip, .sep = .sep)
}

parse_imgt <- function (.filename) {

  .fix.imgt.alleles <- function (.col) {
    sapply(strsplit(.col, " "), function (x) { if (length(x) > 1) {x[[2]]} else { NA } } )
  }

  f = file(.filename, "r")
  l = readLines(f, 2)[2]
  close(f)

  nuc.seq <- 'JUNCTION'
  aa.seq <- NA
  .count <- NA
  vgenes <- 'V-GENE and allele'
  jgenes <- 'J-GENE and allele'
  dgenes <- 'D-GENE and allele'
  vend <- "3'V-REGION end"
  jstart <- "5'J-REGION start"
  dstart <- "D-REGION start"
  dend <- "D-REGION end"
  vd.insertions <- NA
  dj.insertions <- NA
  total.insertions <- NA
  .skip = 0
  .sep = '\t'
  junc_start = .make_names("JUNCTION start")

  df = parse_repertoire(.filename = .filename, .nuc.seq = nuc.seq, .aa.seq = aa.seq, .count = .count, .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes,
                   .vend = vend, .jstart = jstart, .dstart = dstart, .dend = dend,
                   .vd.insertions = vd.insertions, .dj.insertions = dj.insertions,
                   .total.insertions = total.insertions, .skip = .skip, .sep = .sep, .add = junc_start)

  df[[IMMCOL$ve]] = df[[IMMCOL$ve]] - df[[junc_start]]
  df[[IMMCOL$ds]] = df[[IMMCOL$ds]] - df[[junc_start]]
  df[[IMMCOL$de]] = df[[IMMCOL$de]] - df[[junc_start]]
  df[[IMMCOL$js]] = df[[IMMCOL$js]] - df[[junc_start]]

  df[[IMMCOL$ve]][df[[IMMCOL$ve]] < 0] = NA
  df[[IMMCOL$ds]][df[[IMMCOL$ds]] < 0] = NA
  df[[IMMCOL$de]][df[[IMMCOL$de]] < 0] = NA
  df[[IMMCOL$js]][df[[IMMCOL$js]] < 0] = NA

  df[[IMMCOL$v]] = .fix.imgt.alleles(df[[IMMCOL$v]])
  df[[IMMCOL$d]] = .fix.imgt.alleles(df[[IMMCOL$d]])
  df[[IMMCOL$j]] = .fix.imgt.alleles(df[[IMMCOL$j]])

  df[[junc_start]] = NULL

  df
}

# parse_vidjil <- function (.filename) {
#   stop(IMMUNR_ERROR_NOT_IMPL)
# }
#
# parse_rtcr <- function (.filename) {
#   stop(IMMUNR_ERROR_NOT_IMPL)
# }
#
# parse_imseq <- function (.filename) {
#   stop(IMMUNR_ERROR_NOT_IMPL)
# }

parse_airr <- function (.filename) {
  df <- airr::read_rearrangement(.filename)

  df = df %>%
    select(sequence, v_call, d_call, j_call, junction, junction_aa,
           contains("v_germline_end"), contains("d_germline_start"), contains("d_germline_end"),
           contains("j_germline_start"), contains("np1_length"), contains("np2_length"),
           contains("duplicate_count"))

  namekey = c(duplicate_count = IMMCOL$count, junction = IMMCOL$cdr3nt, junction_aa = IMMCOL$cdr3aa,
              v_call = IMMCOL$v, d_call = IMMCOL$d, j_call = IMMCOL$j, v_germline_end = IMMCOL$ve,
              d_germline_start = IMMCOL$ds, d_germline_end = IMMCOL$de, j_germline_start = IMMCOL$js,
              np1_length = "unidins", np2_length = IMMCOL$dnj, sequence = IMMCOL$seq)

  names(df) <- namekey[names(df)]

  if (! ("unidins" %in% colnames(df))) {
    df["unidins"] = NA
  }

  recomb_type = .which_recomb_type(df[[IMMCOL$v]])

  if (!is.na(recomb_type)) {
    if (recomb_type == "VJ") {
      df[IMMCOL$vnj] = df["unidins"]
      df[IMMCOL$vnd] = NA
      df[IMMCOL$dnj] = NA
    } else if (recomb_type == "VDJ") {
      df[IMMCOL$vnj] = NA
      df[IMMCOL$vnd] = df["unidins"]
    }
  }

  for (column in  IMMCOL$order) {
    if (! (column %in% colnames(df))) {
      df[column] = NA
    }
  }

  df = df[IMMCOL$order]
  total = sum(df$Clones)
  df[IMMCOL$prop] = df[IMMCOL$count] / total
  df[IMMCOL$seq] = stringr::str_remove_all(df[[IMMCOL$seq]], "N")
  df = .postprocess(df)
  df
}

parse_immunarch <- function (.filename) {
  df = readr::read_tsv(.filename, skip = 1)
  df = .postprocess(df)
  df
}

parse_10x <- function (.filename) {
  parse_repertoire(.filename, .nuc.seq = "cdr3_nt", .aa.seq = NA, .count = "umis",
                   .vgenes = "v_gene", .jgenes = "j_gene", .dgenes = "d_gene",
                   .vend = NA, .jstart = NA, .dstart = NA, .dend = NA,
                   .vd.insertions = NA, .dj.insertions = NA, .total.insertions = NA,
                   .skip = 0, .sep = ',', .add = "chain")

}

##### Savers #####

save_immunarch <- function(.data, .path, .compress = T) {
  if ( .compress ) {
    filepath = gzfile(paste0(.path, ".tsv.gz"), compression = 9)
  } else {
    filepath = paste0(.path, ".tsv")
  }
  readr::write_lines(paste0("# Exported from immunarch ", packageVersion("immunarch") , " https://immunarch.com"), path = filepath)
  readr::write_tsv(x = .data, path = filepath, append = T, col_names = T)
}

save_vdjtools <- function(.data, .path, .compress = T) {
  if ( .compress ) {
    filepath = gzfile(paste0(.path, ".tsv.gz"), compression = 9)
  } else {
    filepath = paste0(.path, ".tsv")
  }

  old = c(IMMCOL$count,
          IMMCOL$prop,
          IMMCOL$cdr3nt,
          IMMCOL$cdr3aa,
          IMMCOL$v,
          IMMCOL$d,
          IMMCOL$j)

  new = c("#Seq. Count",
          "Percent",
          "N Sequence",
          "AA Sequence",
          "V Segments",
          "D Segment",
          "J Segments")

  names(.data) <- plyr::mapvalues(names(.data),
                                  from = old,
                                  to = as.character(new))

  readr::write_tsv(x = .data, path = filepath)
}

