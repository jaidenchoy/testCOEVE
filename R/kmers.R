#' Kmer analysis of repertoires.
#'
#' @aliases getKmers get.kmers makeKmerTable
#'
#' @export getKmers get.kmers
getKmers <- function (.data, .k, .col = c("aa", "nuc")) {
  seq_col = switch(.col[1], nuc = "CDR3.nucleotide.sequence", aa = "CDR3.amino.acid.sequence", stop("Wrong column name"))

  if (has_class(.data, "list")) {
    res = data.table(split_to_kmers(collect(select_(.data[[1]], seq_col), n = Inf)[[1]], .k = .k))
    colnames(res)[2] = names(.data)[1]
    for (i in 2:length(.data)) {
      tmp = data.table(split_to_kmers(collect(select_(.data[[i]], seq_col), n = Inf)[[1]], .k = .k))
      res = merge(res, tmp, by = "Kmer", all = T, suffixes = c("", ""))
      colnames(res)[ncol(res)] = names(.data)[i]
    }
  } else {
    res = split_to_kmers(collect(select_(.data, seq_col), n = Inf)[[1]], .k = .k)
  }

  add_class(tbl_df(as.data.frame(res)), "immunr_kmer_table")
}

get.kmers = getKmers


#' WIP
#' @export kmerAnalysis kmer.an kmerAnalysis.default kmerAnalysis.immunr_kmer_table
kmerAnalysis <- function (.data) {
  UseMethod("kmerAnalysis")
}

kmer.an = kmerAnalysis

kmerAnalysis.default <- function (.data, .k, .method = c("profile", "gibbs"), .col = c("aa", "nuc"), .ic = F, .remove.stop = T, .use.counts = F) {
  kmers = getKmers(.data, .k = .k, .col = .col)

  kmerAnalysis(kmers, .method = .method, .use.counts = .use.counts, .ic = .ic, .remove.stop = .remove.stop)
}

kmerAnalysis.immunr_kmer_table <- function (.data, .method = c("profile", "gibbs"), .ic = F, .remove.stop = T, .use.counts = F) {
  .method = .method[1]

  if (.method == "profile") {

  } else if (.method == "gibbs") {
    stop("Unimplemented")
  } else {
    stop("Unknown method")
  }

  res
}


#' Kmers processing.
#'
#' @aliases split_to_kmers kmer_profile
#'
#' @param .data Character vector or the output from \code{getKmers}.
#' @param .ic If T then return the position-weight matrix with self-information instead of frequencies.
#'
#' @export split_to_kmers kmer_profile
split_to_kmers <- function (.data, .k) {
  max_len = max(nchar(.data))
  tab = table(unlist(lapply(1:(max_len - .k + 1), function (i) substr(.data, i, i + .k - 1))))
  tab = tab[nchar(names(tab)) == .k]
  add_class(tbl_df(data.frame(Kmer = names(tab), Count = tab, stringsAsFactors = F)), "immunr_kmers")
}

kmer_profile <- function (.data, .ic = F, .remove.stop = T, .norm = T) {
  if (has_class(.data, "immunr_kmers")) {
    seq_vec = .data$Kmer
    cnt_vec = .data$Count
  } else if (length(table(nchar(.data))) > 1) {
    stop("Not all kmers in the input data have the same length.")
  } else {
    seq_vec = .data
    cnt_vec = rep.int(1, length(.data))
  }

  k = nchar(seq_vec[1])
  aas = sort(unique(AA_TABLE))
  if (.remove.stop) {
    aas = aas[3:length(aas)]
    cnt_vec = cnt_vec[grep("[*~]", seq_vec, invert = T)]
    seq_vec = seq_vec[grep("[*~]", seq_vec, invert = T)]
  }

  res = matrix(0, length(aas), k)
  row.names(res) = aas
  for (i in 1:k) {
    tab = tapply(cnt_vec, substr(seq_vec, i, i), sum, simplify = F)
    for (aa in names(tab)) {
      res[aa, i] = res[aa, i] + tab[[aa]]
    }
    if (.norm) { res[,i] = res[,i] / sum(res[,i]) }
    if (.ic) { res[,i] = -res[,i] * log2(res[,i])}
  }

  add_class(res, "immunr_kmer_profile")
}


#' Simple Gibbs sampling on kmers.
#'
#' @export
gibbs_sampling <- function (.data, .motif.len = 5, .niter = 500) {
  .score <- function (.seq, .i, .prof, .background) {
    kmer_aa = strsplit(substr(.seq, seq_i, seq_i + .motif.len - 1), "")[[1]]
    prod(sapply(1:.motif.len, function (kmer_pos) {
      sc = .prof[kmer_aa[kmer_pos], kmer_pos] / .background[kmer_aa[kmer_pos]]
      if (is.nan(sc)) { sc = 0 }
      sc
    }))
  }

  cat("Removed", sum(nchar(.data) < .motif.len), "sequences with the length less than the length of motifs.\n")
  seq_vec = .data[nchar(.data) >= .motif.len]
  background = table(unlist(strsplit(seq_vec, "")))
  background = background / sum(background)

  # Vector of scores for each position in the each input sequence
  score_vec = lapply(seq_vec, function (seq_x) rep(1, nchar(seq_x) - .motif.len + 1) )
  start_pos = sapply(nchar(seq_vec), function (max_pos) sample(1:(max_pos - .motif.len + 1), 1))

  # In the loop:
  pb = set_pb(.niter)
  for (iter in 1:.niter) {
    # Get random kmers
    prev_start_pos = start_pos
    start_pos = sapply(nchar(seq_vec), function (max_pos) sample(1:(max_pos - .motif.len + 1), 1))

    for (out_kmer_i in sample(1:length(seq_vec), length(seq_vec))) {
      max_pos = nchar(seq_vec[out_kmer_i]) - .motif.len + 1
      kmers <- substr(seq_vec[-out_kmer_i], start_pos[-out_kmer_i], start_pos[-out_kmer_i] + .motif.len - 1)
      prof = kmer_profile(kmers[-out_kmer_i])
      for (seq_i in 1:max_pos) {
        score_vec[[out_kmer_i]][seq_i] = .score(seq_vec[out_kmer_i], seq_i, prof, background)
      }
      if (sum(score_vec[[out_kmer_i]]) != 0) {
        poses = c(1:max_pos)[!is.na(score_vec[[out_kmer_i]])]
        start_pos[out_kmer_i] = sample(c(1:max_pos), 1, prob = score_vec[[out_kmer_i]][poses] / sum(score_vec[[out_kmer_i]][poses]))
      }
    }

    add_pb(pb)

    if (sum(prev_start_pos != start_pos) == 0) {
      break
    }
  }
  close(pb)

  data.frame(Motif = substr(seq_vec, start_pos, start_pos + .motif.len - 1),
             Start = start_pos,
             Score = sapply(1:length(score_vec), function (i) { score_vec[[i]][start_pos[i]] }), stringsAsFactors = F)
}
