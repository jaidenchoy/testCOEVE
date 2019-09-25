#' Gene usage
#'
#' @aliases geneUsage get_aliases get_genes
#'
#' @description
#' An utility function to analyse the immune receptor gene usage
#' (IGHD, IGHJ, IDHV, IGIJ, IGKJ, IGKV, IGLJ, IGLV, TRAJ, TRAV, TRBD, etc.)
#' and statistics. For gene details run \code{gene_stats()}.
#'
#' @usage
#' geneUsage(.data, .gene = c("hs.trbv", "HomoSapiens.TRBJ", "macmul.IGHV"),
#' .quant = c(NA, "count"), .ambig = c("exc", "inc", "wei", "maj"),
#' .type = c("segment", "allele", "family"), .norm = F)
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
#' @param .gene A character vector of length one with the name of the gene you want
#' to analyse of the specific species. If you provide a vector of different length, only first element
#' will be used. The string should also contain the species of interest, for example, valid ".gene" arguments
#' are "hs.trbv", "HomoSapiens.TRBJ" or "macmul.IGHV". For details run \code{gene_stats()}.
#' @param .quant Select the column with data to evaluate.
#' Pass NA if you want to compute gene statistics at the clonotype level without re-weighting.
#' Pass "count" to use the "Clones" column to weight genes by abundance of their corresponding clonotypes.
#' @param .ambig An option to handle ambiguous data. You can exclude data for the cases where
#' there is no clear match for gene, include it for every supplied gene, include it with weights,
#' or pick only first from the set. Set it to "exc", "inc", "wei", or "maj", correspondingly.
#' @param .type Set the type of data to evaluate: "segment", "allele", or "family".
#' @param .norm If TRUE than return proportions of genes. If FALSE then return counts of genes.
#'
#' @export geneUsage gen.us
geneUsage <- function (.data, # df, list, MonetDB
                       .gene = c("hs.trbv", "HomoSapiens.TRBJ", "macmul.IGHV"),
                       .quant = c(NA, "count"),
                       .ambig = c("inc", "exc", "wei", "maj"),
                       .type = c("segment", "allele", "family"),
                       .norm = F) {
  .type = .type[1]
  .ambig = .ambig[1]
  .quant = .quant[1]
  .gene = .gene[1]

  if (has_class(.data, 'list')) {
    if (length(.data) == 1) {
      return(geneUsage(.data[[1]]))
    } else {
      genedata = lapply(.data, geneUsage, .gene = .gene, .quant = .quant,
                        .ambig = .ambig, .type = .type, .norm = .norm)
      # print(names(.data))
      # print(names(genedata))
      # Add check if list has 1 element
      for(i in 2:(length(genedata))){
        if (i>2){
          allres = full_join(allres, genedata[[i]], by = "Names")
        } else{
          allres = full_join(genedata[[i-1]], genedata[[i]], by = "Names")
        }
      }
      colnames(allres) = c("Names", names(genedata))

      return(add_class(allres, "immunr_gene_usage"))
    }
  }

  '%!in%' <- function(x,y)!('%in%'(x,y))

  which_gene = strsplit(.gene, ".", T)[[1]][2]
  .gene = get_genes(.gene, .type)

  gene_col = tolower(substr(which_gene, 4, 4))
  if (gene_col == "v") {
    gene_col = IMMCOL$v
  } else if (gene_col == "d") {
    gene_col = IMMCOL$d
  } else if (gene_col == "j") {
    gene_col = IMMCOL$j
  } else {
    stop("The entered gene_col name is invalid")
  }

  if (is.na(.quant)) {
    condition = "n()"
    .quant = IMMCOL$count
  } else {
    condition = paste0("sum(", IMMCOL$count, ")", collapse = "")
  }

  if (.ambig == 'inc') {
    if (.type == "segment") {
      .data[[gene_col]] = return_segments(.data[[gene_col]])
    } else if (.type == "family") {
      .data[[gene_col]] = return_families(.data[[gene_col]])
    }
  }

  dataset = .data %>%
    select_(Gene = gene_col, .quant) %>%
    group_by(Gene) %>%
    summarise_(count = condition) %>%
    collect()

  if (.ambig == 'inc') {
    res = dataset
    names(res) = c("Names", IMMCOL$count)
  } else {
    gene_col_split = strsplit(dataset[["Gene"]], ", ", fixed = T, useBytes = T)
    counts = dataset$count
    names(counts) = dataset[["Gene"]]
    if (.ambig == 'exc'){
      included_names = names(counts) %in% gene_col_split[sapply(gene_col_split, length) == 1]
      chosen = counts[included_names]
      chosen = chosen[names(chosen) %in% .gene]
      other = counts[names(chosen) %!in% .gene]
      res = data.frame(Names = names(chosen),
                       Count = unlist(chosen), stringsAsFactors = F)
    } else if (.ambig == 'maj'){
      maj_counts = data.frame(Names = sapply(gene_col_split, "[[", 1), Counts = counts, stringsAsFactors = FALSE)
      row.names(maj_counts) = NULL
      chosen = filter(maj_counts, Names %in% .gene)
      other = filter(maj_counts, Names %!in% .gene)
      res = group_by(chosen, Names) %>%
        summarise(Counts = sum(Counts))
      Names = "Ambigous"
      Counts = sum(other$Counts)
      Ambigous = data.frame(Names, Counts)
      res = rbind(res, Ambigous)
    } else if (.ambig == 'wei'){
        wei_counts = data.frame(Names = names(counts), Counts = counts/sapply(gene_col_split, length), stringsAsFactors = FALSE)
        chosen = filter(wei_counts, Names %in% .gene)
        other = filter(wei_counts, Names %!in% .gene)
        res = chosen %>%
          mutate(Names = strsplit(as.character(Names), ",")) %>%
          unnest(Names)
        res = group_by(res, Names) %>%
          summarise(Counts = sum(Counts))
        Names = "Ambigous"
        Counts = sum(other$Counts)
        Ambigous = data.frame(Names, Counts)
        res = rbind(res, Ambigous)
      }
  }
  if (length(.gene %in% res[[1]]) < length(.gene)){
    missing = data.frame(.gene[! .gene %in% res[[1]]],
      0, stringsAsFactors = F)
    names(missing) = names(res)
    res = rbind(res, missing)
  }
  row.names(res) = NULL
  f = sort(res[[1]])
  res = res[match(f, res[[1]]), ]
  res = res[c(2:nrow(res), 1), ]
  if (.norm) {
    res[,2] = res[,2] / sum(res[,2])
  }
  add_class(res, "immunr_gene_usage")
}

#' WIP
#' @export gene_stats get_genes
gene_stats <- function () {
  res = GENE_SEGMENTS %>% group_by(alias, species, gene) %>% summarise(n = n()) %>% reshape2::dcast(alias+species ~ gene, value.var = "n")
  res[is.na(res)] = 0
  res
}

get_genes <- function (.gene = c("hs.trbv", "HomoSapiens.trbj", "macmul.ighv"), .type = c("segment", "allele", "family")) {
  .gene = .gene[1]
  .type = .type[1]

  # Species
  .gene = tolower(.gene)
  which_species = strsplit(.gene, ".", T)[[1]][1]
  if (which_species %in% GENE_SEGMENTS$alias) {
    which_species_col = "alias"
  } else {
    if (which_species %in% tolower(GENE_SEGMENTS$species)) {
      which_species_col = "species"
    } else {
      stop("Unknown species name / alias")
    }
  }

  # Genes
  which_gene = strsplit(.gene, ".", T)[[1]][2]
  if (!(which_gene %in% GENE_SEGMENTS[["gene"]][tolower(GENE_SEGMENTS[[which_species_col]]) == which_species])) {
    stop("Unknown gene name")
  }

  # Id's type
  which_type = paste0(.type, "_id")

  sort(unique(GENE_SEGMENTS[[which_type]][tolower(GENE_SEGMENTS[[which_species_col]]) == which_species & GENE_SEGMENTS[["gene"]] == which_gene]))
}

gen.us <- geneUsage
