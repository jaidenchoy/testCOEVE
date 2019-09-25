##### Utility functions #####


.rem_legend <- function (.p) { .p + theme(legend.position="none") }


.rem_xlab <- function (.p) { .p + theme(axis.title.x = element_blank()) }


.rem_ylab <- function (.p) { .p + theme(axis.title.y = element_blank()) }


.colourblind_vector <- function() {
  c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")
}


.colourblind_discrete <- function (.n, .colour = F) {
  cs <- .colourblind_vector()
  if (.colour) {
    scale_colour_manual(values = colorRampPalette(cs)(.n))
  } else {
    scale_fill_manual(values = colorRampPalette(cs)(.n))
  }
}


.tweak_fill <- function (.n) {
  palette_name = ""
  if (.n == 2) { palette_name = "Set1" }
  else if (.n < 4) { palette_name = "YlGnBu" }
  else if (.n < 6) {  palette_name = "RdBu" }
  else if (.n < 12) { palette_name = "Spectral" }
  else { return(scale_fill_hue()) }

  scale_fill_brewer(palette = palette_name)
}

.tweak_col <- function (.n) {
  palette_name = ""
  if (.n == 2) { palette_name = "Set1" }
  else if (.n < 4) { palette_name = "YlGnBu" }
  else if (.n < 6) {  palette_name = "RdBu" }
  else if (.n < 12) { palette_name = "Spectral" }
  else { return(scale_colour_hue()) }

  scale_color_brewer(palette = palette_name)
}


##### The one and only - the ultimate vis() function #####


#' One function to visualise them all
#'
#' @name vis
#'
#' @import ggplot2
#' @importFrom factoextra fviz_cluster fviz_dend fviz_pca_ind
#' @importFrom grDevices colorRampPalette
#'
#' @description Output from every function in immunarch can be visualised with a
#' single function — \code{vis}. The \code{vis} automatically detects
#' the type of the data and draws a proper visualisation. For example, output
#' from the \code{repOverlap} function will be identified as repertoire overlap values
#' and respective visualisation will be chosen without any additional arguments.
#' See "Details" for the list of available visualisations.
#'
#' @param .data Pass the output from any immunarch analysis tool to \code{vis()}.
#'
#' @details
#' List of available visualisations for different kinds of data.
#'
#' Basic analysis:
#'
#' - Exploratory analysis results (from \link{repExplore}) — see \link{vis.immunr_exp_vol};
#'
#' - Clonality statistics (from \link{repClonality}) — see \link{vis.immunr_homeo}.
#'
#' Overlaps and public clonotypes:
#'
#' - Overlaps (from \link{repOverlap}) using heatmaps, circos plots, polar area plots — see \link{vis.immunr_ov_matrix};
#'
#'-  Overlap clustering (from \link{repOverlapAnalysis}) — see \link{vis.immunr_hclust};
#'
#' - Repertoire incremental overlaps (from \link{repOverlap}) — see \link{vis.immunr_top_overlap};
#'
#' - Public repertoire abundance (from \link{pubRep}) — vis \link{vis.immunr_public_repertoire}.
#'
#' Gene usage:
#'
#' - Gene usage statistics (from \link{geneUsage}) using bar plots, box plots, treemaps — see \link{vis.immunr_gene_usage};
#'
#' - Gene usage distances (from \link{geneUsageAnalysis}) using heatmaps, circos plots, polar area plots — see \link{vis.immunr_ov_matrix};
#'
#' - Gene usage clustering (from \link{geneUsageAnalysis}) — see \link{vis.immunr_hclust}.
#'
#' Diversity estimation:
#'
#' - Diversity estimations (from \link{repDiversity}) — see \link{vis.immunr_chao1}.
#'
#' Advanced analysis:
#'
#' - Repertoire dynamics (from \link{trackClonotypes}) — see \link{vis.immunr_dynamics};
#'
#' - Sequence logo plots of amino acid distributions (from \link{kmer_profile}) — see \link{vis_logo};
#'
#' - Kmers distributions (from \link{getKmers}) — see \link{vis.immunr_kmers};
#'
#' - Mutation networks (from \link{mutationNetwork}) — see \link{vis.immunr_mutation_network};
#'
#' - CDR3 amino acid properties, e.g., biophysical (from \link{cdrProp}) — see \link{vis.immunr_cdr_prop}.
#'
#' Additionaly, we provide a wrapper functions for visualisations of common data types:
#'
#' - Any data frames or matrices using heatmaps — see \link{vis_heatmap} and \link{vis_heatmap2};
#'
#' - Any data frames or matrices using circos plots — see \link{vis_circos};
#'
#' - Any data frames or matrices using polar area plots — see \link{vis_radar}.
#'
#' @seealso \link{fixVis} for precise manipulation of plots.
#'
#' @examples
#' \dontrun{
#' # Load the test data
#' data(immdata)
#'
#' # Compute and visualise overlaps
#' ov = repOverlap(immdata$data)
#' vis(ov)
#' }
#'
#' @export vis
vis <- function (.data, ...) {
  UseMethod('vis')
}


##### Overlap & heatmaps, circos plots and polar area plots #####


#' Overlap / gene usage distance visualisation
#'
#' @name vis.immunr_ov_matrix
#'
#' @aliases vis.immunr_ov_matrix vis.immunr_gu_matrix
#'
#' @description Visualise matrices with overlap values or gene usage distances among samples.
#' For details see links below.
#'
#' @param .data Output from \link{repOverlap} or \link{geneUsageAnalysis}.
#'
#' @param .plot A string specifying the plot type:
#'
#' - "heatmap" for heatmaps using \link{vis_heatmap};
#'
#' - "heatmap2" for heatmaps using \link{vis_heatmap2};
#'
#' - "circos" for circos plots using \link{vis_circos};
#'
#' - "radar" for polar area plots using \link{vis_radar}.
#'
#' @param ... Other arguments are passed through to the underlying plotting function:
#'
#' - "heatmap" - passes arguments to \link{vis_heatmap};
#'
#' - "heatmap2" - passes arguments to \link{vis_heatmap2} and \link{heatmap} from the "heatmap3" package;
#'
#' - "circos" - passes arguments to \link{vis_circos} and \link{chordDiagram} from the "circlize" package;
#'
#' - "radar"- passes arguments to \link{vis_radar}.
#'
#' @export vis.immunr_ov_matrix vis.immunr_gu_matrix
vis.immunr_ov_matrix <- function (.data, .plot = c("heatmap", "heatmap2", "circos", "radar"), ...) {
  args = list(...)
  if (!(".title" %in% names(args))) {
    args$.title = "Repertoire overlap"
  }

  .plot = .plot[1]
  heatfun = vis_heatmap
  if (.plot == "heatmap2") {
    heatfun = vis_heatmap2
  } else if (.plot == "circos") {
    heatfun = vis_circos
  } else if (.plot == "radar") {
    heatfun = vis_radar
  }
  do.call(heatfun, c(list(.data = .data), args))
}

vis.immunr_gu_matrix <- function (.data, .plot = c("heatmap", "heatmap2", "circos", "radar"), ...) {
  args = list(...)
  if (!(".title" %in% names(args))) {
    args$.title = "Gene usage"
  }

  .plot = .plot[1]
  heatfun = vis_heatmap
  if (.plot == "heatmap2") {
    heatfun = vis_heatmap2
  } else if (.plot == "circos") {
    heatfun = vis_circos
  } else if (.plot == "radar") {
    heatfun = vis_radar
  }
  do.call(heatfun, c(list(.data = .data), args))
}


#' Visualisation of matrices and data frames using ggplot2-based heatmaps
#'
#' @name vis_heatmap
#'
#' @aliases vis_heatmap
#'
#' @description Fast and easy visualisations of matrices or data frames
#' with functions based on the ggplot2 package.
#'
#' @param .data Input object: a matrix or a data frame.
#'
#' If matrix: column names and row names (if presented) will be used as names for labs.
#'
#' If data frame: the first column will be used for row names and removed from the data.
#' Other columns will be used for values in the heatmap.
#'
#' @param .text If TRUE then plot values in the heatmap cells. If FALSE do not plot values,
#' just plot coloured cells instead.
#'
#' @param .scientific If TRUE then use the scientific notation for numbers (e.g., "2.0e+2").
#'
#' @param .signif.digits Number of significant digits to display on plot.
#'
#' @param .text.size Size of text in the cells of heatmap.
#'
#' @param .labs A character vector of length two with names for x-axis and y-axis, respectively.
#'
#' @param .title The The text for the plot's title.
#'
#' @param .leg.title The The text for the plots's legend. Provide NULL to remove the legend's title completely.
#'
#' @param .legend If TRUE then displays a legend, otherwise removes legend from the plot.
#'
#' @param .na.value Replace NA values with this value. By default they remain NA.
#'
#' @param ... Other passed arguments.
#'
#' @seealso \link{vis}, \link{repOverlap}.
#'
#' @export vis_heatmap
vis_heatmap <- function (.data, .text = T, .scientific = FALSE, .signif.digits = 4, .text.size = 4,
                         .labs = c('Sample', 'Sample'), .title = "Overlap",
                         .leg.title = 'Overlap values', .legend = T,
                         .na.value = NA, ...) {
  if (has_class(.data, 'data.frame')) {
    names <- .data[,1]
    .data <- as.matrix(.data[,-1])
    row.names(.data) <- names
  } else if (is.null(dim(.data))) {
    .data = as.matrix(.data)
  }

  if (is.null(colnames(.data))) {
    colnames(.data) <- paste0('C', 1:ncol(.data))
  }

  if (is.null(row.names(.data))) {
    row.names(.data) <- paste0('C', 1:nrow(.data))
  }

  .data[is.na(.data)] <- .na.value

  tmp <- as.data.frame(.data)
  tmp$name <- row.names(.data)

  m <- reshape2::melt(tmp, id.var = c('name'))
  m[,1] <- factor(m[,1], levels = rev(rownames(.data)))
  m[,2] <- factor(m[,2], levels = colnames(.data))

  m$label <- format(m$value, scientific = .scientific, digits = .signif.digits)

  p <- ggplot(m, aes(x = variable, y = name, fill = value))
  p <- p + geom_tile(aes(fill = value), colour = "white")

  if (.text) { p <- p + geom_text(aes(fill = value, label = label), size = .text.size) }

  p <- p + scale_fill_distiller(palette="RdBu")

  p <- p + ggtitle(.title) +
    guides(fill = guide_colourbar(title=.leg.title)) +
    xlab(.labs[1]) + ylab(.labs[2]) + coord_fixed() +
    theme_linedraw() + theme(axis.text.x  = element_text(angle=90, vjust = .5)) +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))

  if (!.legend) { p = .rem_legend(p) }

  if (is.null(.labs)) { p = .rem_xlab(.rem_ylab(p)) }

  p
}


#' Visualisation of matrices using heatmap3-based heatmaps
#'
#' @name vis_heatmap2
#'
#' @description Visualise matrices with the functions based on the \link[heatmap3]{heatmap3}
#' package with minimum amount of arguments.
#'
#' @param .data Input matrix. Column names and row names (if presented) will be used as names for labs.
#'
#' @param .grid.size Line width of the grid between cells.
#'
#' @param .title The text for the plot's title.
#'
#' @param .labs A character vector of length two with names for x-axis and y-axis, respectively.
#'
#' @param ... Other arguments for the \link[heatmap3]{heatmap3} function.
#'
#' @seealso \link{vis}, \link{repOverlap}.
#'
#' @export vis_heatmap2
vis_heatmap2 <- function (.data, .grid.size = 1.5, .title = NULL, .labs = NULL, ...) {
  require(heatmap3)

  heatmap3::heatmap3(.data, highlightCell = expand.grid(1:nrow(.data), 1:ncol(.data), "black", .grid.size),
                     col = colorRampPalette(c("#67001f", "#d6604d", "#f7f7f7", "#4393c3", "#053061"))(1024),
                     main = .title, xlab = .labs[1], ylab = .labs[2], ...)
}


#' Visualisation of matrices using circos plots
#'
#' @name vis_circos
#'
#' @description Visualise matrices with the \link{chordDiagram} function
#' from the circlize package.
#'
#' @param .data Input matrix.
#'
#' @param .title The The text for the title of the plot.
#'
#' @param ... Other arguments passed to \link{chordDiagram} from the 'circlize' package.
#'
#' @seealso \link{vis}, \link{repOverlap}.
#'
#' @export vis_circos
vis_circos <- function (.data, .title = NULL, ...) {
  require(circlize)
  chordDiagram(.data, ...)
}


#' Visualisation of distance matrices with polar area plots
#'
#' @name vis_radar
#'
#' @description Visualise distance matrices using polar area plots from ggplot2.
#' Due to the nature of polar area plots, it is impossible to visualise all distances on the
#' same plot. So the output of the function will be a number of polar area plots, one for each
#' sample in the input matrix. Each plot visualises a distance from the specific samples
#' to other samples in the dataset.
#'
#' @param .data Input matrix.
#'
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as, age, serostatus or hla.
#'
#' @param .ncol An integer number of columns to display. Provide NA (by default)
#' if you want the function to automatically detect the optimal number of columns.
#'
#' @param .which A character vector with sample names. Samples not presented in
#' the vector will be filtered out and will not be displayed in the resulting plot.
#'
#' @param .errorbars A numeric vector of length two with quantiles for error bars
#' on sectors. Disabled if ".errorbars.off" is TRUE.
#'
#' @param .errorbars.off If TRUE then plot CI bars for distances between each group.
#' Disabled if no group passed to the ".by" argument.
#'
#' @param .expand Passed to ggplot2's \link{scale_y_continuous} argument "expand".
#'
#' @param .title The text for the title of the plot.
#'
#' @param .subtitle The text for the subtitle of the plot. Provide NULL if you want to remove it.
#'
#' @param .legend If TRUE then plots the legend, otherwise removes the legend from the plot.
#'
#' @param .leg.title The text for the legend of the plot. Pass NULL if you want to remove it.
#'
#' @param ... Do not provide anything here.
#'
#' @export vis_radar
vis_radar <- function (.data, .by = NA, .meta = NA,
                       .ncol = NA, .which = NA,
                       .errorbars = c(0.025, 0.975), .errorbars.off = F, .expand = c(.25, 0),
                       .title = NA, .subtitle = NULL,
                       .legend = F, .leg.title = NULL, ...) {
  .data = melt(.data)
  colnames(.data) = c("Source", "Target", "Value")
  .data$Value[is.na(.data$Value)] <- 0

  group_res = process_metadata_arguments(.data, .by, .meta, .defgroupby = "Source", .data.sample.col = "Source")
  group_res2 = process_metadata_arguments(.data, .by, .meta, .defgroupby = "Target", .data.sample.col = "Target")
  .data$Source.gr = group_res$group_column
  .data$Target.gr = group_res2$group_column

  .data_proc = .data %>%
    dplyr::group_by(Source.gr, Target.gr) %>%
    dplyr::summarise(Value = mean(Value, na.rm = T),
                     Value.min = quantile(Value, .errorbars[1], na.rm = T),
                     Value.max = quantile(Value, .errorbars[2], na.rm = T))
  .data_proc = .data_proc[!is.na(.data_proc$Source.gr), ]

  step = length(unique(.data_proc$Source.gr))

  if (!is.na(.which[1])) {
    .data_proc = .data_proc[.data_proc$Source.gr %in% .which, ]
  }

  if (is.na(.leg.title)) {
    if (group_res$is_grouped) {
      .leg.title = "Group"
    } else {
      .leg.title = "Sample"
    }
  }

  if (!group_res$is_grouped) {
    .errorbars.off = T
  }

  if (is.na(.ncol)) {
    .ncol = round(length(unique(.data_proc$Source.gr)) ** .5)
  }

  .data_proc$Value[.data_proc$Source.gr == .data_proc$Target.gr] = 0
  .data_proc$Value.min[.data_proc$Source.gr == .data_proc$Target.gr] = NA
  .data_proc$Value.max[.data_proc$Source.gr == .data_proc$Target.gr] = NA

  ps <- lapply(seq(1, nrow(.data_proc), step), function (l) {
    p = ggplot(.data_proc[l:(l+step-1),], aes(x = Target.gr, y = Value, fill = Target.gr)) +
      geom_bar(colour = 'black', stat = 'identity')

    if (!.errorbars.off) {
      p = p +
        geom_errorbar(aes(ymin = Value.min, ymax = Value.max), width=0.2, position=position_dodge(.9))
    }

    p = p +
      coord_polar() +
      scale_y_continuous(expand = .expand) +
      theme_linedraw() +
      .colourblind_discrete(step)

    if (!.legend) {
      p = p + guides(fill = "none")
    } else {
      p = p + guides(fill = guide_legend(title=.leg.title))
    }

    p = p +
      labs(title = .data_proc$Source.gr[l], subtitle = .subtitle)

    p
  })
  do.call(grid.arrange, c(ps, ncol = .ncol, top = .title))
}


##### Complex overlaps - top overlap & public repertoires #####


#' Visualise incremental overlaps
#'
#' @name vis.immunr_top_overlap
#'
#' @description NOT IMPLEMENTED YET
#'
#' @export vis.immunr_top_overlap
vis.immunr_top_overlap <- function (.data) {
  stop(IMMUNR_ERROR_NOT_IMPL)
}

#' Public repertoire visualisation
#'
#' @name vis.immunr_public_repertoire
#'
#' @importFrom grid rectGrob
#' @importFrom grid gpar
#'
#' @description Visualise public clonotypes.
#'
#' @param .data Public repertoire data — an output from the \link{pubRep} function.
#'
#' @param .x.rep Either indices of samples or character vector of sample names
#' for the x-axis. Must be of the same length as ".y.rep".
#'
#' @param .y.rep Either indices of samples or character vector of sample names
#' for the y-axis. Must be of the same length as ".x.rep".
#'
#' @param .title The text for the title of the plot.
#'
#' @param .ncol An integer number of columns to print in the grid of pairs of repertoires.
#'
#' @param .point.size.modif An integer value that is a modifier of the point size.
#' The larger the number, the larger the points.
#'
#' @param .cut.axes If TRUE then axes limits become shorter.
#'
#' @param .density If TRUE then displays density plot for distributions of clonotypes
#' for each sample. If FALSE then removes density plot from the visualisation.
#'
#' @param .lm If TRUE then fit a linear model and displays an R adjusted coefficient
#' that shows how similar samples are in terms of shared clonotypes.
#'
#' @param .radj.size An integer value, that defines the size of the The text
#' for the R adjusted coefficient.
#'
#' @param .return.grob If TRUE then returh the gridArrange grob instead of the plot.
#'
#' @seealso \link{pubRep}
#'
#' @export vis.immunr_public_repertoire
vis.immunr_public_repertoire <- function (.data, .x.rep = NA, .y.rep = NA,
                                          .title = NA, .ncol = 3,
                                          .point.size.modif = 1, .cut.axes = T,
                                          .density = T, .lm = T, .radj.size = 3.5, .return.grob = F) {
  .shared.rep = .data

  mat <- public_matrix(.shared.rep)

  if (is.na(.x.rep) && is.na(.y.rep)) {
    ps <- list()
    for (i in 1:ncol(mat)) {
      for (j in 1:ncol(mat)) {
        ps <- c(ps, list(vis.immunr_public_repertoire(.shared.rep, i, j, '', .point.size.modif = .point.size.modif,
                                               .cut.axes = .cut.axes, .density = .density, .lm = .lm, .radj.size = .radj.size,
                                               .return.grob = T)))
      }
    }
    grid.arrange(grobs = ps, ncol = .ncol, top = .title)
  } else if (is.na(.x.rep)) {
    ps <- lapply(1:ncol(mat), function (i) {
      vis.immunr_public_repertoire(.shared.rep, i, .y.rep, '', .point.size.modif = .point.size.modif,
                            .cut.axes = .cut.axes, .density = .density, .lm = .lm)
    })
    do.call(grid.arrange, c(ps, ncol = .ncol, top = .title))
  } else if (is.na(.y.rep)) {
    ps <- lapply(1:ncol(mat), function (j) {
      vis.immunr_public_repertoire(.shared.rep, .x.rep, j, '', .point.size.modif = .point.size.modif,
                            .cut.axes = .cut.axes, .density = .density, .lm = .lm)
    })
    do.call(grid.arrange, c(ps, ncol = .ncol, top = .title))
  } else {
    if (!is.character(.x.rep)) { .x.rep <- colnames(mat)[.x.rep] }
    if (!is.character(.y.rep)) { .y.rep <- colnames(mat)[.y.rep] }

    if (.x.rep == .y.rep) {
      return(rectGrob(gp=gpar(col="white")))
    }

    df <- data.frame(cbind(mat[, .x.rep], mat[, .y.rep]))
    df[,1] = df[,1] / sum(df[,1], na.rm = T)
    df[,2] = df[,2] / sum(df[,2], na.rm = T)
    df_full = df
    df <- df[!is.na(df[,1]) & !is.na(df[,2]), ]
    freq <- log10(sqrt(as.numeric(df[, 1]) * df[, 2])) / 2
    names(df) <- c("Xrep", "Yrep")
    names(df_full) <- c("Xrep", "Yrep")

    pnt.cols <- log(df[, 1] / df[, 2])
    suppressWarnings(pnt.cols[pnt.cols > 0] <- pnt.cols[pnt.cols > 0] / max(pnt.cols[pnt.cols > 0]))
    suppressWarnings(pnt.cols[pnt.cols < 0] <- -pnt.cols[pnt.cols < 0] / min(pnt.cols[pnt.cols < 0]))

    if (.cut.axes) {
      mat.lims <- c(min(as.matrix(df_full), na.rm = T), max(as.matrix(df_full), na.rm = T))
    } else {
      mat.lims <- c(min(as.matrix(df_full), na.rm = T), 1)
    }

    empty <- ggplot()+geom_point(aes(1,1), colour="white") +
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )

    min_df = min(floor(log10(min(df_full[,1], na.rm = T))), floor(log10(min(df_full[,2], na.rm = T))))
    max_df = max(trunc(log10(max(df_full[,1], na.rm = T))), trunc(log10(max(df_full[,2], na.rm = T))))
    breaks_values = 10**seq(min_df, 1)
    breaks_labels = format(log10(breaks_values), scientific = F)

    grey_col = "#CCCCCC"

    points = ggplot() +
      geom_point(aes(x = Xrep, y = Yrep, size = freq, fill = pnt.cols), data = df, shape=21) +
      scale_radius(range = c(.point.size.modif, .point.size.modif * 6)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      theme_linedraw() +
      scale_fill_distiller(palette="RdBu") +
      scale_x_log10(breaks = breaks_values, labels = breaks_labels, lim = mat.lims, expand = c(.015, .015)) +
      scale_y_log10(breaks = breaks_values, labels = breaks_labels, lim = mat.lims, expand = c(.015, .015)) +
      theme(legend.position="none") +
      xlab(.x.rep) + ylab(.y.rep)
    if (!is.na(.title)) {
      points = points + ggtitle(.title)
    }
    if (.lm) {
      adj.R.sq = summary(lm(Yrep ~ Xrep, df))$adj.

      points = points +
        geom_smooth(aes(x = Xrep, y = Yrep), method = "lm", data = df, fullrange = T, colour = "grey20", size = .5) +
        geom_text(aes(x = max(df_full, na.rm = T) / 4,
                      y = min(df_full, na.rm = T),
                      label = paste0("R^2(adj.) = ", as.character(round(adj.R.sq, 2)))), size = .radj.size)
      # ggtitle(paste0("R^2(adj.) = ", as.character(round(adj.R.sq, 2))))
    }

    if (.density) {
      df2 = data.frame(Clonotype = df_full[!is.na(df_full[,1]) & is.na(df_full[,2]), 1], Type = "unique", stringsAsFactors = F)
      df2 = rbind(df2, data.frame(Clonotype = df_full[!is.na(df_full[,1]) & !is.na(df_full[,2]), 1], Type = "public", stringsAsFactors = F))
      top_plot = ggplot() +
        geom_density(aes(x = Clonotype, fill = Type), colour = "grey25", data = df2, alpha = .3) +
        scale_x_log10(breaks = 10**(seq(min_df, 0)), lim = mat.lims, expand = c(.12, .015)) +
        theme_bw() + theme(axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks = element_blank(),
                           legend.position = "none") +
        scale_fill_manual(values = colorRampPalette(c(.colourblind.vector()[5], grey_col))(2))

      df2 = data.frame(Clonotype = df_full[!is.na(df_full[,2]) & is.na(df_full[,1]), 2], Type = "unique", stringsAsFactors = F)
      df2 = rbind(df2, data.frame(Clonotype = df_full[!is.na(df_full[,2]) & !is.na(df_full[,1]), 2], Type = "public", stringsAsFactors = F))
      right_plot = ggplot() +
        geom_density(aes(x = Clonotype, fill = Type), colour = "grey25", data = df2, alpha = .3) +
        scale_x_log10(breaks = 10**(seq(min_df, 0)), lim = mat.lims, expand = c(.12, .015)) +
        coord_flip() +
        theme_bw() + theme(axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks = element_blank(),
                           legend.position = "none") +
        scale_fill_manual(values = colorRampPalette(c(.colourblind.vector()[1], grey_col))(2))

      if (!.return.grob) {
        grid.arrange(top_plot, empty, points, right_plot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
      } else {
        arrangeGrob(top_plot, empty, points, right_plot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
      }
    } else {
      points
    }

  }
}


##### Gene usage & histogram plot, boxplot and treemap #####


#' Histograms and boxplots (general case / gene usage)
#'
#' @name vis.immunr_gene_usage
#'
#' @description Visualise distributions of genes using heatmaps, treemaps or other plots.
#'
#' @param .data Output from the \link{geneUsage} function.
#'
#' @param .plot String specifying the plot type:
#'
#' - "hist" for histograms using \link{vis_hist};
#'
#' - "tree" for histograms using \link{vis_treemap};
#'
#' - "heatmap" for heatmaps using \link{vis_heatmap};
#'
#' - "heatmap2" for heatmaps using \link{vis_heatmap2};
#'
#' - "circos" for circos plots using \link{vis_circos}.
#'
#' @param ... Other arguments passed to corresponding functions depending on the plot type:
#'
#' - "hist" — passes arguments to \link{vis_hist};
#'
#' - "box" — passes arguments to \link{vis_box};
#'
#' - "tree" — passes arguments to \link{vis_treemap} and \link{treemap} from the treemap package;
#'
#' - "heatmap" — passes arguments to \link{vis_heatmap};
#'
#' - "heatmap2" — passes arguments to \link{vis_heatmap2} and \link{heatmap} from the "heatmap3" package;
#'
#' - "circos" — passes arguments to \link{vis_circos} and \link{chordDiagram} from the "circlize" package.
#'
#' @examples
#' \dontrun{
#' data(twb)
#' gu = geneUsage(twb)
#'
#' # next two plots are similar
#' vis(gu, .title = "Gene usage on twins")
#' vis_hist(gu, .title = "Gene usage on twins")
#'
#' # next two plots are similar as well
#' vis(gu, "box", .title = "Gene usage on twins")
#' vis_box(gu, .title = "Gene usage on twins")
#' }
#'
#' @seealso \link{geneUsage}
#'
#' @export vis.immunr_gene_usage
vis.immunr_gene_usage <- function (.data, .plot = c("hist", "box", "tree", "heatmap", "heatmap2", "circos"), ...) {
  .plot = .plot[1]
  if (.plot == "hist") {
    vis_hist(.data, ...)
  } else if (.plot == "box") {
    vis_box(.data, .melt = T, .grouping.var = "Gene",
            .labs = c("Gene", "Usage"), .title = "Gene usage",
            .subtitle = NULL, ...) +
      theme(axis.text.x = element_text(angle=90, vjust = .5))
  } else if (.plot == "tree") {
    vis_treemap(.data, ...)
  } else if (.plot == "heatmap") {
    row.names(.data) = .data[[1]]
    .data = t(as.matrix(.data[2:ncol(.data)]))
    vis_heatmap(.data, ...)
  } else if (.plot == "heatmap2") {
    row.names(.data) = .data[[1]]
    .data = t(as.matrix(.data[2:ncol(.data)]))
    vis_heatmap2(.data, ...)
  } else if (.plot == "circos") {
    row.names(.data) = .data[[1]]
    .data = t(as.matrix(.data[2:ncol(.data)]))
    vis_circos(.data, ...)
  } else {
    stop("Error: Unknown value of the .plot parameter. Please provide one of the following: 'hist', 'box', 'tree'.")
  }
}


#' Visualisation of distributions using histograms
#'
#' @name vis_hist
#'
#' @description Visualisation of distributions using ggplot2-based histograms.
#'
#' @param .data Input matrix or data frame.
#'
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as, age, serostatus or hla.
#'
#' @param .title The text for the title of the plot.
#'
#' @param .ncol A number of columns to display. Provide NA (by default) if you want the function
#' to automatically detect the optimal number of columns.
#'
#' @param .points A logical value defining whether points will be visualised or not.
#'
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#'
#' @param .coord.flip If TRUE then swap x- and y-axes.
#'
#' @param .grid If TRUE then plot separate visualisations for each gene segment.
#'
#' @param .labs A character vector of length two with names for x-axis and y-axis, respectively.
#'
#' @param .return.grob If TRUE then returh the gridArrange grob instead of the plot.
#'
#' @param .melt If TRUE then apply \link{melt} to the ".data" before plotting.
#' In this case ".data" is supposed to be a data frame with the first character column reserved
#' for names of genes and other numeric columns reserved to counts or frequencies of genes.
#' Each numeric column should be associated with a specific repertoire sample.
#'
#' @param .legend If TRUE then plots the legend. If FALSE removes the legend from the plot.
#' If NA automatically detects the best way to display legend.
#'
#' @param ... Is not used here.
#'
#' @details
#' If data is grouped, than statistical tests for comparing means of groups will be performed, unless \code{.test = F} is supplied.
#' In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed
#' (R function \code{\link{wilcox.test}} with an argument \code{exact = F}) for testing if there is a difference in mean rank values between two groups.
#' In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) is performed (R function \code{\link{kruskal.test}}), that is equivalent to ANOVA for ranks and it tests whether samples from different groups originated from the same distribution.
#' A significant Kruskal-Wallis test indicates that at least one sample stochastically dominates one other sample.
#' Adjusted for multiple comparisons P-values are plotted on the top of groups.
#' P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method) (also known as Holm-Bonferroni correction).
#' You can execute the command \code{?p.adjust} in the R console to see more.
#'
#' @seealso \link{vis.immunr_gene_usage}, \link{geneUsage}
#'
#' @export vis_hist
vis_hist <- function (.data, .by = NA, .meta = NA, .title = "Gene usage", .ncol = NA,
                      .points = T, .test = T, .coord.flip = F,
                      .grid = T, .labs = c("Gene", NA), .return.grob = F, .melt = T, .legend = NA, ...) {
  res = .data

  if (.melt) {
    res = reshape2::melt(res)
    res = res[1:nrow(res), ]
    if (ncol(.data) == 2) {
      res[[2]] = "Data"
    }
    colnames(res) = c('Gene', 'Sample', 'Freq')
  }

  if (is.na(.labs[2])) {
    .labs[2] = "Frequency"
    if (sum(res$Freq > 1, na.rm = T) > 0) {
      .labs[2] = "Count"
    }
  }

  if (length(unique(res$Sample)) == 1) {
    p <- ggplot() +
      geom_bar(aes(x = Gene, y = Freq, fill = Freq), data = res, stat = 'identity', colour = 'black')

    expand_vec = c(.02, 0)

    p = p + theme_linedraw() +
      theme(axis.text.x = element_text(angle=90, vjust = .5)) +
      labs(x = .labs[1], y = .labs[2], title = .title, fill = .labs[2]) +
      scale_fill_distiller(palette = "RdBu")

    if (.coord.flip) {
      p = p + coord_flip() + scale_y_continuous(expand = c(.005, 0))
    } else {
      p = p + scale_y_continuous(expand = c(.02, 0))
    }

    if (!is.na(.legend)) {
      if (!.legend) {
        p = p + guides(fill = "none")
      }
    }

    p
  }
  else {
    if (.grid) {
      if (is.na(.legend)) {
        .legend = F
      }

      if (is.na(.ncol)) {
        .ncol = round(length(unique(res$Sample)) ** .5)
      }

      res <- split(res, res$Sample)
      ps = list()
      for (i in 1:length(res)) {
        ps[[i]] = vis_hist(res[[i]], .title = names(res)[i], .ncol = NA, .coord.flip = .coord.flip,
                           .grid = F, .labs = c("Gene", NA), .return.grob = F, .melt = F, .legend = .legend, ...)
      }

      if (.return.grob) {
        do.call(gridExtra::arrangeGrob,  c(ps, ncol = .ncol, top = .title))
      } else {
        do.call(gridExtra::grid.arrange, c(ps, ncol = .ncol, top = .title))
      }
    } else {
      res$Value = res$Freq

      vis_bar(.data = res, .by = .by, .meta = .meta,
              .errorbars = c(0.025, 0.975), .errorbars.off = F, .stack = F,
              .points = .points, .test = .test, .signif.label.size = 3.5, .errorbar.width = 0.45,
              .defgroupby = "Sample", .grouping.var = "Gene",
              .labs = .labs,
              .title = .title, .subtitle = NULL,
              .legend=T, .leg.title = NA) +
        theme(axis.text.x = element_text(angle=90, vjust = .5))
    }
  }
}


#' Visualisation of distributions using boxplots
#'
#' @name vis_box
#'
#' @description Visualisation of distributions using ggplot2-based boxplots.
#'
#' @param .data Input matrix or data frame.
#'
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as, age, serostatus or hla.
#'
#' @param .title The text for the title of the plot.
#'
#' @param .ncol Number of columns to display. Provide NA (by default) if you want the function
#' to automatically detect the optimal number of columns.
#'
#' @param .labs Character vector of length two with names for x-axis and y-axis, respectively.
#'
#' @param .melt If TRUE then apply \link{melt} to the ".data" before plotting.
#' In this case ".data" is supposed to be a data frame with the first character column reserved
#' for names of genes and other numeric columns reserved to counts or frequencies of genes.
#' Each numeric column should be associated with a specific repertoire sample.
#'
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#' @param .title The The text for the plot's title.
#'
#' @param .leg.title The The text for the plots's legend. Provide NULL to remove the legend's title completely.
#'
#' @param .legend If TRUE then displays a legend, otherwise removes legend from the plot.
#'
#' @seealso \link{vis.immunr_gene_usage}, \link{geneUsage}
#'
#' @export vis_box
vis_box <- function (.data, .by = NA, .meta = NA, .melt = T,
                     .points=T, .test = T, .signif.label.size = 3.5, .defgroupby = "Sample", .grouping.var = "Group",
                     .labs = c("X", "Y"), .title = "Boxplot (.title argument)",
                     .subtitle = "Subtitle (.subtitle argument)", .legend=NA, .leg.title = "Legend (.leg.title argument)") {
  if (.melt) {
    res = reshape2::melt(.data)
    res = res[1:nrow(res), ]
    if (ncol(.data) == 2) {
      res[[2]] = "Data"
    }
    colnames(res) = c(.grouping.var, 'Sample', 'Value')
  }
  .data = res

  group_res = process_metadata_arguments(.data, .by, .meta, .defgroupby)
  group_column = group_res$name
  .data$Group = group_res$group_column

  if (! (.grouping.var %in% names(.data))) {
    stop("Error: no column '", .grouping.var, "' in the input data frame.")
  }

  .data$Grouping.var = .data[[.grouping.var]]

  .data_proc = .data %>%
    dplyr::group_by(Grouping.var, Group)

  if (is.na(.labs[1])) {
    if (group_res$is_grouped) {
      .labs[1] = "Group"
    } else {
      .labs[1] = "Sample"
    }
  }

  if (is.na(.leg.title)) {
    if (group_res$is_grouped) {
      .leg.title = "Group"
    } else {
      .leg.title = "Sample"
    }
  }

  if (group_res$is_grouped) {
    # p = ggplot(aes(x = Grouping.var, y = Value, fill = Group, color=Group, group = Group), data = .data) +
    #     geom_boxplot(position = "dodge", col = "black")
    p = ggplot(aes(x = Grouping.var, y = Value, fill = Group), data = .data) +
      geom_boxplot(position = "dodge", col = "black")

    if (is.na(.legend)) {
      if (.grouping.var == "Group") {
        p = p + guides(fill=F)
      }
    } else if (.legend == F) {
      p = p + guides(fill=F)
    }

    if (.points) {
      p = p +
        geom_point(color = "black", position=position_jitterdodge(0.05), size = 1)
    }

    if (.test) {
      if (.grouping.var == "Group") {
        comparisons = list()
        for (i in 1:(nrow(.data_proc)-1)) {
          for (j in (i+1):nrow(.data_proc)) {
            # if ((.data$Grouping.var[i] == .data$Grouping.var[i]) && (.data$Group[i] != .data$Group[j])) {
            comparisons = c(comparisons, list(c(i, j)))
            # }
          }
        }

        p_df = compare_means(Value ~ Group, .data, comparisons=comparisons, p.adjust.method = "holm")

        y_max = max(.data$Value)
        p.value.y.coord <- rep(y_max, nrow(p_df))
        step.increase <- (1:nrow(p_df))*(y_max/10)
        p.value.y.coord <- p.value.y.coord + step.increase

        p_df <- p_df %>%
          mutate(y.coord =  p.value.y.coord,
                 p.adj = format.pval(p.adj, digits = 1))

        p = p + ggsignif::geom_signif(data = p_df,
                                      aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
                                      manual= TRUE, tip_length = 0.03, size = .5, inherit.aes = F)
      } else {
        # Seems fine...
        # p_df = compare_means(Value ~ Group, group.by = "Grouping.var", method = "kruskal.test", .data, p.adjust.method = "holm")
        # print(p_df)

        p = p +
          ggpubr::stat_compare_means(aes(label = ..p.adj..), bracket.size = .5, size=.signif.label.size,
                                     label.y = max(.data$Value, na.rm = T) * 1.07)
      }
    }

    p = p + .tweak_fill(length(unique(.data$Group)))

  } else {
    p = ggplot() +
      geom_boxplot(aes(x = Grouping.var, y = Value, fill = Sample), position = "dodge", data = .data, col = "black")

    if (is.na(.legend)) {
      p = p + guides(fill=F)
    }

    p = p + .colourblind_discrete(length(unique(.data$Sample)))
  }

  p + theme_linedraw() + labs(x = .labs[1], y = .labs[2], title = .title, subtitle = .subtitle, fill = group_column)
}

#' Visualisation of data frames and matrices using treemaps
#'
#' @name vis_treemap
#'
#' @importFrom treemap treemap
#'
#' @param .data Input data frame or matrix.
#'
#' @param ... Other arguments passed to \link{treemap}.
#'
#' @export vis_treemap
vis_treemap <- function (.data, ...) {
  melted = reshape2::melt(.data)
  colnames(melted) = c("name", "group", "value")
  treemap::treemap(melted, index=c("group", "name"), vSize="value", ...)
}


##### Clustering #####


#' Visualisation of clustering results
#'
#' @name vis.immunr_hclust
#'
#' @aliases vis.immunr_hclust vis.immunr_kmeans immunr_dbscan
#'
#' @usage
#' vis.immunr_hclust(.data, .rect = F, .plot = c("clust", "best"), .return.grob = F)
#'
#' vis.immunr_kmeans(.data, .point = T, .text = T, .ellipse = T,
#' .point.size = 2, .text.size = 10, .plot = c("clust", "best"), .return.grob = F)
#'
#' vis.immunr_dbscan(.data, .point = T, .text = T, .ellipse = T,
#' .point.size = 2, .text.size = 10, .plot = c("clust", "best"))
#'
#' @param .data Clustering results from \link{repOverlapAnalysis} or \link{geneUsageAnalysis}.
#'
#' @param .rect Passed to \link{fviz_dend} - whether to add a rectangle around groups.
#'
#' @param .point If TRUE then plot sample points. Passed to \link{fviz_cluster}.
#'
#' @param .text If TRUE then plot text labels. Passed to \link{fviz_cluster}.
#'
#' @param .ellipse If TRUE then plot ellipses around all samples. Passed to "ellipse" from \link{fviz_cluster}.
#'
#' @param .point.size Size of points, passed to "pointsize" from \link{fviz_cluster}.
#'
#' @param .text.size Size of text labels, passed to labelsize from \link{fviz_cluster}.
#'
#' @param .plot A character vector of length one or two specifying which plots to visualise.
#' If "clust" than plot only the clustering. If "best" than plot the number of optimal clusters.
#' If both than plot both.
#'
#' @param .return.grob If TRUE then return grob instead of plot.
#'
#' @seealso \link{vis}, \link{repOverlapAnalysis}, \link{geneUsageAnalysis}
#'
#' @export vis.immunr_hclust vis.immunr_kmeans vis.immunr_dbscan
vis.immunr_hclust <- function (.data, .rect = F, .plot = c("clust", "best"), .return.grob = F) {
  p1 = NULL
  if ("clust" %in% .plot) {
    p1 = fviz_dend(.data[[1]], main = "Hierarchical clustering", rect = .rect)
  }

  p2 = NULL
  if ("best" %in% .plot) {
    p2 = .data[[2]]
  }

  if (is.null(p1) | is.null(p2)) {
    if (is.null(p1)) { p2 } else { p1 }
  } else {
    if (.return.grob) {
      gridExtra::arrangeGrob(p1, p2, ncol = 2)
    } else {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    }
  }
}

vis.immunr_kmeans <- function (.data, .point = T, .text = T, .ellipse = T,
                               .point.size = 2, .text.size = 10, .plot = c("clust", "best"),
                               .return.grob = F) {
  p1 = NULL
  if ("clust" %in% .plot) {
    p1 = fviz_cluster(.data[[1]], data = .data[[3]], main = "K-means clustering", geom = c("point", "text")[c(.point, .text)],
                      show.legend.text = F, show.clust.cent = F, repel=T, ellipse = .ellipse, shape = 16,
                      pointsize = .point.size, labelsize = .text.size, label.rectangle = T) +
      guides(fill = "none", shape = "none") +
      theme_linedraw()
  }

  p2 = NULL
  if ("best" %in% .plot) {
    p2 = .data[[2]]
  }

  if (is.null(p1) | is.null(p2)) {
    if (is.null(p1)) { p2 } else { p1 }
    # p = ifelse(is.null(p1), p2, p1) # doesnt work for some reason
  } else {
    if (.return.grob) {
      gridExtra::arrangeGrob(p1, p2, ncol = 2)
    } else {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    }
  }
}

vis.immunr_dbscan <- function (.data, .point = T, .text = T, .ellipse = T,
                               .point.size = 2, .text.size = 10, .plot = c("clust", "best")) {
  fviz_cluster(.data[[1]], data = .data[[2]], main = "DBSCAN clustering", geom = c("point", "text")[c(.point, .text)],
               show.legend.text = F, show.clust.cent = F, repel=T, ellipse = .ellipse, shape = 16,
               pointsize = .point.size, labelsize = .text.size, label.rectangle = T) +
    guides(fill = "none", shape = "none") +
    theme_linedraw()
}


##### Dimension reduction #####


#' PCA / MDS / tSNE visualisation (mainly overlap / gene usage)
#'
#' @name vis.immunr_mds
#'
#' @aliases vis.immunr_mds vis.immunr_pca vis.immunr_tsne
#'
#' @param .data Output from \code{immunr_pca} or \code{immunr_mds}.
#'
#' @export vis.immunr_mds vis.immunr_pca vis.immunr_tsne
vis.immunr_mds <- function (.data, .by = NA, .meta = NA,
                            .point = T, .text = T, .ellipse = T,
                            .point.size = 2, .text.size = 4) {
  if (!.point & !.text) {
    stop("Error: Please provide at least one of the arguments: .point and .text")
  }
  group_res = process_metadata_arguments(data.frame(Sample = row.names(.data$x), stringsAsFactors = F), .by, .meta)
  group_column = group_res$name
  if (!group_res$is_grouped) {
    .ellipse = F
  }

  fviz_pca_ind(.data, habillage = group_res$group_column, geom = c("point", "text")[c(.point, .text)],
               repel=T, addEllipses = .ellipse, mean.point = F, pointshape = 16,
               pointsize = .point.size, labelsize = .text.size, label.rectangle = T, show.legend.text = F) +
    theme_linedraw() +
    guides(fill = "none", shape = "none") +
    labs(x = "Dim1", y = "Dim2", title = "Multidimensional scaling", col = group_column)
}

vis.immunr_pca <- function (.data, .by = NA, .meta = NA,
                            .point = T, .text = T, .ellipse = T,
                            .point.size = 2, .text.size = 4) {
  if (!.point & !.text) {
    stop("Error: Please provide at least one of the arguments: .point and .text")
  }
  group_res = process_metadata_arguments(data.frame(Sample = row.names(.data$x), stringsAsFactors = F), .by, .meta)
  group_column = group_res$name
  if (!group_res$is_grouped) {
    .ellipse = F
  }

  fviz_pca_ind(.data, habillage = group_res$group_column, geom = c("point", "text")[c(.point, .text)],
               repel=T, addEllipses = .ellipse, mean.point = F, pointshape = 16,
               pointsize = .point.size, labelsize = .text.size, label.rectangle = T, show.legend.text = F) +
    theme_linedraw() +
    guides(fill = "none", shape = "none") +
    labs(title = "Principal component analysis", col = group_column)
}

vis.immunr_tsne <- function (.data, .by = NA, .meta = NA,
                             .point = T, .text = T, .ellipse = T,
                             .point.size = 2, .text.size = 4) {
  .data = data.frame(.data)
  colnames(.data) = c("Dim1", "Dim2")
  .data$Sample = row.names(.data)
  group_res = process_metadata_arguments(.data, .by, .meta)
  group_column = group_res$name
  .data$Group = group_res$group_column

  if (!group_res$is_grouped) {
    .ellipse = F
  }

  ggscatter(data = .data, x = "Dim1", y = "Dim2", color="Group", ellipse = .ellipse, size = .point.size,
            point = .point, label = ifelse(.text, "Sample", NULL), repel = T, label.rectangle = T, show.legend.text = F) +
    theme_linedraw() +
    guides(fill = "none", shape = "none") +
    labs(title = "t-Distributed Stochastic Neighbor Embedding", col = group_column)
}

#
# ToDo: ideally, geneUsageAnalysis returns raw objects, and vis.immunr_mds transforms it into the fviz-compatible objects
#


##### Clonality analysis #####


vis_bar_stacked <- function (.data, .by = NA, .meta = NA,
                             .errorbars=c(0.025, 0.975), .errorbars.off = F, .stack = NA,
                             .grouping.var = NA,
                             .labs = c(NA, "Y"),
                             .title = "Barplot (.title argument)", .subtitle = "Subtitle (.subtitle argument)",
                             .legend=NA, .leg.title = NA) {

  group_res = process_metadata_arguments(.data, .by, .meta)
  .data$Group = group_res$group_column

  if (! (.grouping.var %in% names(.data))) {
    stop("Error: no column '", .grouping.var, "' in the input data frame.")
  }
  .data$Grouping.var = .data[[.grouping.var]]

  position_bar = "dodge"
  if (is.na(.stack)) {
    if (group_res$is_grouped) {
      if (.errorbars.off) {
        position_bar = "stack"
      } else {
        position_bar = "dodge"
      }
    } else {
      position_bar = "stack"
    }
  } else if (.stack) {
    if (group_res$is_grouped) {
      if (.errorbars.off) {
        position_bar = "stack"
      } else {
        warning("Warning: Provide .errorbars.off = T in order to display stacked barplots.")
        position_bar = "dodge"
      }
    } else {
      position_bar = "stack"
    }
  }

  if (is.na(.labs[1])) {
    if (group_res$is_grouped) {
      .labs[1] = "Group"
    } else {
      .labs[1] = "Sample"
    }
  }

  if (group_res$is_grouped) {
    perc = .data %>%
      dplyr::group_by(Grouping.var, Group) %>%
      dplyr::summarise(Value.mean = mean(Value, na.rm = T),
                       Value.min = quantile(Value, .errorbars[1], na.rm = T),
                       Value.max = quantile(Value, .errorbars[2], na.rm = T))

    p <- ggplot() +
      geom_bar(aes(x = Group, y = Value.mean, fill = Grouping.var), data = perc, colour = 'black', stat = 'identity', position=position_bar)

    if (!.errorbars.off) {
      p = p + geom_errorbar(aes(x = Group, fill = Grouping.var, ymin = Value.min, ymax = Value.max),
                            data = perc, colour = 'black', width=.2, position=position_dodge(.9))
    }
  } else {
    p <- ggplot() +
      geom_bar(aes(x = Sample, y = Value, fill = Grouping.var), data = .data, colour = 'black', stat = 'identity', position = position_bar)
  }

  p + theme_linedraw() +
    labs(x = .labs[1], y = .labs[2], title = .title, subtitle = .subtitle, fill = .leg.title) +
    .colourblind_discrete(length(unique(.data$Grouping.var)))
}


#' Visualise clonality
#'
#' @name vis.immunr_clonal_prop
#'
#' @aliases vis.immunr_clonal_prop vis.immunr_homeo vis.immunr_top_prop vis.immunr_tail_prop
#'
#' @description An utility function to visualise the output from \code{\link{repClonality}}.
#'
#' @importFrom reshape2 melt
#'
#' @param .data Output from \code{\link{repClonality}}.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as, age, serostatus or hla.
#'
#' @param .errorbars A numeric vector of length two with quantiles for error bars
#' on sectors. Disabled if ".errorbars.off" is TRUE.
#'
#' @param .errorbars.off If TRUE then plot CI bars for distances between each group.
#' Disabled if no group passed to the ".by" argument.
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#'
#' @details
#' If data is grouped, than statistical tests for comparing means of groups will be performed, unless \code{.test = F} is supplied.
#' In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed
#' (R function \code{\link{wilcox.test}} with an argument \code{exact = F}) for testing if there is a difference in mean rank values between two groups.
#' In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) is performed (R function \code{\link{kruskal.test}}), that is equivalent to ANOVA for ranks and it tests whether samples from different groups originated from the same distribution.
#' A significant Kruskal-Wallis test indicates that at least one sample stochastically dominates one other sample.
#' Adjusted for multiple comparisons P-values are plotted on the top of groups.
#' P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method) (also known as Holm-Bonferroni correction).
#' You can execute the command \code{?p.adjust} in the R console to see more.
#'
#' @seealso \link{repClonality} \link{vis}
#'
#' @export vis.immunr_clonal_prop vis.immunr_homeo vis.immunr_top_prop vis.immunr_tail_prop
vis.immunr_clonal_prop <- function (.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = F, .points=T, .test = T, .signif.label.size=3.5) {
  perc_value = round(.data[1, 2][1])
  .data = data.frame(Sample = row.names(.data), Value = .data[,1])

  p = vis_bar(.data = .data, .by = .by, .meta = .meta,
          .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = F,
          .points = .points, .test = .test, .signif.label.size = .signif.label.size,
          .defgroupby = "Sample", .grouping.var = "Group",
          .labs = c(NA, "Clonal proportion"),
          .title = "Clonal proportion", .subtitle = paste0("Number of clonotypes occupying the ", perc_value, "% of repertoires"),
          .legend=NA, .leg.title = NA)

  if (is.na(.by)) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }
  p
}


vis.immunr_homeo <- function (.data, .by = NA, .meta = NA, .errorbars=c(0.025, 0.975), .errorbars.off = F, .stack = NA, .test = T, .points=T) {
  melted <- reshape2::melt(.data)
  colnames(melted) <- c('Sample', 'Clone.group', 'Value')

  if (is.na(.by[1]) && is.na(.stack)) {
    .stack = T
  }

  p = vis_bar(melted, .by = .by, .meta = .meta,
              .errorbars=.errorbars, .errorbars.off = .errorbars.off, .stack = .stack,
              .grouping.var = "Clone.group", .points = .points,
              .labs = c(NA, "Relative abundance, perc"),
              .title = "Relative abundance", .subtitle = "Summary proportion of clonotypes with specific frequencies",
              .legend=NA, .leg.title = "Clonotype group", .test=.test) +
    scale_y_continuous(labels=scales::percent)

  if (is.na(.by)[1]) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  } else {
    p = p + theme(axis.text.x  = element_text(angle=30, vjust = 1, hjust=1))
  }
  p
}

vis.immunr_top_prop <- function (.data, .by = NA, .meta = NA, .errorbars=c(0.025, 0.975), .errorbars.off = F, .stack = NA, .points=T, .test=T) {
  tmp = .data
  if (is.null(dim(tmp))) {
    tmp = t(as.matrix(tmp))
    .data = t(as.matrix(.data))
  }
  for (i in 2:ncol(.data)) {
    tmp[,i] = .data[,i] - .data[,i-1]
  }
  res = tmp
  .head = as.numeric(colnames(.data))
  colnames(res) <- paste0('[', c(1, .head[-length(.head)] + 1), ':', .head, ')')
  res <- as.data.frame(res)
  res$Sample <- row.names(res)
  res <- reshape2::melt(res)
  colnames(res) = c("Sample", "Clone.index", "Value")

  if (is.na(.by[1]) && is.na(.stack)) {
    .stack = T
  }

  p = vis_bar(res, .by = .by, .meta = .meta, .points = .points, .test = .test,
                  .errorbars=.errorbars, .errorbars.off = .errorbars.off, .stack = .stack,
                  .grouping.var = "Clone.index",
                  .labs = c(NA, "Occupied repertoire space"),
                  .title = "Top clonal proportion", .subtitle = "Summary proportion of clonotypes with specific indices",
                  .legend=NA, .leg.title = "Clonotype indices") +
    scale_y_continuous(labels=scales::percent)

  if (is.na(.by)) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }
  p
}

vis.immunr_tail_prop <- function (.data, .by = NA, .meta = NA, .errorbars=c(0.025, 0.975), .errorbars.off = F, .stack = NA, .points=T, .test=T) {
  tmp = .data
  if (is.null(dim(tmp))) {
    tmp = t(as.matrix(tmp))
    .data = t(as.matrix(.data))
  }

  colnames_new = c()
  prev_value = "1"
  start_index = 1
  if (colnames(.data)[1] == "1") {
    colnames_new = "1"
    prev_value = "2"
    start_index = 2
  }
  for (i in start_index:ncol(.data)) {
    colnames_new = c(colnames_new, paste0(prev_value, " - ", colnames(.data)[i]))
    prev_value = as.character(as.integer(colnames(.data)[i]) + 1)
  }
  colnames(tmp) = colnames_new

  for (i in 2:ncol(.data)) {
    tmp[,i] = .data[,i] - .data[,i-1]
  }
  res = tmp

  res <- as.data.frame(res)
  res$Sample <- row.names(res)
  res <- reshape2::melt(res)
  colnames(res) = c("Sample", "Counts", "Value")

  if (is.na(.by[1]) && is.na(.stack)) {
    .stack = T
  }

  p = vis_bar(res, .by = .by, .meta = .meta,
                  .errorbars=.errorbars, .errorbars.off = .errorbars.off, .stack = .stack,
                  .grouping.var = "Counts", .points = .points, .test = .test,
                  .labs = c(NA, "Occupied repertoire space"),
                  .title = "Tail clonal proportion", .subtitle = "Summary proportion of clonotypes with specific counts",
                  .legend=NA, .leg.title = "Clonotype counts") +
    scale_y_continuous(labels=scales::percent)

  if (is.na(.by)) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }
  p
}


##### Diversity estimation & dodged bar plots #####


#'
#' @name vis_bar
#'
#' @param .stack If TRUE and .errorbars.off is TRUE then plot stacked bar plots for each Group or Sample
#'
#' @export vis_bar
vis_bar <- function (.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = F, .stack = F,
                     .points=T, .test = T, .signif.label.size = 3.5, .errorbar.width = 0.2, .defgroupby = "Sample", .grouping.var = "Group",
                     .labs = c("X", "Y"), .title = "Barplot (.title argument)",
                     .subtitle = "Subtitle (.subtitle argument)", .legend=NA, .leg.title = "Legend (.leg.title argument)") {
  group_res = process_metadata_arguments(.data, .by, .meta, .defgroupby)
  group_column = group_res$name
  .data$Group = group_res$group_column

  if (! (.grouping.var %in% names(.data))) {
    stop("Error: no column '", .grouping.var, "' in the input data frame.")
  }

  .data$Grouping.var = .data[[.grouping.var]]

  .data_proc = .data %>%
    dplyr::group_by(Grouping.var, Group) %>%
    dplyr::summarise(Value.mean = mean(Value, na.rm = T),
                     Value.min = quantile(Value, .errorbars[1], na.rm = T),
                     Value.max = quantile(Value, .errorbars[2], na.rm = T))

  if (is.na(.labs[1])) {
    if (group_res$is_grouped) {
      .labs[1] = "Group"
    } else {
      .labs[1] = "Sample"
    }
  }

  if (is.na(.leg.title)) {
    if (group_res$is_grouped) {
      .leg.title = "Group"
    } else {
      .leg.title = "Sample"
    }
  }

  position_name = "dodge"
  if (group_res$is_grouped) {
    if (is.na(.stack)) {
      .stack = F
    }
  } else {
    if (is.na(.stack)) {
      .stack = F
    } else if (.stack) {
      .errorbars.off = T
      .points = F
      .test = F
    }
  }

  if (group_res$is_grouped) {
    if (.stack) {
      p = ggplot() +
        geom_bar(aes(x = Group, y = Value.mean, fill = Grouping.var),
                 stat = "identity", position = "stack", data = .data_proc, col = "black")
    } else {
      p = ggplot(aes(x = Grouping.var, y = Value, fill = Group, color=Group, group = Group), data = .data) +
        geom_bar(aes(x = Grouping.var, y = Value.mean),
                 stat = "identity", position = position_name, data = .data_proc, col = "black")
    }

    if (is.na(.legend)) {
      if (.grouping.var == "Group") {
        p = p + guides(fill=F)
      }
    } else if (.legend == F) {
      p = p + guides(fill=F)
    }

    if (!.errorbars.off) {
       p = p +
         geom_errorbar(aes(x = Grouping.var, y = Value.mean, ymin = Value.min, ymax = Value.max, color = Group),
                       color = "black", data = .data_proc, width=.errorbar.width, position=position_dodge(.9))
    }

    if (.points) {
      p = p +
        geom_point(color = "black", position=position_jitterdodge(0.05), size = 1)
    }

    if (.test) {
      if (.grouping.var == "Group") {
        comparisons = list()
        for (i in 1:(nrow(.data_proc)-1)) {
          for (j in (i+1):nrow(.data_proc)) {
            # if ((.data$Grouping.var[i] == .data$Grouping.var[i]) && (.data$Group[i] != .data$Group[j])) {
            comparisons = c(comparisons, list(c(i, j)))
            # }
          }
        }

        p_df = compare_means(Value ~ Group, .data, comparisons=comparisons, p.adjust.method = "holm")

        y_max = max(.data$Value)
        p.value.y.coord <- rep(y_max, nrow(p_df))
        step.increase <- (1:nrow(p_df))*(y_max/10)
        p.value.y.coord <- p.value.y.coord + step.increase

        p_df <- p_df %>%
          mutate(y.coord =  p.value.y.coord,
                 p.adj = format.pval(p.adj, digits = 1))

        p = p + ggsignif::geom_signif(data = p_df,
                                      aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
                                      manual= TRUE, tip_length = 0.03, size = .5, inherit.aes = F)
      } else {
        # Seems fine...
        # p_df = compare_means(Value ~ Group, group.by = "Grouping.var", method = "kruskal.test", .data, p.adjust.method = "holm")
        # print(p_df)

        p = p +
          ggpubr::stat_compare_means(aes(label = ..p.adj..), bracket.size = .5, size=.signif.label.size,
                                     label.y = max(.data$Value, na.rm = T) * 1.07)
      }
    }

    p = p + .tweak_fill(length(unique(.data$Group)))

  } else {
    if (.stack) {
      p <- ggplot() +
        geom_bar(aes(x = Sample, y = Value, fill = Grouping.var), data = .data, colour = 'black', stat = 'identity', position = "stack")

      if (!is.na(.legend)) {
        if (!.legend) {
          p = p + guides(fill=F)
        }
      }

      p =  p + .colourblind_discrete(length(unique(.data$Grouping.var)))
    } else {
      p = ggplot() +
        geom_bar(aes(x = Grouping.var, y = Value, fill = Sample), position = position_name, stat = "identity", data = .data, col = "black")

      if (is.na(.legend)) {
        p = p + guides(fill=F)
      }

      p = p + .colourblind_discrete(length(unique(.data$Sample)))
    }
  }


  p + theme_linedraw() + labs(x = .labs[1], y = .labs[2], title = .title, subtitle = .subtitle, fill = group_column)
}


#' Visualise diversity.
#'
#' @name vis.immunr_chao1
#'
#' @aliases vis.immunr_chao1 vis.immunr_dxx vis.immunr_rarefaction vis.immunr_div vis.immunr_ginisimp vis.immunr_invsimp vis.immunr_hill
#' @description An utility function to visualise the output from \code{\link{repDiversity}}.
#'
#' @importFrom reshape2 melt
#'
#' @param .data Output from \code{\link{repDiversity}}.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as, age, serostatus or hla.
#'
#' @param .errorbars A numeric vector of length two with quantiles for error bars
#' on sectors. Disabled if ".errorbars.off" is TRUE.
#'
#' @param .errorbars.off If TRUE then plot CI bars for distances between each group.
#' Disabled if no group passed to the ".by" argument.
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#'
#' @details
#' If data is grouped, than statistical tests for comparing means of groups will be performed, unless \code{.test = F} is supplied.
#' In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed
#' (R function \code{\link{wilcox.test}} with an argument \code{exact = F}) for testing if there is a difference in mean rank values between two groups.
#' In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) is performed (R function \code{\link{kruskal.test}}), that is equivalent to ANOVA for ranks and it tests whether samples from different groups originated from the same distribution.
#' A significant Kruskal-Wallis test indicates that at least one sample stochastically dominates one other sample.
#' Adjusted for multiple comparisons P-values are plotted on the top of groups.
#' P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method) (also known as Holm-Bonferroni correction).
#' You can execute the command \code{?p.adjust} in the R console to see more.
#'
#' @seealso \link{repDiversity} \link{vis}
#' @export vis.immunr_chao1 vis.immunr_dxx vis.immunr_rarefaction vis.immunr_div vis.immunr_ginisimp vis.immunr_invsimp vis.immunr_hill
vis.immunr_chao1 <- function (.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = F, .points=T, .test = T, .signif.label.size=3.5) {
  .data = data.frame(Sample = row.names(.data), Value = .data[,1])

  p = vis_bar(.data = .data, .by = .by, .meta = .meta,
          .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = F,
          .points = .points, .test = .test, .signif.label.size = .signif.label.size,
          .defgroupby = "Sample", .grouping.var = "Group",
          .labs = c(NA, "Chao1"),
          .title = "Chao1", .subtitle = "Sample diversity estimation using Chao1",
          .legend=NA, .leg.title = NA)

  if (is.na(.by)[1]) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }
  p
}

vis.immunr_hill <- function (.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = F,
                             .points=T, .add.points = T, .point.size=1.5, .add.point.size = 1, .line.size = 0.75, .leg.title = NA) {
  group_res = process_metadata_arguments(.data, .by, .meta)
  .data$Group = group_res$group_column

  if (is.na(.leg.title)) {
    if (group_res$is_grouped) {
      .leg.title = "Group"
    } else {
      .leg.title = "Sample"
    }
  }

  grouped_data = .data %>%
    dplyr::group_by(Group, Q) %>%
    dplyr::summarise(Value.mean = mean(Value, na.rm = T),
                     Value.min = quantile(Value, .errorbars[1], na.rm = T),
                     Value.max = quantile(Value, .errorbars[2], na.rm = T))

  position_value = "identity"
  if (group_res$is_grouped) {
    position_value = position_jitter(width=0.1)
  }

  p = ggplot() +
    geom_line(aes(x = Q, y = Value.mean, col = Group), size = .line.size, data = grouped_data, stat = "identity")

  if (.points) {
    p = p + geom_point(aes(x = Q, y = Value, col = Group), data = .data,
                       position=position_value, alpha = .5, size = .point.size)
  }

  if (group_res$is_grouped) {
    if (.add.points) {
      p = p +
        geom_point(aes(x = Q, y = Value.mean, col = Group), data = grouped_data, size = .add.point.size)
    }

    if (!.errorbars.off) {
      p = p +
        geom_errorbar(aes(x = Q, ymin = Value.min, ymax = Value.max, group = Group, col = Group),
                      data = grouped_data, width=.2, size = .line.size, position=position_dodge(0.1))
    }
  }

  p = p +
    theme_linedraw() +
    labs(x = "Q values", y = "Diversity estimation", title = "Hill numbers", subtitle = "Sample diversity estimation using Hill numbers", color = .leg.title)

  p
}

vis.immunr_div <- function (.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = F,
                            .points=T, .test = T, .signif.label.size=3.5,
                            .legend = NA) {
  p = vis_bar(.data = .data, .by = .by, .meta = .meta,
          .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = F,
          .points = .points, .test = .test, .signif.label.size = .signif.label.size,
          .defgroupby = "Sample", .grouping.var = "Group",
          .labs = c(NA, "Effective number of clonotypes"),
          .title = "True diversity", .subtitle = "Sample diversity estimation using the true diversity index",
          .legend=.legend, .leg.title = NA)

  if (is.na(.by)) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }
  p
}

vis.immunr_ginisimp <- function (.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = F, .points=T, .test = T, .signif.label.size=3.5) {
  p = vis_bar(.data = .data, .by = .by, .meta = .meta,
          .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = F,
          .points = .points, .test = .test, .signif.label.size = .signif.label.size,
          .defgroupby = "Sample", .grouping.var = "Group",
          .labs = c(NA, "Gini-Simpson index"),
          .title = "Gini-Simpson index", .subtitle = "Sample diversity estimation using the Gini-Simpson index",
          .legend=NA, .leg.title = NA)

  if (is.na(.by)) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }
  p
}

vis.immunr_invsimp <- function (.data, .by = NA, .meta = NA,
                                .errorbars = c(0.025, 0.975), .errorbars.off = F,
                                .points=T, .test = T, .signif.label.size=3.5,
                                .legend = NA) {
  p = vis_bar(.data = .data, .by = .by, .meta = .meta,
          .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = F,
          .points = .points, .test = .test, .signif.label.size = .signif.label.size,
          .defgroupby = "Sample", .grouping.var = "Group",
          .labs = c(NA, "Inverse Simpson index"),
          .title = "Inverse Simpson index", .subtitle = "Sample diversity estimation using the inverse Simpson index",
          .legend=.legend, .leg.title = NA)

  if (is.na(.by)) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }
  p
}

vis.immunr_dxx <- function (.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = F,
                            .points=T, .test = T, .signif.label.size=3.5, .legend = NA) {
  perc_value = round(.data[1, 2][1])
  .data = data.frame(Sample = row.names(.data), Value = .data[,1])

  p = vis_bar(.data = .data, .by = .by, .meta = .meta,
          .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = F,
          .points = .points, .test = .test, .signif.label.size = .signif.label.size,
          .defgroupby = "Sample", .grouping.var = "Group",
          .labs = c(NA, paste0("D", perc_value)),
          .title = paste0("D", perc_value, " diversity index"), .subtitle = paste0("Number of clonotypes occupying the ", perc_value, "% of repertoires"),
          .legend=.legend, .leg.title = NA)

  if (is.na(.by)) {
    p = p + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }
  p
}

vis.immunr_rarefaction <- function (.data, .by = NA, .meta = NA, .mean = T, .errors = T, .log = F, .labels = T) {
  .muc.res = .data

  group_res = process_metadata_arguments(.data, .by, .meta)
  data_groups = group_res$groups
  group_column = group_res$name
  names(data_groups) = group_res$group_names
  is_grouped = group_res$is_grouped
  .muc.res$Group = group_res$group_column

  .muc.res$Type <- factor(.muc.res$Type, levels = c('interpolation', 'extrapolation'), ordered = T)

  p <- ggplot() + xlab('Sample size') + ylab('Estimated diversity') + ggtitle("Rarefaction analysis") +
    theme_linedraw()

  if (.mean && is_grouped) {
    data_proc = .muc.res %>%
      group_by(Group, Size) %>%
      summarise(MinVal = quantile(Mean, .025, na.rm = T),
                MaxVal = quantile(Mean, .975, na.rm = T),
                MeanVal = mean(Mean, na.rm = T))

    p = p + geom_line(aes(x = Size, y = MeanVal, color=Group), data = data_proc)

    if (.errors) {
      p = p + geom_ribbon(aes(x = Size, group = Group, ymin=MinVal, ymax=MaxVal, fill=Group), data = data_proc, alpha=0.15)
    }
  } else {
    colnames(.muc.res)[2] = "Q1"
    colnames(.muc.res)[4] = "Q2"
    if (.errors) {
      p = p + geom_ribbon(aes(x = Size, group = Sample, ymin=Q1, ymax=Q2), col = "grey80", data = .muc.res, alpha=0.2)
    }

    p = p + geom_line(aes(x = Size, y = Mean, colour = Group, linegroup = Sample, linetype = Type), size=1, data = .muc.res)

    if (.labels) {
      tmp <- .muc.res[tapply(1:nrow(.muc.res), .muc.res$Sample, tail, 1), ]
      tmp = tmp[order(tmp$Group),]
      p <- p + ggrepel::geom_label_repel(aes(x = Size, y = Mean, label = Sample, fill = Group), data = tmp)
    }
  }

  if (.log) {
    p <- p + scale_x_log10()
  }

  p + .tweak_col(length(unique(.muc.res$Group))) + .tweak_fill(length(unique(.muc.res$Group)))
}


##### Exploratory analysis #####


#' Visualise exploratory analysis.
#'
#' @name vis.immunr_exp_vol
#'
#' @aliases vis.immunr_exp_vol vis.immunr_exp_count vis.immunr_exp_len
#' @description An utility function to visualise the output from \code{\link{repExplore}}.
#'
#' @importFrom reshape2 melt
#'
#' @param .data Output from \code{\link{repExplore}}.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as, age, serostatus or hla.
#'
#' @param .errorbars A numeric vector of length two with quantiles for error bars
#' on sectors. Disabled if ".errorbars.off" is TRUE.
#'
#' @param .errorbars.off If TRUE then plot CI bars for distances between each group.
#' Disabled if no group passed to the ".by" argument.
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#'
#' @details
#' If data is grouped, than statistical tests for comparing means of groups will be performed, unless \code{.test = F} is supplied.
#' In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed
#' (R function \code{\link{wilcox.test}} with an argument \code{exact = F}) for testing if there is a difference in mean rank values between two groups.
#' In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) is performed (R function \code{\link{kruskal.test}}), that is equivalent to ANOVA for ranks and it tests whether samples from different groups originated from the same distribution.
#' A significant Kruskal-Wallis test indicates that at least one sample stochastically dominates one other sample.
#' Adjusted for multiple comparisons P-values are plotted on the top of groups.
#' P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method) (also known as Holm-Bonferroni correction).
#' You can execute the command \code{?p.adjust} in the R console to see more.
#'
#' @seealso \link{repExplore} \link{vis}
#' @export vis.immunr_exp_vol vis.immunr_exp_count vis.immunr_exp_len
vis.immunr_exp_vol <- function (.data, .by = NA, .meta = NA,
                                .errorbars = c(0.025, 0.975), .errorbars.off = F,
                                .points=T, .test = T,
                                .signif.label.size=3.5) {
  .data = rename_column(.data, "Volume", "Value")

  vis_bar(.data = .data, .by = .by, .meta = .meta,
          .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = F,
          .points = .points, .test = .test, .signif.label.size = .signif.label.size,
          .defgroupby = "Sample", .grouping.var = "Group",
          .labs = c(NA, "Clonotypes"),
          .title = "Number of clonotypes", .subtitle = "Number of unique clonotypes in the input data",
          .legend=NA, .leg.title = NA)
}

vis.immunr_exp_count <- function (.data, .by = NA, .meta = NA, .logx=T, .logy=T, .labs = c("Abundance", "Number of clonotypes"),
                                  .title = "Distribution of clonotype abundances", .subtitle = NULL,
                                  .legend=T, .leg.title = NA) {
  group_res = process_metadata_arguments(.data, .by, .meta)
  .data$Group = group_res$group_column

  if (is.na(.leg.title)) {
    if (group_res$is_grouped) {
      .leg.title = "Group"
    } else {
      .leg.title = "Sample"
    }
  }

  p = ggplot() + geom_line(aes(x = Clone.num, y = Clonotypes, col = Group, group = Sample), data = .data, stat = "identity")

  if (.logx) {
    p = p + scale_x_log10()
  }
  if (.logy) {
    p = p + scale_y_log10()
  }

  p = p +
    theme_linedraw() +
    labs(x = .labs[1], y = .labs[2], title = .title, subtitle = .subtitle, color = .leg.title)

  p
}

vis.immunr_exp_len <- function (.data, .by = NA, .meta = NA,
                                .errorbars = c(0.025, 0.975), .errorbars.off = F,
                                .points = T, .test = T,
                                .signif.label.size = 3.5) {
  .data = rename_column(.data, "Count", "Value")

  vis_bar(.data = .data, .by = .by, .meta = .meta,
          .errorbars = c(0.025, 0.975), .errorbars.off = .errorbars.off, .stack = F,
          .points = .points, .test = .test, .signif.label.size = .signif.label.size,
          .defgroupby = "Sample", .grouping.var = "Length",
          .labs = c("CDR3 length", "Clonotypes"),
          .title = "Distribution of CDR3 lengths", .subtitle = NULL,
          .legend=T, .leg.title = NA)
}


##### Kmer analysis & sequence logo plot #####


#' Most frequent kmers visualisation.
#'
#' @name vis.immunr_kmer_table
#'
#' @description
#' Plot a distribution (bar plot) of the most frequent kmers in a data.
#'
#' @param .kmers Data frame with two columns "Kmers" and "Count" or a list with such data frames. See Examples.
#' @param .head Number of the most frequent kmers to choose for plotting from each data frame.
#' @param .position Character vector of length 1. Position of bars for each kmers. Value for the \code{ggplot2} argument \code{position}.
#'
#' @seealso \code{get.kmers}
#'
#' @examples
#' \dontrun{
#' # Load necessary data and package.
#' library(gridExtra)
#' load('immdata.rda')
#' # Get 5-mers.
#' imm.km <- get.kmers(immdata)
#' # Plots for kmer proportions in each data frame in immdata.
#' p1 <- vis.kmer.histogran(imm.km, .position = 'stack')
#' p2 <- vis.kmer.histogran(imm.km, .position = 'fill')
#' grid.arrange(p1, p2)
#' }
#' @export
vis.immunr_kmer_table <- function (.data, .head = 100, .position = c('stack', 'dodge', 'fill'), .log = F) {
  .position = switch(substr(.position[1], 1, 1), s = "stack", d = "dodge", f = "fill")
  # .data[is.na(.data)] <- 0

  max_counts = apply(.data[,-1], 1, max, na.rm = T)
  max_indices = head(order(max_counts, decreasing = T), .head)
  n_samples = ncol(.data) - 1
  .data = .data[max_indices,]

  melted = reshape2::melt(.data, id.vars = "Kmer")
  colnames(melted) = c("Kmer", "Sample", "Count")
  p <- ggplot() + geom_bar(aes(x = Kmer, y = Count, fill = Sample), col = "black",
                           data = melted, stat = 'identity', position = .position) +
    theme_linedraw()
  if (.position[1] == 'stack' || .position[1] == 'dodge') {
    p <- p + ylab('Count') + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  } else {
    p <- p + ylab('Proportions') + theme(axis.text.x  = element_text(angle=90, vjust = .5))
  }

  if (.log) {
    p = p + scale_y_log10(expand = c(0, 0))
  } else {
    p = p + scale_y_continuous(expand = c(0, 0))
  }

  p + scale_x_discrete(expand = c(0, 0)) + .tweak_fill(n_samples) + ggtitle("Kmers distribution")
}

#' Logo - plots for amino acid and nucletide profiles.
#'
#' @name vis_logo
#'
#' @description
#' Plot logo-like graphs for visualising of nucleotide or amino acid motif sequences / profiles.
#'
#' @param .data Output from the \code{kmer.profile} function.
#' @param .replace.zero.with.na if T then replace all zeros with NAs, therefore letters with
#' zero frequency wont appear at the plot.
#' @param .jitter.width,.jitter.height,.dodge.width Parameters to \code{position_jitterdodge}
#' for aligning text labels of letters.
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' d <- kmer_profile(c('CASLL', 'CASSQ', 'CASGL'))
#' vis_logo(d)
#' }
#' @export vis.immunr_kmer_profile vis_logo
vis_logo <- function (.data, .replace.zero.with.na = T, .jitter.width = .01, .jitter.height = .01, .dodge.width = .15) {
  .data = reshape2::melt(.data)
  if (.replace.zero.with.na) {
    .data$value[.data$value == 0] = NA
  }
  colnames(.data) = c('AA', "Position", "Freq")
  ggplot(aes(x = Position, y = Freq, fill = AA, colour = AA), data = .data) +
    geom_point(colour = 'black') +
    geom_text(aes(label = AA), size = 5,
              position = position_jitterdodge(jitter.width = .jitter.width,
                                              jitter.height = .jitter.height,
                                              dodge.width = .dodge.width)) +
    xlab("Position") + ylab("Proportion") +
    theme_linedraw()
}

vis.immunr_kmer_profile <- function (...) {
  vis_logo(...)
}


##### Other & WIP visualisations #####


#' Repertoire dynamics
#'
#' @name vis.immunr_dynamics
#'
#' @export vis.immunr_dynamics
vis.immunr_dynamics <- function (.data, .by = NA, .meta = NA, .plot = c("area", "line"), .log = T) {
  .plot = .plot[1]

  proc = reshape2::dcast(.data, Sequence ~ Sample, value.var = c("Data"))
  logic = rowSums(is.na(proc)) == 0
  get_seq = proc$Sequence[logic]
  .data = subset(.data, .data$Sequence %in% get_seq)
  # print(reshape2::dcast(.data, Sequence ~ Sample, value.var = c("Data")))

  if (.plot == "area") {
    p = ggplot() + geom_area(aes(x = Sample, y = Hi - Lo, group = Sequence, fill = Sequence), color = "black", data = .data) + theme_linedraw()
  } else {
    p = ggplot() +
      geom_line(aes(x = Sample, y = Data, group = Sequence, color = Sequence), data = .data) +
      geom_errorbar(aes(x = Sample, ymin = Lo, ymax = Hi, y = Data, group = Sequence, color = Sequence), data = .data) +
      theme_linedraw()
  }

  if (.log) {
    p = p + scale_y_log10()
  }
  p
}

#' Work in progress
#'
#' @name vis.immunr_mutation_network
#'
#' @export vis.immunr_mutation_network
vis.immunr_mutation_network <- function (.data) {
  stop(IMMUNR_ERROR_NOT_IMPL)
}


#' Work in progress
#'
#' @name vis.immunr_cdr_prop
#'
#' @export vis.immunr_cdr_prop
vis.immunr_cdr_prop <- function (.data, by = NA, .meta = NA, .plot = c("box", "hist")) {
  stop(IMMUNR_ERROR_NOT_IMPL)
}
