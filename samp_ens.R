library(shiny)
library(ggplot2)
library(dplyr)
library(magrittr)
library(openxlsx)
library(RColorBrewer)
library(MASS)
library(plotly)

#' Summarize entropy values across a list of recordings
#'
#' @param lista A list of entropy data frames.
#' @param signature A character string identifying the analysis configuration.
#' @param label A character string used as a group label in the output.
#' @return A data frame with columns: r, SampEn, vars, group.
summarize_entropies <- function(lista, signature, label) {
  data.frame(r, SampEn = means[[1]], vars = vars[[1]], group = group, stringsAsFactors = FALSE)
}

#' Find intersection points between two SampEn curves
#'
#' Detects the x-positions where curve y2 crosses curve y1, and returns the
#' crossing coordinates for intersections within the physiologically relevant
#' tolerance range (0.03, 3.4).
#'
#' @param x Numeric vector. The shared independent variable (tolerance r).
#' @param y1 Numeric vector. SampEn values of the first recording.
#' @param y2 Numeric vector. SampEn values of the second recording.
#' @param stopme Logical. Reserved for debugging; currently unused.
#' @return A list of crossing point records. Each element is either a named list
#'   with entries \code{r}, \code{y1}, \code{y2} (coordinates at the crossing),
#'   or the string \code{"pusta"} for crossings outside the valid range.
crosses <- function(x, y1, y2, stopme) {
  if (length(y1) != length(y2)) {
    default_length <- min(length(y1), min(length(y2)))
    x <- x[1:default_length]
    y1 <- y1[1:default_length]
    y2 <- y2[1:default_length]
  }
  diffs <- y2 - y1 > 0
  crosski <- which(diff(diffs, na.rm = TRUE) != 0)
  lapply(crosski, function(itm) {
    if (x[itm] > 0.03 & x[itm] < 3.4) {
      list("r" = x[itm], "y1" = y1[itm], "y2" = y2[itm])
    } else {
      "pusta"
    }
  })
}

#' Load and preprocess SampEn data from a folder of recordings
#'
#' Reads comma- or tab-separated SampEn files, one per recording. For each
#' recording the tolerance axis (first column, in milliseconds) is normalised
#' by the recording's SDNN value read from \file{PoincrePlot.xlsx}, yielding
#' the dimensionless ratio \eqn{r = \rho / \mathrm{SDNN}}. Sample entropy is
#' computed as \eqn{\mathrm{SampEn} = \ln A - \ln B}, where A and B are the
#' template-match counts stored in columns 2 and 3 of each file.
#'
#' @param folder_list A nested list of file paths as returned by
#'   \code{\link{get_data_paths}}.
#' @param sep Character. Column separator used in the data files (default
#'   \code{"\t"}).
#' @param label Character. Passed through for downstream labelling; currently
#'   unused inside the function body.
#' @return A data frame with columns \code{r} (normalised tolerance),
#'   \code{SampEn}, and \code{group} (filename of the recording), containing
#'   only complete cases.
extract_data <- function(folder_list, sep = '\t', label) {
  wynik <- lapply(folder_list, function(subfolder) {
    lapply(subfolder, function(data_file) {
      a <- read.csv(file.path(data_file), sep = sep)
      if (ncol(a) == 1) {
        # Fallback: re-read with tab separator when auto-detection fails
        a <- read.csv(file.path(data_file), sep = '\t')
      }
      group <- strsplit(data_file, "/")[[1]]
      group <- group[length(group)]
      print(max(a[[1]] / 1000))
      HRV <- openxlsx::read.xlsx('PoincrePlot.xlsx')
      SDNN <- HRV[sub("_P.rea", "", HRV$file) == sub(".mat", "", group), ][["SDNN"]]
      print(paste("group", SDNN))
      data.frame(r = a[[1]] / SDNN, SampEn = (log(a[[2]]) - log(a[[3]])), group = group)
    })
  })
  wynik <- wynik[[1]] %>%
    dplyr::bind_rows()
  wynik[complete.cases(wynik), ]
}

#' Compute pairwise SampEn curve intersections and produce a summary plot
#'
#' Loads all recordings via \code{\link{extract_data}}, iterates over every
#' ordered pair of recordings to detect their SampEn curve crossings, and
#' returns both the crossing data and a ggplot2 figure showing all curves on a
#' log-scaled tolerance axis.
#'
#' @param folder_list A nested list of file paths as returned by
#'   \code{\link{get_data_paths}}.
#' @param sep Character. Column separator for the data files (default
#'   \code{"\t"}).
#' @param x_lim Numeric vector of length 2, or \code{NULL}. Optional x-axis
#'   limits for the plot.
#' @param y_lim Numeric vector of length 2, or \code{NULL}. Optional y-axis
#'   limits for the plot.
#' @param label Character. Group label prefix (default \code{"p"}).
#' @return A named list with two elements:
#'   \describe{
#'     \item{num_results}{A named list of crossing records, one entry per
#'       recording pair.}
#'     \item{plot}{A \code{ggplot} object showing all SampEn curves.}
#'   }
get_figures_and_results <- function(folder_list, sep = '\t', x_lim = NULL, y_lim = NULL, label = "p") {
  wynik_summix <- extract_data(folder_list, sep = sep, label = label)
  if (length(which(wynik_summix$r == 0)) > 0) {
    wynik_summix <- wynik_summix[-which(wynik_summix$r == 0), ]
  }
  p <- ggplot2::ggplot(data = wynik_summix, ggplot2::aes(r, SampEn, group = group)) +
    annotate("rect", xmin = 0.1, xmax = 0.25, ymin = -Inf, ymax = Inf, fill = "gray80") +
    geom_line() +
    scale_x_log10(
      limit = 'if'(is.null(x_lim), NULL, x_lim),
      breaks = function(x) sort(unique(c(scales::log_breaks()(x), 0.25))),
      labels = scales::label_number(accuracy = 0.01),
      name = expression(rho)
    ) +
    'if'(is.null(y_lim), NULL, ylim(y_lim)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  separate_lines <- split(wynik_summix, wynik_summix$group)
  results_list <- list()
  for (idx in 1:length(separate_lines)) {
    stopme <- FALSE
    for (idy in idx:length(separate_lines)) {
      intermediate_result <- crosses(x = separate_lines[[idx]]$r,
                                     y1 = separate_lines[[idx]]$SampEn,
                                     y2 = separate_lines[[idy]]$SampEn,
                                     stopme)
      if (length(intermediate_result) != 0) {
        results_list[paste0(separate_lines[[idx]]$group[[1]], ", ", separate_lines[[idy]]$group[[1]])] <- list(intermediate_result)
      }
    }
  }
  list(num_results = results_list, plot = p)
}

#' Collect file paths for all recordings in a data folder
#'
#' Recursively lists all files under \code{folder}, excluding any path that
#' contains the string \code{"README"}.
#'
#' @param folder Character. Path to the root data directory.
#' @return A list containing one character vector of absolute file paths.
get_data_paths <- function(folder) {
  folders <- list.dirs(path = folder, full.names = TRUE, recursive = TRUE) %>%
    lapply(function(subfolder) {
      if (!grepl("README", subfolder)) {
        list.files(subfolder, full.names = TRUE, recursive = FALSE)
      } else {
        return(NULL)
      }
    })

  folders[[1]][!grepl("README", folders[[1]])] %>%
    list()
}

#' Summarise crossing points into a tabular form
#'
#' Filters out out-of-range crossings (marked \code{"pusta"}) and assembles a
#' data frame containing, for each pair of recordings, up to two intersection
#' points (the first and second crossing) characterised by their tolerance
#' value \code{r} and the mean SampEn of the two curves at that point.
#'
#' @param num_results A named list of crossing records as returned by
#'   \code{\link{get_figures_and_results}}.
#' @return A data frame with columns: \code{crossing lines}, \code{r1},
#'   \code{SampEn} (at first crossing), \code{r2}, \code{SampEn} (at second
#'   crossing). Pairs with only one valid crossing have \code{NA} in the
#'   first-crossing columns.
summarize_crosses <- function(num_results) {
  result <- list()
  for (itm in names(num_results)) {
    intermediate <- Filter(function(elem) { elem[[1]] != "pusta" }, num_results[[itm]])
    if (length(intermediate) != 0) {
      result[[itm]] <- intermediate
    }
  }
  wynik <- data.frame()
  for (itm in names(result)) {
    if (length(result[[itm]]) >= 2) {
      intermediate_wynik <- data.frame(
        itm,
        result[[itm]][[2]]$r,
        (result[[itm]][[2]]$y1 + result[[itm]][[2]]$y2) / 2,
        result[[itm]][[1]]$r,
        (result[[itm]][[1]]$y1 + result[[itm]][[1]]$y2) / 2,
        stringsAsFactors = FALSE
      )
      names(intermediate_wynik) <- c("name", "r1", "sampEn1", "r2", "sampEn2")
      wynik <- rbind(wynik, intermediate_wynik, stringsAsFactors = FALSE)
    } else {
      # Single crossing: no return crossing is recorded
      intermediate_wynik <- data.frame(
        itm, NA, NA,
        result[[itm]][[1]]$r,
        (result[[itm]][[1]]$y1 + result[[itm]][[1]]$y2) / 2,
        stringsAsFactors = FALSE
      )
      names(intermediate_wynik) <- c("name", "r1", "sampEn1", "r2", "sampEn2")
      wynik <- rbind(wynik, intermediate_wynik, stringsAsFactors = FALSE)
    }
  }
  colnames(wynik) <- c("crossing lines", "r1", "SampEn", "r2", "SampEn")
  wynik
}

#' Draw an interactive 3-D scatter plot of crossing tolerance vs. HRV descriptor
#'
#' For each crossing point the HRV descriptor value (e.g. SDNN) of both
#' recordings is retrieved and plotted against the crossing tolerance \code{r}
#' in a three-dimensional scatter plot.
#'
#' @param table_crosses Data frame of crossings as returned by
#'   \code{\link{summarize_crosses}}.
#' @param table_hr Either a data frame or a path to an Excel file containing
#'   HRV descriptor values, with recording identifiers in the first column.
#' @param analyzed_column Character. Name of the HRV descriptor column to
#'   use on the y- and z-axes (default \code{"SDNN"}).
#' @return A \code{plotly} 3-D scatter plot object.
draw_3d_scatter <- function(table_crosses, table_hr, analyzed_column = "SDNN") {
  if (!is.data.frame(table_hr)) {
    table_hr <- read_table_hr(table_hr)
  }
  crosses_names <- read_crosses_names(table_crosses)
  r_data <- get_density_data(table_crosses, table_hr, crosses_names, analyzed_column)
  xax <- list(title = "r")
  yax <- list(title = paste0(analyzed_column, "1"))
  zax <- list(title = paste0(analyzed_column, "2"))
  plotly::plot_ly(x = r_data[[1]],
                  y = r_data[[2]],
                  z = r_data[[3]]) %>%
    plotly::add_markers(hovertemplate = paste(
      " r: %{x}<br>",
      paste0(analyzed_column, "1: %{y}<br>"),
      paste0(analyzed_column, "2: %{z}<br>"),
      "<extra></extra>"
    )) %>%
    layout(scene = list(xaxis = xax, yaxis = yax, zaxis = zax))
}

#' Draw a 2-D density heatmap of crossing tolerance vs. HRV descriptor ratio
#'
#' Estimates the joint density of the crossing tolerance \code{r} and the ratio
#' of the HRV descriptor for the two recordings at each crossing, and renders
#' it as a filled contour image. Supports both static PNG output and an
#' interactive plotly contour plot.
#'
#' @param table_crosses Data frame of crossings as returned by
#'   \code{\link{summarize_crosses}}.
#' @param table_hr Either a data frame or a path to an Excel file containing
#'   HRV descriptor values.
#' @param analyzed_column Character. Name of the HRV descriptor column
#'   (default \code{"SDNN"}).
#' @param plotly Logical. If \code{TRUE}, return an interactive plotly contour
#'   plot; otherwise produce a base-graphics image (default \code{FALSE}).
#' @param contours_no Integer. Number of grey levels in the colour ramp for the
#'   static plot (default \code{10}).
#' @param filename Character or \code{NULL}. If provided, the static plot is
#'   saved to this PNG file; otherwise it is drawn to the active graphics device.
#' @param ylab Character or \code{NULL}. Custom y-axis label. When \code{NULL}
#'   defaults to \code{"<analyzed_column>1/<analyzed_column>2"}.
#' @return Invisibly \code{NULL} for the static variant; a \code{plotly} object
#'   for the interactive variant.
draw_heatmap <- function(table_crosses,
                         table_hr,
                         analyzed_column = "SDNN",
                         plotly = FALSE,
                         contours_no = 10,
                         filename = NULL,
                         ylab = NULL) {
  if (!is.data.frame(table_hr)) {
    table_hr <- read_table_hr(table_hr)
  }
  crosses_names <- read_crosses_names(table_crosses)
  r_data <- get_density_data(table_crosses, table_hr, crosses_names, analyzed_column)
  # Greyscale colour ramp from white to black
  rf <- colorRampPalette(c("grey100", "grey0"))
  r <- rf(contours_no)
  x <- r_data[[1]]
  y <- r_data[[2]] / r_data[[3]]
  df <- data.frame(x, y)
  df <- subset(df, df$y <= 5)
  df <- df[complete.cases(df), ]
  df <- df[!is.infinite(df$y) & df$y > 0, ]
  if (nrow(df) < 2) {
    return(plotly::plotly_empty())
  }
  # 2-D kernel density estimation on a 500 x 500 grid
  k <- kde2d(df$x, df$y, n = 500)
  if (plotly) {
    x_axis <- list(title = "r")
    y_axis <- list(title = paste0(analyzed_column, "1/", analyzed_column, "2"))
    plot_ly(x = k$x,
            y = k$y,
            z = t(k$z),
            type = "contour",
            showscale = FALSE,
            hovertemplate = paste(
              "%{xaxis.title.text}: %{x}<br>",
              "%{yaxis.title.text}: %{y}<br>",
              "density: %{z}<br>",
              "<extra></extra>"
            )) %>%
      layout(xaxis = x_axis, yaxis = y_axis)
  } else {
    default_ylab <- paste0(analyzed_column, "1/", analyzed_column, "2")
    if (!is.null(filename)) {
      png(filename, width = 1400, height = 1200, res = 300)
      par(mgp = c(2.5, 1, 0))
      image(k, col = r, xlab = "r", ylab = if (is.null(ylab)) default_ylab else ylab)
      dev.off()
    } else {
      image(k, col = r, xlab = "r", ylab = if (is.null(ylab)) default_ylab else ylab)
    }
  }
}

#' Read an HRV descriptor table from an Excel file
#'
#' @param filename Character. Path to the \code{.xlsx} file.
#' @return A data frame containing the HRV descriptor table.
read_table_hr <- function(filename) {
  openxlsx::read.xlsx(filename)
}

#' Extract the pair of recording names from a crossings table
#'
#' Parses the \code{"crossing lines"} column of a crossings table, which stores
#' entries of the form \code{"recordingA.mat, recordingB.mat"}, and returns a
#' two-column matrix of clean recording identifiers (without the \code{.mat}
#' extension or leading spaces).
#'
#' @param table_crosses Data frame of crossings as returned by
#'   \code{\link{summarize_crosses}}.
#' @return A character matrix with one row per crossing pair and two columns
#'   (\code{[,1]}: first recording, \code{[,2]}: second recording).
read_crosses_names <- function(table_crosses) {
  all_names <- c(one = c(), two = c())
  for (name in table_crosses[[1]]) {
    individual_names <- strsplit(name, ",")[[1]]
    individual_names <- c(sub(" ", "", sub(".mat", "", individual_names[[1]])),
                          sub(" ", "", sub(".mat", "", individual_names[[2]])))
    all_names <- rbind(all_names, individual_names)
  }
  rownames(all_names) <- NULL
  all_names
}

#' Retrieve tolerance and HRV descriptor values for each crossing pair
#'
#' For every row in the crossings table, looks up the HRV descriptor value for
#' each of the two recordings involved and assembles a long-format data frame
#' with one row per crossing event (two rows per pair, corresponding to the
#' two crossing points stored in \code{r1} and \code{r2}).
#'
#' @param table_crosses Data frame of crossings as returned by
#'   \code{\link{summarize_crosses}}.
#' @param table_hr Data frame of HRV descriptors with recording identifiers in
#'   the first column.
#' @param crosses_names Character matrix as returned by
#'   \code{\link{read_crosses_names}}.
#' @param analyzed_column Character. Name of the HRV descriptor column to
#'   retrieve.
#' @return A data frame with columns \code{r},
#'   \code{<analyzed_column>1}, and \code{<analyzed_column>2}, containing
#'   only complete cases.
get_density_data <- function(table_crosses, table_hr, crosses_names, analyzed_column) {
  r_s <- data.frame(r = c(), indep_1 = c(), indep_2 = c())
  for (line in seq(nrow(table_crosses))) {
    indep_var_location_1 <- grepl(crosses_names[line, ][[1]], table_hr[[1]])
    indep_var_location_2 <- grepl(crosses_names[line, ][[2]], table_hr[[1]])
    indep_var_1 <- table_hr[indep_var_location_1, analyzed_column]
    indep_var_2 <- table_hr[indep_var_location_2, analyzed_column]
    r_s <- rbind(r_s, c(table_crosses[line, 'r1'], indep_var_1, indep_var_2))
    r_s <- rbind(r_s, c(table_crosses[line, 'r2'], indep_var_1, indep_var_2))
  }
  r_s <- r_s[complete.cases(r_s), ]
  names(r_s) <- c('r', paste0(analyzed_column, "1"), paste0(analyzed_column, "2"))
  r_s
}
