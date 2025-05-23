library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(rlang)

## set command line arguments
args <- commandArgs(trailingOnly = TRUE)

# stop if no command line argument
if(length(args)==0){
  print("Nothing to do...")
  stop("Requires command line argument.")
}

## Read in data
basename <<- args[1]
zipfile <<- args[2]
minseqs <<- args[3]
plotting <<- args[4]

#' qc_positions
#'
#' @description
#' `qc_positions` converts the bins in positional fastqc data to numerical format
#'
#' @details
#' This function expands the positional column by expanding the numbers in the range of n-n+k
#' essentially the fastqc module bins the range of the data collapsing and in order to convert it
#' to a numerical range we must expand the dataframe. No new data is made with this function, it
#' merely expands the existing data to the cycle range included in the bin.
qc_positions <- function(df) {
  df$Position <- as.character(df$Position)
  df <- df %>%
    mutate(expansions = map(Position, function(Base) {
      if (grepl("-", Base)) {
        range <- as.numeric(str_extract_all(Base, "\\d+")[[1]])
        if (length(range) == 2) {
          return(seq(range[1], range[2]))
        }
      }
      return(as.numeric(Base))
    })) %>%
    unnest(expansions) %>%
    mutate(Position = as.numeric(expansions)) %>%
    select(-expansions)
  df <- df %>% drop_na()
}

#' module_data
#'
#' @description
#' `module_data` can be used for extracting the dataframes of data from the
#'  fastqc module
#'
#' @details
#' This function handles the main portion of the script extracting each of
#' the data portions within the fastqc module independently - this way we
#' can parse them into a R dataframe.
module_data <- function(df, module, basename) {
  if (module == "adapter") {
    start <- "^>>Adapter Content"
    columnnames <- c(
      "Position", "Illumina Universal Adapter",
      "Illumina Small RNA 3' Adapter",
      "Illumina Small RNA 5' Adapter",
      "Nextera Transposase Sequence",
      "PolyA", "PolyG"
    )
  } else if (module == "quality") {
    start <- "^>>Per base sequence quality"
    columnnames <- c(
      "Position",
      "Mean", "Median",
      "Lower.Quartile",
      "Upper.Quartile",
      "10th.Percentile",
      "90th.Percentile"
    )
  } else if (module == "per_quality") {
    start <- "^>>Per sequence quality scores"
    columnnames <- c("Quality", "Count")
  } else {
    print("Unrecognized module!")
    return()
  }

  end_module <- "^>>END_MODULE"

  start_line <- grep(start, df)
  end_line <- grep(end_module, df)


  for (i in end_line) {
    if (i > start_line) {
      end_line <- i
      break
    }
  }

  module_data <- df[(start_line + 1):(end_line - 1)]
  module_data <- as.data.frame(do.call(
    rbind,
    strsplit(
      module_data[2:(end_line - start_line)],
      "\t"
    )
  ))
  colnames(module_data) <- columnnames
  if (module == "per_quality") {
    module_data <- module_data %>% mutate_at(columnnames, as.numeric)
    module_data$basename <- basename
    return(module_data)
  }
  module_data <- qc_positions(module_data)
  module_data <- module_data %>% mutate_at(columnnames, as.numeric)
  module_data$basename <- basename
  module_data <- module_data %>%
    mutate(across(where(is.list), ~ sapply(.x, paste, collapse = ",")))
  return(module_data)
}


#' Check total seq count
#'
#' @description
#' `check_seq_count` returns TRUE if the number of seqs exceeds the minimum
#' it returns false otherwise
#'
#' @details
#' This function is used for parsing fastqc data to obtain the number of
#' sequences passed to the module. It can be used for filtering out files with
#' limited reads
check_seq_count <- function(df) {
  total_seqs <- "^Total Sequences"
  l_total_seq <- grep(total_seqs, df)
  seqs <- as.numeric(strsplit(df[l_total_seq], "\t")[[1]][2])
  if (seqs < minseqs) {
    cat("Not enough reads! ", seqs, "\n")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

qc_histogram <- function(df, basename, plot_group = NULL) {
  if (!is.null(plot_group)) {
    plot_group <- facet_wrap(paste0("~", plot_group))
  }
  p <- df %>% ggplot(aes(x = Quality, y = Count)) +
    annotate("rect", xmin = 30, xmax = 40, ymin = -Inf, ymax = Inf, fill = "green", alpha = 0.2) +
    annotate("rect", xmin = 20, xmax = 30, ymin = -Inf, ymax = Inf, fill = "yellow", alpha = 0.2) +
    annotate("rect", xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) +
    xlim(0, 40) +
    geom_col(width = 2, fill = TesseraGrey) +
    labs(x = "QScore", y = "Count", title = "Quality Score Histogram") +
    theme_light() +
    theme(
      legend.position = "none", text = element_text(size = 11),
      axis.text.x = element_text(angle = 90, hjust = 0)
    ) +
    plot_group
  ggsave(paste0(basename, ".seqplot.png")), p, h = 5, w = 7)
  return(p)
}


qc_boxplot <- function(df, basename, plot_group = NULL) {
  # Bin calculation
  df$Bin <- as.factor(ceiling(df$Position / 10) * 10)
  df <- df %>% filter(!is.na(Bin))
  if (!is.null(plot_group)) {
    group_syms <- syms(c("Bin", plot_group))
    df_summary <- df %>%
      group_by(!!!group_syms) %>%
      summarize(
        lower = first(Lower.Quartile),
        middle = first(Median),
        upper = first(Upper.Quartile),
        ymin = first(`10th.Percentile`),
        ymax = first(`90th.Percentile`),
        .groups = "drop"
      )
  } else {
    df_summary <- df %>%
      group_by(Bin) %>%
      summarize(
        lower = first(Lower.Quartile),
        middle = first(Median),
        upper = first(Upper.Quartile),
        ymin = first(`10th.Percentile`),
        ymax = first(`90th.Percentile`),
        .groups = "drop"
      )
  }
  p <- ggplot(df_summary, aes(x = Bin)) +
    annotate("rect", ymin = 30, ymax = 40, xmin = -Inf, xmax = Inf, fill = "green", alpha = 0.2) +
    annotate("rect", ymin = 20, ymax = 30, xmin = -Inf, xmax = Inf, fill = "yellow", alpha = 0.2) +
    annotate("rect", ymin = 0, ymax = 20, xmin = -Inf, xmax = Inf, fill = "red", alpha = 0.2) +
    ylim(0, 40) +
    geom_boxplot(
      aes(
        lower = lower,
        middle = middle,
        upper = upper,
        ymin = ymin,
        ymax = ymax,
        fill = if (!is.null(plot_group)) .data[[plot_group]] else NULL
      ),
      col = TesseraDarkGrey, fill = TesseraGrey, stat = "identity"
    ) +
    labs(x = "Base", y = "QScore", title = "Positional QScore") +
    theme_light() +
    theme(
      legend.position = if (!is.null(plot_group)) "right" else "none",
      text = element_text(size = 11),
      axis.text.x = element_text(angle = 90, hjust = 0)
    )
  if (!is.null(plot_group)) {
    p <- p + facet_wrap(as.formula(paste("~", plot_group)))
  }
  
  ggsave(paste0(basename, ".qcboxplot.png")), p, h = 5, w = 7)
  return(p)
}


qc_adapter_plot <- function(df, basename, plot_group = NULL) {
  p <- df %>% ggplot(aes(x = `Position`, y = `Illumina Universal Adapter`)) +
    annotate("rect", ymin = 0, ymax = 10, xmin = -Inf, xmax = Inf, fill = "green", alpha = 0.2) +
    annotate("rect", ymin = 10, ymax = 40, xmin = -Inf, xmax = Inf, fill = "yellow", alpha = 0.2) +
    annotate("rect", ymin = 40, ymax = 100, xmin = -Inf, xmax = Inf, fill = "red", alpha = 0.2) +
    ylim(0, 100) +
    theme_light() +
    labs(x = "Base", y = "Illumina Adapter %", title = "Adapter Content") +
    theme(text = element_text(size = 11),
          axis.text.x = element_text(angle = 90, hjust = 0))

  if (is.null(plot_group)) {
    p <- p +
      geom_smooth(stat = "smooth", se = FALSE, color = "black") +
      theme(legend.position = "none")
  } else {
    p <- p +
      geom_smooth(stat = "smooth", se = FALSE, aes(color = .data[[plot_group]])) +
      theme(legend.position = "right")}
  ggsave(paste0(basename, ".adapterplot.png")), p, h = 5, w = 7)
  return(p)
}

#' extract_qc_data
#'
#' @description
#' `extract_qc_data` is mainly a wrapper for unziping and extracting the
#' fastqc data from the zip file.
#'
#' @details
#' This fucnction passes the text data into the parser functions which handle
#' the main logic
extract_qc_data <- function(basename, zipfile, plotting = FALSE) {
  cat("Processing file:", zipfile, "\n")
  #   basename <- as.character(strsplit(basename(zipfile), "\\.zip")[1])
  data_name <- paste0(basename, "/fastqc_data.txt")
  zip_contents <- unzip(zipfile, list = TRUE)

  if (data_name %in% zip_contents$Name) {
    # parse lines
    qc_data <- readLines(unz(zipfile, data_name))
    if (!check_seq_count(qc_data, minseqs = minseqs)) {
      print("EXIT")
      return()
    }
    basename <- sapply(strsplit(basename, "_gDNA"), `[`, 1)

    # process each of the modules
    modules <- c("adapter", "quality", "per_quality")
    dfs <- setNames(
        lapply(modules, function(x) module_data(qc_data, x, basename)),
        modules
        )

    mapply(function(df, module) {
        write.csv(df, paste0(basename, "_", module, ".csv"))
    }, dfs, modules)


    if(plotting){
      p1 <- qc_histogram(df$per_quality, basename)
      p2 <- qc_adapter_plot(df$adapter, basename)
      p3 <- qc_boxplot(df$quality, basename)
      p4 <- plot_grid(
        p1, p2, p3,
        nrow = 2,
        align = "hv"
      )
      
      title <- ggdraw() +
        draw_label(
          basename,
          fontface = "bold",
          x = 0,
          hjust = 0
        ) +
        theme(
          plot.margin = margin(0, 0, 0, 7),
          plot.background = element_rect(fill = "white", color = NA)
        )
      p5 <- plot_grid(
        title, p4,
        ncol = 1,
        rel_heights = c(0.1, 1)
      )

      ggsave(paste0(basename, "_qc_plot.png")), p5)
    }
    print(paste0(basename, " complete."))
  } else {
    cat("fastqc_data.txt not found in", data_name, "\n")
  }
}

# exec
extract_qc_data(basename, zipfile, plotting)