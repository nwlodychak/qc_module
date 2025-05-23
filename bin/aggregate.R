library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(rlang)
source("qc_module.R")

## set command line arguments
args <- commandArgs(trailingOnly = TRUE)

# stop if no command line argument
if(length(args)==0){
  print("Nothing to do...")
  stop("Requires command line argument.")
}

## Read in data
samplesheet <<- read.csv(args[1])
plot_group <<- args[2]

pat <<- c(".*_adapter_df.csv$", ".*per_quality_df.csv$", ".*global_quality_df.csv$")

#' combine_by_factor
#'
#' @description
#' `combine_by_factor` converts all the dataframes into a single dataframe that can be collapsed by factor
#'
#' @details
#' This function collects the dataframes from each of the previous outputs in the prior step. For each of the
#' following plots, adatper, quality and qc_histogram, we can aggregate them together into a single plot
combine_by_factor <- function(submission_id, value, plot_group = NULL){
  for (p in pat) {
    files <- list.files(path = file.path("output", submission_id), pattern = paste0("^", submission_id, p, sep = ""), full.names = TRUE)
    df <- files %>% map_df(read_csv, show_col_types = FALSE)
    plot_df <- sample_sheet %>% left_join(df)
    
    
    plot_df[[value]] <- as.factor(plot_df[[value]])
    if(p == ".*_adapter_df.csv$"){
      p1 <- qc_adapter_plot(plot_df, value, out = "plots", plot_group = value)
    }
    else if (p == ".*per_quality_df.csv$") {
      p2 <- qc_histogram(plot_df, value, out = "plots", plot_group = plot_group)
    }
    else {
      p3 <- qc_boxplot(plot_df, value, out = "plots", plot_group = plot_group)
    }
  }
  
  pvoid <- ggplot() + theme_void()
  p4 <- plot_grid(
    p1, p2, p3, pvoid,
    nrow = 2,
    align = "hv"
  )
  
  title <- ggdraw() +
    draw_label(
      submission_id,
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
  ) + theme(plot.background = element_rect(fill = "white"))
  ggsave(file.path("plots", submission_id ,paste0(submission_id, "_qc_plot.png")), p5, h = 10, w = 15)
}

combine_by_factor(submission_id, "prep")