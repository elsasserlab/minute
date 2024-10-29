suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))


#' Create a Minute summary barplot grouped by pooled + replicates. Each pooled
#' sample is a bar in this figure, corresponding replicates are points.
#' 
#' The figure is faceted according to scaling group.
#'
#' @param scaling_file scalinginfo.txt out of the minute pipeline
#'
#' @return A ggplot
minute_scaled_grouped_barplot <- function(scaling_file) {
  scaling <- read.table(scaling_file, sep="\t", header = T, comment.char = "")
  scaling <- calculate_ratios_and_groups(scaling)
  
  ggplot(data = scaling) + 
    aes(x = replace_delims_with_spaces(rep_grp), y = msr, color = scaling_group, fill = scaling_group) +
    geom_point(data = scaling[scaling$is_pool == FALSE, ]) +
    geom_bar(data = scaling[scaling$is_pool == TRUE, ], stat = "identity", alpha = 0.5) +
    style_minute_barplot() +
    theme(legend.position = "none") +
    scale_x_discrete(labels = scales::label_wrap(20)) +
    labs(subtitle = "Points - Replicates; Bars - Pooled")
}


#' Create a Minute summary barplot excluding the pooled samples. 
#' 
#' The figure is faceted according to scaling group.
#'
#' @param scaling_file scalinginfo.txt out of the minute pipeline
#'
#' @return A ggplot
minute_scaled_replicates_barplot <- function(scaling_file) {
  scaling <- read.table(scaling_file, sep="\t", header = T, comment.char = "")
  scaling <- calculate_ratios_and_groups(scaling)
  ggplot(data = scaling) + 
    aes(x = replace_delims_with_spaces(rep_grp), y = msr, fill = as.factor(rep_number)) +
    geom_bar(
      data = scaling[scaling$is_pool == FALSE, ],
      stat = "identity",
      alpha = 0.9,
      position = position_dodge2(preserve = "single"),
      color = "#555555",
      linewidth = 0.2
    ) +
    style_minute_barplot() +
    theme(legend.position = "bottom") +
    scale_x_discrete(labels = scales::label_wrap(20)) +
    labs(fill = "Replicate") +
    scale_fill_discrete()
}


#' Stacked barplot with number of reads per replicate pool, including Input
#'
#' @param scaling_file scalinginfo.txt out of the minute pipeline
#'
#' @return A ggplot showing barcode representation
barcode_representation_barplot <- function(scaling_file, percent = FALSE) {
  scaling <- read.table(scaling_file, sep="\t", header = T, comment.char = "")
  scaling <- calculate_ratios_and_groups(scaling)

  df_input <- scaling %>%
    select(scaling_group, n_input_reads, condition, rep_grp, is_pool) %>%
    filter(scaling_group == scaling$scaling_group[[1]]) %>%
    mutate(scaling_group = "Input") %>%
    rename(X.reads = n_input_reads)

  df_combined <- rbind(
    scaling %>% select(scaling_group, X.reads, condition, rep_grp, is_pool),
    df_input
  ) %>% filter(!is_pool) %>%
    group_by(scaling_group) %>%
    mutate(perc = (X.reads / sum(X.reads)) * 100) %>%
    ungroup()

  if (percent == TRUE) {
    stacked_replicate_groups_plot(df_combined, "perc") +
      style_barcode_representation() +
      labs(y="Reads (%)")
  } else {
    stacked_replicate_groups_plot(df_combined, "X.reads") +
      scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix="M")) +
      style_barcode_representation() +
      labs(y="Final mapped reads (millions)")
  }
}


#' Calculate input-normalized ratio and group scaled fraction out of a table
#' containing scalinginfo.txt contents.
#' 
#' It annotates the replicate groups and pools for easier plotting.
#'
#' @param scaling Scaling table
#'
#' @return A tibble with new columns inr, msr, rep_group and is_pool
calculate_ratios_and_groups <- function(scaling) {
  extract_rep_number <- function(str) {
    as.integer(gsub(".+[_rep|pooled]([0-9]*)\\.[[:alnum:]]+", "\\1", str))
  }
  scaling %>% 
    mutate(inr = X.reads / n_input_reads,
           rep_grp = gsub("_(pooled|rep[1-9]).*", "", sample_name),
           is_pool = grepl("pooled", sample_name),
           rep_number = extract_rep_number(sample_name),
           reference = gsub(".+\\.([[:alnum:]]+)", "\\1", sample_name)) %>%
    group_by(scaling_group) %>% 
    mutate(msr = inr / first(inr)) %>%
    rowwise() %>%
    mutate(condition = gsub(paste0(scaling_group, "_"), "", rep_grp)) %>%
    ungroup()
}


#' Appearance ggplot functions common to minute scaled barplots
#'
#' @return A list of ggproto objects
style_minute_barplot <- function() {
  list(theme_classic(base_size = 8),
       facet_wrap(~scaling_group, scales = "free_x", ncol = 2),
       geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.4),
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)),
       labs(title = "MINUTE-ChIP scaled global read levels",
            x = "Sample",
            y = "Minute-ChIP Scaled Fraction"))
}


#' GGplot style options for the barcode representation plot
#'
#' @return List of ggproto elements
style_barcode_representation <- function() {
  list(
    coord_flip(),
    theme_minimal(base_size = 8),
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks.x = element_line(),
          axis.ticks.length.x = grid::unit(2, "mm"),
          legend.position = "bottom"),
    labs(title = "Barcode representation",
         x = "",
         colour = "Condition",
         fill = "Condition")
  )
}

#' Replace dashes, underscores and points with spaces for labelling
#' @param s String to clean up
#' @return A string
replace_delims_with_spaces <- function(s) {
  gsub("_|-|\\.", " ", s)
}

#' Split a string into shorter strings delimited by \n
#'
#' @param s String to split
#' @param width Max width
#' @return A string
wrap_label <- function(s, width) {
  s <- sapply(s, function(x) {paste0(strwrap(x, width = width), collapse = "\n")})
  names(s) <- NULL
  s
}

#' Main part of the stacked barcode representation plot.
#' Overlays the pool values with the replicate groups to have separate lines
#' within each  bar chunk
#'
#' @param df_combined Table with values
#' @param value_column Column to use (perc or X.reads)
#' @return A ggplot object
stacked_replicate_groups_plot <- function(df_combined, value_column) {
  # Avoid the summarise warning, that prints on the snakemake workflow
  # using .groups is experimental, so it is probably better to silence
  options(dplyr.summarise.inform = FALSE)

  df_combined <- df_combined %>%
    mutate(condition = replace_delims_with_spaces(condition))

  df_blocks <- df_combined %>%
    group_by(scaling_group, rep_grp, condition) %>%
    summarise(total = sum(.data[[value_column]]))

  ggplot(
    df_combined,
    aes(x=scaling_group, y=!!sym(value_column), fill=wrap_label(condition, 20), color=wrap_label(condition, 20))
  ) +
    geom_bar(stat="identity", color="white", linewidth=0.2, alpha = 0.8) +
    geom_bar(
      data=df_blocks,
      aes(x=scaling_group, y=total),
      stat="identity",
      color="white",
      alpha = 0.5,
      linewidth = 1.2
    )
}

#' Calculate number of groups in a scalinginfo.txt file
#' 
#' @param scaling_file Path to scalinginfo.txt file
get_scaling_groups_number <- function(scaling_file) {
  scaling <- read.table(scaling_file, sep="\t", header = T, comment.char = "")
  length(levels(as.factor(scaling$scaling_group)))
}

scalinginfo <- snakemake@input[[1]]
ngroups <- get_scaling_groups_number(scalinginfo)

single_width <- 5

# Account for longer names
single_height <- 9

panel_width <- single_width * 2
panel_height <- ceiling(ngroups / 2) * single_height

ggsave(snakemake@output[[1]],
       plot = minute_scaled_replicates_barplot(scalinginfo),
       width = panel_width,
       height = panel_height,
       dpi = 150,
       units = "cm")

ggsave(snakemake@output[[2]],
       plot = minute_scaled_replicates_barplot(scalinginfo),
       width = panel_width,
       height = panel_height,
       dpi = 300,
       units = "cm")

ggsave(snakemake@output[[3]],
       plot = minute_scaled_grouped_barplot(scalinginfo),
       width = panel_width,
       height = panel_height,
       dpi = 150,
       units = "cm")

ggsave(snakemake@output[[4]],
       plot = minute_scaled_grouped_barplot(scalinginfo),
       width = panel_width,
       height = panel_height,
       dpi = 300,
       units = "cm")

ggsave(snakemake@output[[5]],
       plot = barcode_representation_barplot(scalinginfo),
       width = 12,
       height = 7,
       dpi = 150,
       units = "cm",
       bg = "white")

ggsave(snakemake@output[[6]],
       plot = barcode_representation_barplot(scalinginfo),
       width = 12,
       height = 7,
       dpi = 300,
       units = "cm")

ggsave(snakemake@output[[7]],
       plot = barcode_representation_barplot(scalinginfo, percent = TRUE),
       width = 12,
       height = 7,
       dpi = 150,
       units = "cm",
       bg = "white")

ggsave(snakemake@output[[8]],
       plot = barcode_representation_barplot(scalinginfo, percent = TRUE),
       width = 12,
       height = 7,
       dpi = 300,
       units = "cm")
