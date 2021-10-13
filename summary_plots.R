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
    aes(x = rep_grp, y = msr, color = scaling_group, fill = scaling_group) + 
    geom_point(data = scaling[scaling$is_pool == FALSE, ]) +
    geom_bar(data = scaling[scaling$is_pool == TRUE, ], stat = "identity", alpha = 0.5) +
    style_minute_barplot(max_y = max(scaling$msr)) +
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
    aes(x = sample_name, y = msr, color = scaling_group, fill = scaling_group) + 
    geom_bar(data = scaling[scaling$is_pool == FALSE, ], stat = "identity", alpha = 0.5) +
    style_minute_barplot(max_y = max(scaling$msr)) 
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
  scaling %>% 
    mutate(inr = X.reads / n_input_reads,
           rep_grp = gsub("_(pooled|rep[1-9]).*", "", sample_name),
           is_pool = grepl("pooled", sample_name)) %>% 
    group_by(scaling_group) %>% 
    mutate(msr = inr / first(inr))
}


#' Appearance ggplot functions common to minute scaled barplots
#'
#' @param max_y Maximum value on the y axis.
#'
#' @return A list of ggproto objects
style_minute_barplot <- function(max_y) {
  list(theme_classic(base_size = 10),
       facet_wrap(~scaling_group, scales = "free_x", ncol = 2),
       geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.4),
       theme(legend.position = "None",
             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)),
       scale_y_continuous(breaks = seq(0, max_y, by = 0.2)),
       labs(title = "MINUTE-ChIP scaled global read levels",
            x = "Sample",
            y = "Minute-ChIP Scaled Fraction"))
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

single_width <- 6

# Account for longer names
single_height <- 8 

panel_width <- single_width * 2
panel_height <- ceiling(ngroups / 2) * single_height


ggsave(snakemake@output[[1]],
       plot = minute_scaled_replicates_barplot(scalinginfo),
       width = panel_width,
       height = panel_height,
       dpi = 300,
       units = "cm")

ggsave(snakemake@output[[2]],
       plot = minute_scaled_grouped_barplot(scalinginfo),
       width = panel_width,
       height = panel_height,
       dpi = 300,
       units = "cm")

