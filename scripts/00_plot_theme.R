# =============================================================================
# 00_plot_theme.R — shared ggplot theme + genotype colors (source() from every R script)
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })

GENO_LEVELS <- c("WT", "Q157R")                       # WT is always the reference level
GENO_COLORS <- c(WT = "#3B6FB6", Q157R = "#C8553D")

theme_lr <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = base_size + 1),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "top")
}

save_plot <- function(p, path_noext, width = 6, height = 4) {
  ggsave(paste0(path_noext, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_noext, ".png"), p, width = width, height = height, dpi = 300)
  message("saved ", path_noext, ".{pdf,png}")
}
