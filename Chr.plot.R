#!/usr/bin/env Rscript

# ------------------------
# Load libraries
# ------------------------
library(data.table)   # Efficient data import
library(ggplot2)      # Plotting
library(scales)       # For scaling axes
library(methods)      # Required for Rscript execution in some environments

# ------------------------
# Parse input
# ------------------------
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]  # Sample name prefix

# ------------------------
# Load coverage deviation file
# ------------------------
df4 <- fread(paste0(prefix, "_cov_dev_1Mb.tsv"), header = TRUE, data.table = FALSE)
colnames(df4) <- c("chr", "pos", "nread")
df4 <- df4[!grepl("unloc|^chrUn|^HAP1|^scaffold", df4$chr, ignore.case = TRUE), ]
df4$color <- "#386cb0"  # Default colour (blue)

# ------------------------
# Prepare chromosome sorting keys
# ------------------------
extract_num  <- function(x) gsub("[^0-9]", "", x)
extract_suf  <- function(x) gsub("chr|Chr|[0-9]", "", x)
suffix_order <- c("A","B","D","T","X","Z","")

df4$chr_number   <- extract_num(df4$chr)
df4$chr_suffix   <- extract_suf(df4$chr)
df4$chr_sort_key <- paste0(sprintf("%02d", as.numeric(df4$chr_number)),
                           match(df4$chr_suffix, suffix_order))
ord <- df4[order(df4$chr_sort_key), "chr"]
df4$chr <- factor(df4$chr, levels = unique(ord))

# ------------------------
# Identify wild relative chromosomes (ending in 'T') for comparison groups
# ------------------------
valid_T_chrs <- character(0)
T_chr_groups <- list()

for (cname in levels(df4$chr)) {
  if (grepl("T$", cname)) {
    sub <- df4[df4$chr == cname, ]
    nread_vals <- sub$nread
    rle_vals <- rle(nread_vals > 0.5)
    ends <- cumsum(rle_vals$lengths)
    starts <- c(1, head(ends, -1) + 1)
    long_runs <- which(rle_vals$values & rle_vals$lengths >= 5)
    if (length(long_runs) > 0) {
      valid_T_chrs <- c(valid_T_chrs, cname)
      grp <- gsub("[^0-9]", "", cname)
      T_chr_groups[[grp]] <- unique(c(T_chr_groups[[grp]], cname))
    }
  }
}

# ------------------------
# Function to flag candidate introgressions by colouring runs of interest red
# ------------------------
color_runs <- function(df, cname) {
  rows <- which(df$chr == cname)
  sub  <- df[rows, ]
  red  <- integer(0)

  if (grepl("T$", cname)) {
    # For wild relative chromosomes
    i <- 1
    while (i <= (nrow(sub) - 4)) {
      # Start red run if 5 consecutive bins have nread > 0.2
      if (all(sub$nread[i:(i+4)] > 0.2)) {
        run_start <- i
        i <- i + 5
        # End the run after 10 consecutive bins ≤ 0.2
        while (i <= nrow(sub)) {
          if (sub$nread[i] <= 0.2) {
            if (i + 9 <= nrow(sub) && all(sub$nread[i:(i+9)] <= 0.2)) break
          }
          i <- i + 1
        }
        run_end <- i - 1
        red <- c(red, run_start:run_end)
      } else {
        i <- i + 1
      }
    }
  } else {
    # For wheat chromosomes: mark drops < 0.2 for 5+ bins (if matched in T group)
    grp <- gsub("[^0-9]", "", cname)
    if (grp %in% names(T_chr_groups)) {
      cond <- sub$nread < 0.2
      rle_c <- rle(cond)
      ends <- cumsum(rle_c$lengths)
      starts <- c(1, head(ends, -1) + 1)
      good <- which(rle_c$values & rle_c$lengths >= 5)
      red <- unlist(lapply(good, function(i) starts[i]:ends[i]))
    }
  }

  if (length(red) > 0) df$color[rows[red]] <- "red"
  df
}

# ------------------------
# Apply red colouring logic to all chromosomes
# ------------------------
for (cname in levels(df4$chr)) {
  df4 <- color_runs(df4, cname)
}

# ------------------------
# Plot coverage deviation with red highlighting for candidate introgressions
# ------------------------
p <- ggplot(df4, aes(pos, nread)) +
  geom_point(aes(color = color), size = 0.5, alpha = 0.6) +  # Adjust alpha/size as needed
  scale_color_identity() +
  facet_wrap(~ chr, ncol = 4, strip.position = "left") +
  scale_x_continuous(
    labels = function(x) x / 1e6,
    breaks = seq(0, 9e8, by = 1e8),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 2),
    oob = squish,
    expand = c(0, 0)
  ) +
  xlab("\nGenomic Position (Mb)") +
  ylab("Coverage deviation\n") +
  theme(
    panel.grid.major.x = element_line(size = 0.08, linetype = "solid", colour = "gray68"),
    panel.grid.major.y = element_line(size = 0.08, linetype = "solid", colour = "gray68"),
    panel.background = element_rect(fill = "white"),
    panel.spacing = unit(1, "lines"),
    strip.text.y.right = element_text(angle = 90, face = "bold", size = 10),
    panel.border = element_rect(colour = "black", size = 0.4, fill = NA),
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  )

# ------------------------
# Save plot as PDF
# ------------------------
ggsave(paste0(prefix, "_Mut.pdf"), p,
       width = 15, height = 12, units = "in", limitsize = FALSE)
