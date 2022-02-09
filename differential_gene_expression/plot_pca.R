# plot pca for maha
## load packages
library(dplyr)
library(tidyr)
library(tidyverse)
library(sleuth)
library(ggplot2)
library(ggfortify)
library(cowplot)

base_dir <- "~/Downloads/maha/"
setwd(base_dir)

sample_id <- dir(file.path(base_dir, "mRNA"))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, 'mRNA', id))
s2c <- data.frame(sample = sample_id, 
                  host = rep(c("CC7", "H2", "CC7"), each = 8), 
                  symbiont = rep(c("A", "B", "B"), each = 8),
                  temp = rep(c("25", "32"), each = 4))
s2c$condition <- paste(s2c$host, s2c$symbiont, s2c$temp, sep = "")
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c_cc7 <- s2c[c(1:8, 17:24),]

so_all <- sleuth::sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so_all <- sleuth::sleuth_fit(so_all, ~host + symbiont + temp, "full")
so_all <- sleuth::sleuth_fit(so_all, ~1, "reduced")
so_all <- sleuth::sleuth_lrt(so_all, "reduced", "full")

so_cc7 <- sleuth::sleuth_prep(s2c_cc7, extra_bootstrap_summary = TRUE)
so_cc7 <- sleuth::sleuth_fit(so_cc7, ~condition, "full")
so_cc7 <- sleuth::sleuth_fit(so_cc7, ~1, "reduced")
so_cc7 <- sleuth::sleuth_lrt(so_cc7, "reduced", "full")


so_kt_fil_norm <- kallisto_table(so_all, use_filtered = TRUE)
tpm_fn <- select(so_kt_fil_norm, target_id, tpm, sample) %>% spread(target_id, tpm)
prc <- prcomp(tpm_fn[, -1], scale. = TRUE)
autoplot(prc)
tpm_fn$condition <- rep(c("CC7 25", "CC7 32", "H2 25", "H2 32", "CC7-B01 25", "CC7-B01 32"), each = 4)
p1 <- autoplot(prc, data = tpm_fn, colour = "condition")
cols <- c("#99d8ca", "#2aa15f", "#fcae6b", "#e6550b", "#9fcae0", "#3082bd")
p1 <- p1 + 
  scale_color_manual(values = cols) +
  theme_classic()
p1

cc7_kt_fn <- kallisto_table(so_cc7, use_filtered = TRUE)
cc7_tpm_fn <- select(cc7_kt_fn, target_id, tpm, sample) %>% spread(target_id, tpm)
cc7_prc <- prcomp(cc7_tpm_fn[,-1], scale. = TRUE)
autoplot(cc7_prc)
cc7_tpm_fn$condition <- rep(c("CC7 25", "CC7 32", "CC7-B01 25", "CC7-B01 32"), each = 4)
p2 <- autoplot(cc7_prc, data = cc7_tpm_fn, colour = "condition")
cc7_cols <- c("#99d8ca", "#2aa15f", "#fcae6b", "#e6550b")
p2 <- p2 + 
  scale_color_manual(values = cc7_cols) +
  theme_classic()
p2

ppr <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2])
ppc <- plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2])

ppr

ppc

sessionInfo()
