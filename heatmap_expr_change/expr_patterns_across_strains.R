#!/usr/bin/env Rscript

"> expr_patterns_across_strains.R <

For the set of 731 symbiosis-related genes (Cui et al., 2019), check whether
the three tested strains (CC7 / CC7-B01 / H2) behave differently.

Behaviour is coded three ways.

concordant: genes that are expressed post-symbiosis and post-heat stress in the 
            same direction (upreg|upreg + downreg|downreg)
not_diff:   genes that are not differentially expressed post-heat stress
discordant: genes that are expressed post-symbiosis and post-heat stress in the 
            different direction (upreg|downreg + downreg|upreg)
" -> doc

# install.packages('chisq.posthoc.test')
library(chisq.posthoc.test)

# build matrix from raw files
up_df <- read.table('clustermap.symbiosis_genes.up.tsv')
up_concordant <- colSums(up_df > 0)
up_concordant
up_discordant <- colSums(up_df < 0)
up_discordant

down_df <- read.table('clustermap.symbiosis_genes.down.tsv')
down_concordant <- colSums(down_df < 0)
down_concordant
down_discordant <- colSums(down_df > 0)
down_discordant

# rbind concordants and discordants
overall_matrix <- rbind(up_concordant + down_concordant,
                        up_discordant + down_discordant)
# not_diff is simply 731 minus the other two numbers, as all categories sum to 731
overall_matrix <- rbind(overall_matrix, 731 - colSums(overall_matrix))
rownames(overall_matrix) <- c('concordant', 'discordant', 'not_diff')
overall_matrix

# finally, perform chi-square test with bonferroni correction (default)
chisq.posthoc.test(overall_matrix, method='bonferroni')

sessionInfo()
