###############################################################################
# Microbiome analysis pipeline (CLEANED, same logic; BORAL and t-SNE removed)
###############################################################################

############################
# 0) Working directories
############################
setwd("C:/Users/orsia/Dropbox/Microbiome ms/stat")
setwd("/Users/lovas-kissadam/Library/CloudStorage/Dropbox/Ákos dropbox/stat")
getwd()

############################
# 1) Packages
############################
suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
  library(stringr)
  library(patchwork)
  library(cowplot)
  library(scales)
library(lme4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(geometry)
library(rstatix)
library(stringr)
library(patchwork)
library(cowplot)
library(scales)
library(forcats)
library(mvabund)
library(ggrepel)
library(readr)
})

############################
# 2) Load data
############################
mvb <- read.csv2("merged_wider_javított2.csv", stringsAsFactors = FALSE)

###############################################################################
# SECTION A — PERMANOVA, dispersion, SIMPER, NMDS (2D & 3D)
###############################################################################

mvb$.__tot <- rowSums(mvb[, 3:1111])
mvb <- subset(mvb, .__tot > 0, select = -.__tot)

pmv <- adonis2(mvb[,-c(1:2)]^0.25 ~ as.factor(BIRD), data = mvb, permutations = 999, method = "bray")
print(pmv)

Hell_Seeds <- decostand(mvb[,-c(1:2)], method = "hellinger")

bd <- betadisper(vegdist(Hell_Seeds), group = mvb$BIRD)
print(bd)
print(anova(bd))
print(TukeyHSD(bd))

sim <- simper(mvb[, -c(1:2)]^0.25, group = mvb$BIRD, permutations = 999)

sim_sum <- summary(sim)
full_df <- do.call(rbind, lapply(names(sim_sum), function(comp) {
  df <- as.data.frame(sim_sum[[comp]])
  df$Species <- rownames(sim_sum[[comp]])
  df$Comparison <- comp
  rownames(df) <- NULL
  df[, c("Comparison","Species", setdiff(names(df), c("Comparison","Species")))]
}))
write.csv(full_df, "SIMPER_full_results.csv", row.names = FALSE)

sig_list <- lapply(sim, function(x) {
  df <- as.data.frame(x)
  df[df$p <= 0.05, ]
})
o_nz <- names(Filter(function(dd) nrow(dd) > 0, sig_list))
if (length(o_nz)) {
  sig_df <- do.call(rbind, lapply(o_nz, function(name) cbind(Comparison = name, sig_list[[name]])))
  write.csv(sig_df, "SIMPER_significant_results.csv", row.names = FALSE)
}

set.seed(42)
seedNMDS <- metaMDS(comm = Hell_Seeds, distance = "bray", k = 3, trymax = 200)
stressplot(seedNMDS)
print(seedNMDS$stress)

cols <- c("orange", "green")
plot(seedNMDS)
orditorp(seedNMDS, display = "sites", label = mvb$sample, cex = 1, pch = "+", pcol = "grey", air = 0.01)
points(seedNMDS, cex = 1, pch = 16, col = cols[as.factor(mvb$BIRD)])
legend("topright", pch = 16, col = cols, legend = c("Gallinago gallinago", "Tringa glareola"), title = "Bird Species")
box()

###############################################################################
# SECTION B — Diversity indices and differences
###############################################################################

df_div <- read.csv2("merged_wider_javított2.csv", check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample" %in% names(df_div)) df_div$sample <- paste0("sample", seq_len(nrow(df_div)))
df_div <- df_div %>% relocate(BIRD, sample)

taxa_cols <- setdiff(names(df_div), c("BIRD","sample"))
X <- as.data.frame(df_div[, taxa_cols, drop = FALSE])
X_num <- as.data.frame(lapply(X, function(col) if (is.numeric(col)) col else suppressWarnings(as.numeric(col))))
stopifnot(sum(is.na(as.matrix(X_num))) == 0)

row_fun <- function(v) {
  s <- sum(v)
  if (s <= 0) return(c(shan = NA_real_, lambda = NA_real_, simpson = NA_real_, inv = NA_real_, rich = 0))
  p <- v / s
  lambda <- sum(p^2)
  c(shan = vegan::diversity(v, index = "shannon"), lambda = lambda, simpson = 1 - lambda, inv = 1 / lambda, rich = vegan::specnumber(v))
}
M <- t(apply(X_num, 1, row_fun)) %>% as.data.frame()

diversity_df <- df_div %>% bind_cols(M) %>% rename(shannon = shan, simpson_1mD = simpson, invsimpson = inv, lambda = lambda, richness = rich)

div2 <- diversity_df %>% mutate(BIRD = str_trim(BIRD)) %>% filter(BIRD %in% c("GALGAL","TRIGLA")) %>% mutate(BIRD = recode(BIRD, "GALGAL" = "Gallinago gallinago", "TRIGLA" = "Tringa glareola"))

w_shan <- wilcox.test(shannon ~ BIRD, data = div2)
w_s1mD <- wilcox.test(simpson_1mD ~ BIRD, data = div2)
w_inv  <- wilcox.test(invsimpson ~ BIRD, data = div2)
w_rich <- wilcox.test(richness ~ BIRD, data = div2, exact = FALSE)

effects <- bind_rows(
  rstatix::wilcox_effsize(div2, shannon ~ BIRD) %>% mutate(metric = "Shannon"),
  rstatix::wilcox_effsize(div2, simpson_1mD ~ BIRD) %>% mutate(metric = "Simpson (1 − D)"),
  rstatix::wilcox_effsize(div2, invsimpson ~ BIRD) %>% mutate(metric = "Inverse Simpson (1/D)"),
  rstatix::wilcox_effsize(div2, richness ~ BIRD) %>% mutate(metric = "Species richness")
) %>% select(metric, effsize, magnitude)

tests_tbl <- dplyr::bind_rows(
  rstatix::wilcox_test(div2, shannon ~ BIRD) %>% dplyr::mutate(metric = "Shannon"),
  rstatix::wilcox_test(div2, simpson_1mD ~ BIRD) %>% dplyr::mutate(metric = "Simpson (1 − D)"),
  rstatix::wilcox_test(div2, invsimpson ~ BIRD) %>% dplyr::mutate(metric = "Inverse Simpson (1/D)"),
  rstatix::wilcox_test(div2, richness ~ BIRD) %>% dplyr::mutate(metric = "Species richness")
) %>% rstatix::add_significance() %>% mutate(alternative = "two.sided") %>%
  select(metric, group1, group2, n1, n2, statistic, p, p.signif, alternative) %>%
  left_join(effects %>% rename(r = effsize), by = "metric") %>%
  mutate(W = round(statistic, 1), p = signif(p, 3), r = ifelse(is.na(r), NA, round(r, 3)), stars = dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")) %>%
  select(metric, group1, group2, n1, n2, W, p, stars, r, magnitude)
print(tests_tbl)

# ---- Diversity panel figure (Shannon, Simpson (1-D), Inverse Simpson, Richness) ----
p_lab <- function(p) paste0("p = ", format(p, digits = 3, scientific = TRUE))
mk <- function(var, ylab, pval) {
  ggplot(div2, aes(BIRD, .data[[var]], fill = BIRD)) +
    geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = 21) +
    geom_jitter(width = 0.12, alpha = 0.4) +
    labs(x = "Species", y = ylab) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.5, label = p_lab(pval)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
}

p1 <- mk("shannon",     "Shannon diversity",     w_shan$p.value)
p2 <- mk("simpson_1mD", "Simpson (1 − D)",       w_s1mD$p.value)
p3 <- mk("invsimpson",  "Inverse Simpson (1/D)", w_inv$p.value)
p4 <- mk("richness",    "Species richness",      w_rich$p.value)

final_fig <- (p1 | p2) / (p3 | p4)
final_fig

ggsave("diversity_indices_panel.png", final_fig, width = 9, height = 7, dpi = 300)

###############################################################################
# SECTION C — Within-species Bray–Curtis dissimilarity + permutation test
###############################################################################

comm <- mvb[, -c(1:2)]
grp  <- mvb$BIRD
dist_mat <- as.matrix(vegdist(comm^0.25, method = "bray"))

within_stats <- do.call(rbind, lapply(unique(grp), function(sp) {
  idx <- which(grp == sp)
  dvals <- dist_mat[idx, idx]
  dvals <- dvals[upper.tri(dvals)]
  data.frame(Species = sp, Mean_within_dissimilarity = mean(dvals), SD_within_dissimilarity = sd(dvals), N_pairs = length(dvals))
}))
print(within_stats)

set.seed(123)
nperm <- 999
perm_results <- do.call(rbind, lapply(unique(grp), function(sp) {
  idx <- which(grp == sp)
  obs <- mean(as.dist(dist_mat[idx, idx]))
  perm_means <- replicate(nperm, { rand_idx <- sample(1:nrow(comm), length(idx)); mean(as.dist(dist_mat[rand_idx, rand_idx])) })
  p_value <- mean(perm_means <= obs)
  data.frame(Species = sp, Obs_mean = obs, P_value = p_value)
}))
print(perm_results)

###############################################################################
# SECTION D — Relative abundance stacked bars (grouped by bird)
###############################################################################

df <- read.csv2("merged_wider_javított2.csv", stringsAsFactors = FALSE) %>%
  mutate(BIRD = recode(BIRD, "GALGAL" = "Gallinago gallinago", "TRIGLA" = "Tringa glareola", .default = BIRD)) %>%
  mutate(across(-c(BIRD, sample), ~ as.numeric(as.character(.))))

long_df <- df %>% pivot_longer(cols = -c(BIRD, sample), names_to = "Species", values_to = "Abundance") %>% filter(!is.na(Abundance) & Abundance > 0) %>% group_by(sample) %>% mutate(RelAbundance = Abundance / sum(Abundance)) %>% ungroup()

rare_species <- long_df %>% group_by(Species) %>% summarise(max_ab = max(RelAbundance), .groups = "drop") %>% filter(max_ab < 0.20) %>% pull(Species)

long_df <- long_df %>% mutate(Group = if_else(Species %in% rare_species, "<20% abundance", Species))
sample_levels <- long_df %>% distinct(BIRD, sample) %>% arrange(BIRD, sample) %>% pull(sample)
long_df <- long_df %>% mutate(sample = factor(sample, levels = sample_levels), BIRD = factor(BIRD))

ggplot(long_df, aes(sample, RelAbundance, fill = Group)) +
  geom_col(width = 0.9) +
  facet_grid(~ BIRD, scales = "free_x", space = "free_x") +
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Relative abundance", fill = "Taxa", title = "Microbiome Composition Grouped by Bird Species") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.background = element_blank(), strip.placement = "outside")

###############################################################################
# SECTION E — Core microbiome per bird
###############################################################################

df <- read.csv2("merged_wider_javított2.csv", check.names = FALSE) %>% mutate(BIRD = str_trim(BIRD), sample = str_trim(sample)) %>% mutate(BIRD = recode(BIRD, "GALGAL" = "Gallinago gallinago", "TRIGLA" = "Tringa glareola", .default = BIRD)) %>% mutate(across(-c(BIRD, sample), ~ suppressWarnings(as.numeric(as.character(.)))))

long_full <- df %>% pivot_longer(cols = -c(BIRD, sample), names_to = "Species", values_to = "Abundance") %>% mutate(Abundance = tidyr::replace_na(Abundance, 0), Species = str_trim(Species)) %>% group_by(sample) %>% mutate(total_abundance = sum(Abundance), RelAbundance = ifelse(total_abundance > 0, Abundance / total_abundance, 0)) %>% ungroup() %>% select(-total_abundance)

prevalence_threshold <- 0.80
abundance_threshold <- 0.001

prev_summary <- long_full %>% mutate(present = RelAbundance >= abundance_threshold) %>% group_by(BIRD, Species) %>% summarise(prevalence = mean(present), mean_abundance = mean(RelAbundance), .groups = "drop")

core_microbiome <- prev_summary %>% filter(prevalence >= prevalence_threshold) %>% arrange(BIRD, desc(prevalence), desc(mean_abundance))

plot_df <- core_microbiome %>% group_by(BIRD) %>% arrange(desc(mean_abundance), .by_group = TRUE) %>% mutate(Species = factor(Species, levels = rev(unique(Species)))) %>% ungroup()

gg_core <- ggplot(plot_df, aes(x = Species, y = mean_abundance, fill = prevalence)) + geom_col(width = 0.7) + geom_text(aes(label = scales::percent(prevalence, accuracy = 1)), hjust = -0.15, size = 3.3, color = "black") + coord_flip(clip = "off") + scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) + scale_fill_gradient(limits = c(prevalence_threshold, 1), low = "#9ecae1", high = "#3b5bdb", oob = scales::squish, name = "Prevalence") + facet_wrap(~ BIRD, ncol = 1, scales = "free_y") + labs(x = NULL, y = "Mean relative abundance") + theme_cowplot(font_size = 12) + theme(plot.margin = margin(5.5, 30, 5.5, 5.5), strip.text = element_text(face = "bold", size = 14))

gg_core

ggsave("core_microbiome_per_bird.pdf", gg_core, width = 180, height = 160, units = "mm")
ggsave("core_microbiome_per_bird.png", gg_core, width = 1600, height = 1400, units = "px", dpi = 300)
