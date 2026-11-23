#!/usr/bin/env Rscript
# Step 0b - Microbiome differential abundance (infected vs control), on the cascade genus table.
# Primary: ALDEx2 (CLR, Monte Carlo Dirichlet instances, Welch t, effect size).
# Sensitivity: ANCOM-BC2, and ALDEx2 on the strict cascade variant.
# PERMANOVA (Aitchison) is descriptive only: with 3 vs 3 the smallest possible p is 0.1.
# Run from the repository root with the integration_diatrea env:
#   Rscript scripts/step0b_microbiome_DA.R
suppressPackageStartupMessages({
  library(ALDEx2); library(ANCOMBC); library(phyloseq); library(vegan); library(ggplot2)
})
set.seed(20260930)

in_dir  <- "results/step0a_microbiome_cascade"
out_dir <- "results/step0b_microbiome_DA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
samples <- c("control_1", "control_2", "control_3", "infected_1", "infected_2", "infected_3")
group   <- factor(rep(c("control", "infected"), each = 3), levels = c("control", "infected"))
min_prevalence <- 3        # present in >= 3 of the 6 samples
mc_samples     <- 1000     # ALDEx2 Monte Carlo instances (>= 128 requested)
bh_cut         <- 0.1
effect_cut     <- 1
# genera removed by the original pipeline's fixed list; here they are kept and flagged
original_removal_list <- c("Cutibacterium", "Abiotrophia", "Akkermansia", "Alloprevotella",
  "Bacteroides", "Bifidobacterium", "Blautia", "Dialister", "Eggerthella", "Faecalibacterium",
  "Gardnerella", "Methanobrevibacter", "Neisseria", "Prevotella", "Ruminococcus", "Rothia",
  "Veillonella", "Dehalococcoides", "Anaerococcus", "Malassezia", "Candida", "Streptococcus")
likely_contaminant <- c("Cutibacterium", "Malassezia", "Prevotella")   # skin / human-associated

read_table <- function(f) {
  t <- read.delim(file.path(in_dir, f), check.names = FALSE)
  keep <- rowSums(t[, samples] > 0) >= min_prevalence
  list(all = t, kept = t[keep, ])
}

run_aldex <- function(t) {
  counts <- as.matrix(t[, samples]); rownames(counts) <- t$genus
  clr <- aldex.clr(counts, as.character(group), mc.samples = mc_samples, denom = "all",
                   verbose = FALSE)
  tt <- aldex.ttest(clr, paired.test = FALSE)          # Welch t (we.*) and Wilcoxon (wi.*)
  ef <- aldex.effect(clr, verbose = FALSE)
  res <- cbind(genus = rownames(tt), tt, ef)
  # ALDEx2 differences are (second level - first level) = infected - control, in log2 CLR units
  res
}

## ------------------------------------------------------------ primary table
prim <- read_table("abundance_table_genus_cascade.tsv")
t <- prim$kept
cat("Genera in cascade table:", nrow(prim$all), "| prevalent (>=", min_prevalence, "of 6):",
    nrow(t), "\n"); print(table(t$kingdom))

ald <- run_aldex(t)
stopifnot(identical(ald$genus, t$genus))
mb_diff <- ald[ald$genus == "Microbacterium", "diff.btw"]
cat("Sign of diff.btw for Microbacterium (expected > 0, higher in infected):", round(mb_diff, 2), "\n")
stopifnot(mb_diff > 0)

## ------------------------------------------------------------ ANCOM-BC2
counts <- as.matrix(t[, samples]); rownames(counts) <- t$genus
meta <- data.frame(group = group, row.names = samples)
ps <- phyloseq(otu_table(counts, taxa_are_rows = TRUE), sample_data(meta))
anc <- ancombc2(data = ps, fix_formula = "group", group = "group", p_adj_method = "BH", prv_cut = 0,
                lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, alpha = bh_cut,
                pseudo_sens = TRUE, n_cl = 4, verbose = FALSE)
ar <- anc$res
ar <- data.frame(genus = ar$taxon,
                 ancombc2_log2fc = ar$lfc_groupinfected / log(2),   # natural log -> log2
                 ancombc2_p = ar$p_groupinfected, ancombc2_q = ar$q_groupinfected,
                 ancombc2_diff = ar$diff_groupinfected,
                 ancombc2_passed_ss = ar$passed_ss_groupinfected)
sz <- anc$zero_ind
if (!is.null(sz)) {
  names(sz)[1] <- "genus"
  sz$ancombc2_structural_zero <- ifelse(sz[[2]] & !sz[[3]], "absent in control",
                                 ifelse(!sz[[2]] & sz[[3]], "absent in infected",
                                 ifelse(sz[[2]] & sz[[3]], "absent in both", "")))
  sz <- sz[, c("genus", "ancombc2_structural_zero")]
}

## ------------------------------------------------------------ strict cascade (sensitivity)
strict <- read_table("abundance_table_genus_cascade_strict.tsv")$kept
ald_s <- run_aldex(strict)
ald_s <- data.frame(genus = ald_s$genus, strict_log2_effect = ald_s$diff.btw,
                    strict_effect_size = ald_s$effect, strict_BH = ald_s$we.eBH)

## ------------------------------------------------------------ combine and apply the rule
da <- data.frame(
  genus = t$genus, kingdom = t$kingdom, phylum = t$phylum, tax = t$tax,
  feature_type = ifelse(grepl("_uncl\\.$", t$genus), "higher-rank group (_uncl.)", "genus"),
  mean_count_control = rowMeans(t[, samples[1:3]]), mean_count_infected = rowMeans(t[, samples[4:6]]),
  n_present_control = rowSums(t[, samples[1:3]] > 0), n_present_infected = rowSums(t[, samples[4:6]] > 0),
  log2_effect_infected_vs_control = ald$diff.btw,     # median CLR difference, log2
  aldex2_effect_size = ald$effect, aldex2_overlap = ald$overlap,
  p = ald$we.ep, BH = ald$we.eBH, wilcoxon_p = ald$wi.ep, wilcoxon_BH = ald$wi.eBH)
da <- merge(da, ar, by = "genus", all.x = TRUE)
if (!is.null(sz)) da <- merge(da, sz, by = "genus", all.x = TRUE)
da <- merge(da, ald_s, by = "genus", all.x = TRUE)

# ANCOM-BC2 gives no LFC for structural zeros (taxon absent from every sample of one group) and
# treats them as differentially abundant; they count as agreeing when the absent group matches
# the ALDEx2 direction.
struct_agrees <- !is.na(da$ancombc2_structural_zero) &
  ((da$ancombc2_structural_zero == "absent in control" & da$log2_effect_infected_vs_control > 0) |
   (da$ancombc2_structural_zero == "absent in infected" & da$log2_effect_infected_vs_control < 0))
consistent <- struct_agrees | (!is.na(da$ancombc2_log2fc) & !struct_agrees &
  sign(da$ancombc2_log2fc) == sign(da$log2_effect_infected_vs_control))
da$rule_BH   <- da$BH < bh_cut
da$rule_effect_ancom <- abs(da$aldex2_effect_size) > effect_cut & consistent
da$differential <- (da$rule_BH | da$rule_effect_ancom) & da$feature_type == "genus"
da$criterion <- ifelse(!da$rule_BH & !da$rule_effect_ancom, "",
                ifelse(da$rule_BH & da$rule_effect_ancom, "ALDEx2 BH<0.1 and |effect|>1 + ANCOM-BC2 direction",
                ifelse(da$rule_BH, "ALDEx2 BH<0.1", "|ALDEx2 effect|>1 + ANCOM-BC2 direction")))
da$direction <- ifelse(da$log2_effect_infected_vs_control > 0, "up_in_infected", "down_in_infected")
da$method <- "ALDEx2 (primary); ANCOM-BC2 (direction agreement)"
da$ancombc2_agrees <- consistent
# confidence tiers inside the differential set (the rule above defines "differential")
anc_sig <- !is.na(da$ancombc2_q) & da$ancombc2_q < bh_cut & da$ancombc2_passed_ss %in% TRUE
da$confidence_tier <- ifelse(!da$differential, "",
  ifelse(da$rule_BH & consistent, "1: ALDEx2 BH<0.1 + ANCOM-BC2 agrees",
  ifelse(da$rule_BH, "1b: ALDEx2 BH<0.1, ANCOM-BC2 disagrees",
  ifelse(anc_sig, "2: |effect|>1 + ANCOM-BC2 q<0.1 (passed sensitivity)",
  ifelse(struct_agrees, "3a: |effect|>1 + ANCOM-BC2 structural zero (one group only)",
         "3b: |effect|>1 + ANCOM-BC2 same direction, not significant")))))
da$strict_same_direction <- sign(da$strict_log2_effect) == sign(da$log2_effect_infected_vs_control)
da$contaminant_flag <- ifelse(da$genus %in% likely_contaminant, "likely skin/human-associated contaminant",
                       ifelse(da$genus %in% original_removal_list, "in original pipeline removal list", ""))
da <- da[order(!da$differential, -abs(da$aldex2_effect_size)), ]

first <- c("genus", "phylum", "kingdom", "log2_effect_infected_vs_control", "p", "BH", "method",
           "direction", "differential", "criterion", "confidence_tier")
da <- da[, c(first, setdiff(names(da), first))]
write.csv(da, file.path(out_dir, "microbiome_DA.csv"), row.names = FALSE)

## ------------------------------------------------------------ PERMANOVA (descriptive)
clr_mat <- t(apply(as.matrix(t[, samples]) + 0.5, 2, function(x) log(x) - mean(log(x))))
perm <- adonis2(dist(clr_mat) ~ group, permutations = permute::how(complete = TRUE, nperm = 720))
n_perm <- nrow(permute::allPerms(6)) + 1
cat("\n== PERMANOVA on Aitchison distance (descriptive)\n"); print(perm)
cat("All", n_perm, "permutations of 6 labels; each 3-vs-3 split is produced by 3!*3!*2 = 72 of them,",
    "so the smallest attainable p =", 72 / 720, "\n")
write.csv(as.data.frame(perm), file.path(out_dir, "permanova_aitchison.csv"))

## ------------------------------------------------------------ report
g <- da[da$feature_type == "genus", ]
cat("\n== Differential abundance (", nrow(g), "genera tested; _uncl. groups kept in CLR but not reported)\n")
cat("ALDEx2 BH <", bh_cut, ":", sum(g$rule_BH), "\n")
cat("|ALDEx2 effect| >", effect_cut, "with ANCOM-BC2 in the same direction:", sum(g$rule_effect_ancom), "\n")
cat("ANCOM-BC2 q <", bh_cut, ":", sum(g$ancombc2_q < bh_cut, na.rm = TRUE),
    "| of which passed sensitivity:", sum(g$ancombc2_q < bh_cut & g$ancombc2_passed_ss, na.rm = TRUE), "\n")
cat("Overlap ALDEx2 |effect|>1 vs ANCOM-BC2 q<0.1:\n")
print(table(ALDEx2_effect_gt1 = abs(g$aldex2_effect_size) > effect_cut,
            ANCOMBC2_q_lt_0.1 = g$ancombc2_q < bh_cut, useNA = "ifany"))
cat("Direction agreement ALDEx2 vs ANCOM-BC2 (all genera):", round(mean(g$ancombc2_agrees, na.rm = TRUE), 3), "\n")
d <- g[g$differential, ]
cat("\nDIFFERENTIAL genera:", nrow(d), "\n"); print(table(d$kingdom, d$direction))
print(table(d$confidence_tier, d$direction))
cat("Tiers 1-2 by kingdom:\n"); print(table(d$kingdom[grepl("^(1|2)", d$confidence_tier)]))
cat("Differential genera with the same direction in the strict cascade:",
    sum(d$strict_same_direction, na.rm = TRUE), "of", nrow(d), "(", sum(is.na(d$strict_same_direction)),
    "not prevalent in strict table )\n")
options(width = 220)
show <- c("genus", "kingdom", "log2_effect_infected_vs_control", "aldex2_effect_size", "BH",
          "ancombc2_log2fc", "ancombc2_q", "criterion", "contaminant_flag")
for (k in c("Bacteria", "Fungi", "Archaea")) {
  cat("\n--", k, "\n"); x <- d[d$kingdom == k, show]
  if (nrow(x)) print(format(x, digits = 3), row.names = FALSE) else cat("none\n")
}
cat("\nPositive control and flagged genera:\n")
print(format(g[g$genus %in% c("Fusarium", likely_contaminant, "Streptococcus", "Candida"), show],
             digits = 3), row.names = FALSE)

## ------------------------------------------------------------ figure
g$class <- ifelse(g$differential, g$direction, "not differential")
p <- ggplot(g, aes(aldex2_effect_size, log2_effect_infected_vs_control, colour = class, shape = kingdom)) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey60") +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = c(up_in_infected = "#B2182B", down_in_infected = "#2166AC",
                                 `not differential` = "grey70")) +
  ggrepel::geom_text_repel(data = subset(g, grepl("^(1|2)", confidence_tier) &
                             (abs(aldex2_effect_size) > 4 | kingdom == "Fungi")),
                           aes(label = genus), size = 2.6, max.overlaps = 30, show.legend = FALSE) +
  labs(caption = "Labels: confidence tiers 1-2 with |effect| > 4, and fungi",
       x = "ALDEx2 effect size (infected vs control)", y = "ALDEx2 log2 difference (CLR)",
       colour = NULL, shape = NULL) + theme_bw()
ggsave(file.path(out_dir, "microbiome_DA_effect_plot.png"), p, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "microbiome_DA_effect_plot.pdf"), p, width = 8, height = 6)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
