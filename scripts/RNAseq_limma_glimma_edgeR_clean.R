
library(tidyverse)
library(edgeR)
library(Glimma)
library(biomaRt)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(gplots)

library(clusterProfiler)
library(enrichplot)


# ============================================================
# LOAD COUNTS
# ============================================================
counts <- read.delim("data/DSAD-ApoE125_readcounts.txt", 
                     header = TRUE, 
                     skip = 1)
head(counts)

# Remove corrupted rows
bad_rows <- which(is.na(counts$Geneid) | !grepl("^ENSG", counts$Geneid) | 
                    rowSums(is.na(counts[, grepl("^R", colnames(counts))])) > 0)
length(bad_rows)
counts_clean <- counts[-bad_rows, ]

# Drop annotation columns, clean sample column names
counts_clean <- counts_clean %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)
colnames(counts_clean) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(counts_clean))

count_matrix <- as.matrix(counts_clean[, -1])
rownames(count_matrix) <- counts_clean$Geneid

# ============================================================
# LOAD AND CLEAN METADATA
# ============================================================
meta_raw <- read.csv("data/EMTAB14179_metadata.csv", check.names = FALSE)
colnames(meta_raw)

meta <- meta_raw %>%
  dplyr::select(
    sample_name = `Source Name`,
    age          = `Characteristics[age]`,
    brain_part   = `Characteristics[organism part]`,
    disease      = `Characteristics[disease]`,
    sex          = `Characteristics[sex]`,
    file_name    = `Comment[SUBMITTED_FILE_NAME]`,
    ENA_run      = `Comment[ENA_RUN]`,
    braak        = Braak,
    pmi          = PMI,
    adrc         = ADRC
  ) %>%
  mutate(
    sample_name = str_replace_all(sample_name, " ", "_"),
    brain_part  = str_replace_all(brain_part, " ", "_"),
    disease = case_when(
      disease == "normal" ~ "CTL",
      disease == "Down syndrome with Alzheimer's disease" ~ "DSAD",
      disease == "Alzheimer's disease" ~ "AD",
      TRUE ~ disease
    ),
    braak = as.factor(braak),
    adrc  = as.factor(adrc)
  )

meta <- meta %>%
  mutate(R_number = str_extract(file_name, "R[0-9]+")) %>%
  mutate(R_number = case_when(
    sample_name == "CTL_2_Cbl" ~ "R54_2",
    sample_name == "DSAD_21_Ctx" ~ "R115_2",
    sample_name == "sAD_9_Cbl" ~ "R61_2",
    TRUE ~ R_number
  )) %>%
  mutate(individual_id = str_remove(sample_name, "_Cbl$|_Ctx$")) %>%
  mutate(group_tissue = paste0(disease, "_", 
                               ifelse(brain_part == "prefrontal_cortex", "PC", "Cer"))) %>%
  dplyr::select(-file_name) %>%
  distinct()

# Check alignment with count matrix
setdiff(colnames(count_matrix), meta$R_number)
setdiff(meta$R_number, colnames(count_matrix))

# Reorder meta to match count_matrix column order
meta <- meta[match(colnames(count_matrix), meta$R_number), ]
all(colnames(count_matrix) == meta$R_number)  # should be TRUE

# ============================================================
# BUILD DGEList
# ============================================================
dge <- DGEList(counts = count_matrix)
dge$samples$group        <- meta$disease
dge$samples$brain_part   <- meta$brain_part
dge$samples$sex          <- meta$sex
dge$samples$age          <- meta$age
dge$samples$individual_id <- meta$individual_id
dge$samples$braak        <- meta$braak
dge$samples$pmi          <- meta$pmi
dge$samples$adrc         <- meta$adrc
dge$samples$group_tissue <- as.factor(meta$group_tissue)

head(dge$samples)
dim(dge)

# ============================================================
# GENE ANNOTATION (biomaRt)
# Cached to disk after the first successful query, since Ensembl's
# public BioMart service is occasionally unavailable and this
# annotation doesn't change between runs.
# ============================================================
annotation_cache <- "data/gene_annotation.rds"

if (file.exists(annotation_cache)) {
  genes_bm <- readRDS(annotation_cache)
} else {
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  genes_bm <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
                    filters = "ensembl_gene_id",
                    values = rownames(dge),
                    mart = mart)
  saveRDS(genes_bm, annotation_cache)
}

nrow(genes_bm)
nrow(dge)

genes_bm <- genes_bm[!duplicated(genes_bm$ensembl_gene_id), ]
genes_bm <- genes_bm[match(rownames(dge), genes_bm$ensembl_gene_id), ]

nrow(genes_bm) == nrow(dge)  # should be TRUE
sum(is.na(genes_bm$ensembl_gene_id))  # how many unmatched

dge$genes <- genes_bm
dge


# ============================================================
# FILTER LOW-EXPRESSION GENES
# ============================================================
keep.exprs <- filterByExpr(dge, group = dge$samples$group)
dge <- dge[keep.exprs, , keep.lib.sizes = FALSE]
dim(dge)

# ============================================================
# NORMALIZATION (with before/after comparison)
# ============================================================
dge_unnorm <- dge
dge_unnorm$samples$norm.factors <- 1

dge <- normLibSizes(dge, method = "TMM")
dge$samples$norm.factors

par(mfrow = c(1, 2))
boxplot(cpm(dge_unnorm, log = TRUE), las = 2, main = "Before TMM")
boxplot(cpm(dge, log = TRUE), las = 2, main = "After TMM")
par(mfrow = c(1, 1))

# ============================================================
# CONVERT KEY VARIABLES TO FACTORS
# ============================================================
dge$samples$group      <- as.factor(dge$samples$group)
dge$samples$brain_part <- as.factor(dge$samples$brain_part)
dge$samples$sex        <- as.factor(dge$samples$sex)
# braak, adrc, group_tissue already factors - confirmed by str()
# pmi, age stay numeric - correct as-is

str(dge$samples)

# ============================================================
# FULL-DATASET MDS (sanity check: brain region dominates)
# ============================================================
lcpm <- cpm(dge, log = TRUE)

col.group <- dge$samples$group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

col.brain_part <- dge$samples$brain_part
levels(col.brain_part) <- brewer.pal(max(nlevels(col.brain_part), 3), "Set2")
col.brain_part <- as.character(col.brain_part)

par(mfrow = c(1, 2))
plotMDS(lcpm, labels = dge$samples$group, col = col.group)
title("A. Sample disease")

plotMDS(lcpm, labels = dge$samples$brain_part, col = col.brain_part, dim.plot = c(1,2))
title("B. Brain part")
par(mfrow = c(1, 1))

# Confirm: brain_part splits dim1 cleanly, group does not
mds <- plotMDS(lcpm, plot = FALSE)
table(dge$samples$group, sign(mds$x))
table(dge$samples$brain_part, sign(mds$x))

col.group_tissue <- dge$samples$group_tissue
levels(col.group_tissue) <- brewer.pal(max(nlevels(col.group_tissue), 3), "Set2")
col.group_tissue <- as.character(col.group_tissue)
plotMDS(lcpm, labels = dge$samples$group_tissue, col = col.group_tissue, dim.plot = c(1,2))

# ============================================================
# CORTEX SUBSET: AD vs CTL
# (DSAD excluded — different site/age, addressed separately
# in cross-tissue interaction analysis)
# ============================================================
dge_pc <- dge[, dge$samples$brain_part == "prefrontal_cortex" &
                dge$samples$group %in% c("AD", "CTL")]
dge_pc$samples <- droplevels(dge_pc$samples)

keep.exprs.pc <- filterByExpr(dge_pc, group = dge_pc$samples$group)
dge_pc <- dge_pc[keep.exprs.pc, , keep.lib.sizes = FALSE]
dge_pc <- normLibSizes(dge_pc, method = "TMM")

dim(dge_pc)
table(dge_pc$samples$group)

# ============================================================
# CONFOUNDER CHECKS: AD vs CTL cortex
# (both groups USC — site not a concern here)
# ============================================================
table(dge_pc$samples$group, dge_pc$samples$sex)

aggregate(age ~ group, data = dge_pc$samples,
          FUN = function(x) c(mean=mean(x), range=range(x)))

aggregate(pmi ~ group, data = dge_pc$samples,
          FUN = function(x) c(mean=mean(x, na.rm=TRUE), range=range(x, na.rm=TRUE)))

table(dge_pc$samples$group, dge_pc$samples$braak)
# Note: braak and group are nearly collinear (CTL=0-1, AD=3-6)
# — confirms braak should NOT be added as covariate alongside group

# Does PMI correlate with expression variation?
lcpm_pc_tmp <- cpm(dge_pc, log = TRUE)
mds_tmp <- plotMDS(lcpm_pc_tmp, plot = FALSE)
cor(dge_pc$samples$pmi, mds_tmp$x, use = "complete.obs")



# ============================================================
# CORTEX MDS: AD vs CTL
# ============================================================
par(mfrow = c(1,1), mar = c(5,4,4,8), xpd = FALSE)
lcpm_pc <- cpm(dge_pc, log = TRUE)

group        <- dge_pc$samples$group
sex          <- dge_pc$samples$sex
group_colors <- brewer.pal(nlevels(group), "Set1")
col.group    <- group_colors[group]
sex_shapes   <- c(16, 17)
pch.sex      <- sex_shapes[sex]

mds <- plotMDS(lcpm_pc, plot = TRUE, pch = pch.sex, col = col.group, cex = 1.5)
title("Cortex: AD vs CTL")

text(mds$x, mds$y, labels = dge_pc$samples$age, pos = 2, cex = 0.7)

par(xpd = TRUE)
xlim <- range(mds$x); ylim <- range(mds$y)
legend(x = xlim[2]+0.3, y = ylim[2],
       legend = levels(group), col = group_colors, pch = 16,
       title = "Disease", bty = "n")
legend(x = xlim[2]+0.3, y = ylim[2]-1.2,
       legend = levels(sex), pch = sex_shapes,
       title = "Sex", bty = "n", col = "black")
par(mfrow = c(1,1), mar = c(5,4,4,2), xpd = FALSE)

# ============================================================
# CORTEX PRIMARY DE: AD vs CTL (age and sex adjusted)
# ============================================================

design_pc <- model.matrix(~0 + group + age + sex, data = dge_pc$samples)
colnames(design_pc) <- gsub("group", "", colnames(design_pc))

contr_pc <- makeContrasts(AD_vs_CTL = AD - CTL, levels = colnames(design_pc))

par(mfrow = c(1,2))
v_pc <- voom(dge_pc, design_pc, plot = TRUE)
par(mfrow = c(1,1))

vfit_pc <- lmFit(v_pc, design_pc)
vfit_pc <- contrasts.fit(vfit_pc, contrasts = contr_pc)
efit_pc <- eBayes(vfit_pc)
plotSA(efit_pc, main = "Cortex: Mean-variance trend")

summary(decideTests(efit_pc))

# Top 10 — report tendency even if no FDR-significant genes
topTable(efit_pc, coef = "AD_vs_CTL", n = 10, sort.by = "P") %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, P.Value, adj.P.Val, description)

# ============================================================
# SENSITIVITY CHECK: PMI-adjusted (AD/CTL both USC — clean adjustment)
# ============================================================
dge_pc_pmi <- dge_pc[, !is.na(dge_pc$samples$pmi)]
dge_pc_pmi$samples <- droplevels(dge_pc_pmi$samples)

keep <- filterByExpr(dge_pc_pmi, group = dge_pc_pmi$samples$group)
dge_pc_pmi <- dge_pc_pmi[keep, , keep.lib.sizes = FALSE]
dge_pc_pmi <- normLibSizes(dge_pc_pmi, method = "TMM")

design_pc_pmi <- model.matrix(~0 + group + age + sex + pmi, data = dge_pc_pmi$samples)
colnames(design_pc_pmi) <- gsub("group", "", colnames(design_pc_pmi))

v_pc_pmi  <- voom(dge_pc_pmi, design_pc_pmi, plot = FALSE)
vfit_pmi  <- lmFit(v_pc_pmi, design_pc_pmi)
vfit_pmi  <- contrasts.fit(vfit_pmi,
                           makeContrasts(AD_vs_CTL = AD - CTL, levels = colnames(design_pc_pmi)))
efit_pmi  <- eBayes(vfit_pmi)

summary(decideTests(efit_pmi))
# Expected: still null — confirms result robust to PMI adjustment

# ============================================================
# EXPLORATORY: SEX-STRATIFIED AD vs CTL (cortex)
# ============================================================

# --- FEMALES ---
dge_pc_f <- dge_pc[, dge_pc$samples$sex == "female"]
dge_pc_f$samples <- droplevels(dge_pc_f$samples)

keep <- filterByExpr(dge_pc_f, group = dge_pc_f$samples$group)
dge_pc_f <- dge_pc_f[keep, , keep.lib.sizes = FALSE]
dge_pc_f <- normLibSizes(dge_pc_f, method = "TMM")

dim(dge_pc_f)
table(dge_pc_f$samples$group)

design_f <- model.matrix(~0 + group + age, data = dge_pc_f$samples)
colnames(design_f) <- gsub("group", "", colnames(design_f))

v_f    <- voom(dge_pc_f, design_f, plot = FALSE)
vfit_f <- lmFit(v_f, design_f)
vfit_f <- contrasts.fit(vfit_f,
                        makeContrasts(AD_vs_CTL = AD - CTL, levels = colnames(design_f)))
efit_f <- eBayes(vfit_f)

summary(decideTests(efit_f))

top_f <- topTable(efit_f, coef = "AD_vs_CTL", n = 10, sort.by = "P") %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, P.Value, adj.P.Val)
top_f

# --- MALES ---
dge_pc_m <- dge_pc[, dge_pc$samples$sex == "male"]
dge_pc_m$samples <- droplevels(dge_pc_m$samples)

keep <- filterByExpr(dge_pc_m, group = dge_pc_m$samples$group)
dge_pc_m <- dge_pc_m[keep, , keep.lib.sizes = FALSE]
dge_pc_m <- normLibSizes(dge_pc_m, method = "TMM")

dim(dge_pc_m)
table(dge_pc_m$samples$group)

design_m <- model.matrix(~0 + group + age, data = dge_pc_m$samples)
colnames(design_m) <- gsub("group", "", colnames(design_m))

v_m    <- voom(dge_pc_m, design_m, plot = FALSE)
vfit_m <- lmFit(v_m, design_m)
vfit_m <- contrasts.fit(vfit_m,
                        makeContrasts(AD_vs_CTL = AD - CTL, levels = colnames(design_m)))
efit_m <- eBayes(vfit_m)

summary(decideTests(efit_m))

top_m <- topTable(efit_m, coef = "AD_vs_CTL", n = 10, sort.by = "P") %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, P.Value, adj.P.Val)
top_m

# --- COMPARE TOP HITS BETWEEN SEXES ---
cat("Genes in female top 10 but not male top 10:\n")
setdiff(top_f$external_gene_name, top_m$external_gene_name)

cat("Genes in male top 10 but not female top 10:\n")
setdiff(top_m$external_gene_name, top_f$external_gene_name)

cat("Genes shared between both top 10s:\n")
intersect(top_f$external_gene_name, top_m$external_gene_name)



# ============================================================
# CEREBELLUM SUBSET: AD vs CTL
# ============================================================
dge_cer <- dge[, dge$samples$brain_part == "cerebellum" & 
                 dge$samples$group %in% c("AD", "CTL")]
dge_cer$samples <- droplevels(dge_cer$samples)

keep.exprs.cer <- filterByExpr(dge_cer, group = dge_cer$samples$group)
dge_cer <- dge_cer[keep.exprs.cer, , keep.lib.sizes = FALSE]
dge_cer <- normLibSizes(dge_cer, method = "TMM")

dim(dge_cer)
table(dge_cer$samples$group)

# Confounder checks
table(dge_cer$samples$group, dge_cer$samples$sex)
aggregate(age ~ group, data = dge_cer$samples,
          FUN = function(x) c(mean=mean(x), range=range(x)))
aggregate(pmi ~ group, data = dge_cer$samples,
          FUN = function(x) c(mean=mean(x, na.rm=TRUE), range=range(x, na.rm=TRUE)))
table(dge_cer$samples$group, dge_cer$samples$braak)

# ============================================================
# CEREBELLUM MDS: AD vs CTL
# ============================================================
par(mfrow = c(1,1), mar = c(5,4,4,8), xpd = FALSE)
lcpm_cer <- cpm(dge_cer, log = TRUE)

group        <- dge_cer$samples$group
sex          <- dge_cer$samples$sex
group_colors <- brewer.pal(nlevels(group), "Set1")
col.group    <- group_colors[group]
sex_shapes   <- c(16, 17)
pch.sex      <- sex_shapes[sex]

mds_cer <- plotMDS(lcpm_cer, plot = TRUE, pch = pch.sex, col = col.group, cex = 1.5)
title("Cerebellum: AD vs CTL")
text(mds_cer$x, mds_cer$y, labels = dge_cer$samples$age, pos = 2, cex = 0.7)

par(xpd = TRUE)
xlim <- range(mds_cer$x); ylim <- range(mds_cer$y)
legend(x = xlim[2]+0.3, y = ylim[2],
       legend = levels(group), col = group_colors, pch = 16,
       title = "Disease", bty = "n")
legend(x = xlim[2]+0.3, y = ylim[2]-2.5,
       legend = levels(sex), pch = sex_shapes,
       title = "Sex", bty = "n", col = "black")
par(mfrow = c(1,1), mar = c(5,4,4,2), xpd = FALSE)

# ============================================================
# CEREBELLUM PRIMARY DE: AD vs CTL (age + sex adjusted)
# ============================================================
design_cer <- model.matrix(~0 + group + age + sex, data = dge_cer$samples)
colnames(design_cer) <- gsub("group", "", colnames(design_cer))

contr_cer <- makeContrasts(AD_vs_CTL = AD - CTL, levels = colnames(design_cer))

par(mfrow = c(1,2))
v_cer <- voom(dge_cer, design_cer, plot = TRUE)
par(mfrow = c(1,1))

vfit_cer <- lmFit(v_cer, design_cer)
vfit_cer <- contrasts.fit(vfit_cer, contrasts = contr_cer)
efit_cer <- eBayes(vfit_cer)
plotSA(efit_cer, main = "Cerebellum: Mean-variance trend")

summary(decideTests(efit_cer))

top_cer <- topTable(efit_cer, coef = "AD_vs_CTL", n = 10, sort.by = "P") %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, P.Value, adj.P.Val, description)
top_cer

# ============================================================
# COMPARE CORTEX vs CEREBELLUM TOP HITS: AD vs CTL
# (do the same genes show tendency in both regions?)
# ============================================================
top_pc <- topTable(efit_pc, coef = "AD_vs_CTL", n = 10, sort.by = "P") %>%
  dplyr::select(external_gene_name, logFC, P.Value, adj.P.Val)

cat("Genes in cortex top 10 but not cerebellum:\n")
setdiff(top_pc$external_gene_name, top_cer$external_gene_name)

cat("Genes in cerebellum top 10 but not cortex:\n")
setdiff(top_cer$external_gene_name, top_pc$external_gene_name)

cat("Genes shared across both regions:\n")
intersect(top_pc$external_gene_name, top_cer$external_gene_name)




# ============================================================
# BRAAK GRADIENT ANALYSIS: CORTEX (CTL + AD, age + sex adjusted)
# ============================================================

# Use CTL + AD cortex samples with non-NA Braak
dge_braak <- dge[, dge$samples$brain_part == "prefrontal_cortex" &
                   dge$samples$group %in% c("AD", "CTL") &
                   !is.na(dge$samples$braak)]
dge_braak$samples <- droplevels(dge_braak$samples)

# Convert braak to numeric for continuous regression
dge_braak$samples$braak_num <- as.numeric(as.character(dge_braak$samples$braak))

dim(dge_braak)
table(dge_braak$samples$braak_num)

keep <- filterByExpr(dge_braak, group = dge_braak$samples$group)
dge_braak <- dge_braak[keep, , keep.lib.sizes = FALSE]
dge_braak <- normLibSizes(dge_braak, method = "TMM")

# Design: Braak as continuous covariate, adjusted for age and sex
# No group term — Braak replaces it as the continuous severity measure
design_braak <- model.matrix(~braak_num + age + sex, data = dge_braak$samples)
design_braak

par(mfrow = c(1,2))
v_braak <- voom(dge_braak, design_braak, plot = TRUE)
par(mfrow = c(1,1))

vfit_braak <- lmFit(v_braak, design_braak)
efit_braak_pc <- eBayes(vfit_braak)

summary(decideTests(efit_braak_pc))

# Top genes associated with Braak severity
topTable(efit_braak_pc, coef = "braak_num", n = 10, sort.by = "P") %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, P.Value, adj.P.Val, description)


# ============================================================
# BRAAK GRADIENT ANALYSIS: CEREBELLUM (CTL + AD, age + sex adjusted)
# ============================================================

# Use CTL + AD cerebellum samples with non-NA Braak
dge_braak <- dge[, dge$samples$brain_part == "cerebellum" &
                   dge$samples$group %in% c("AD", "CTL") &
                   !is.na(dge$samples$braak)]
dge_braak$samples <- droplevels(dge_braak$samples)

# Convert braak to numeric for continuous regression
dge_braak$samples$braak_num <- as.numeric(as.character(dge_braak$samples$braak))

dim(dge_braak)
table(dge_braak$samples$braak_num)

keep <- filterByExpr(dge_braak, group = dge_braak$samples$group)
dge_braak <- dge_braak[keep, , keep.lib.sizes = FALSE]
dge_braak <- normLibSizes(dge_braak, method = "TMM")

# Design: Braak as continuous covariate, adjusted for age and sex
# No group term — Braak replaces it as the continuous severity measure
design_braak <- model.matrix(~braak_num + age + sex, data = dge_braak$samples)
design_braak

par(mfrow = c(1,2))
v_braak <- voom(dge_braak, design_braak, plot = TRUE)
par(mfrow = c(1,1))

vfit_braak <- lmFit(v_braak, design_braak)
efit_braak_cer <- eBayes(vfit_braak)

summary(decideTests(efit_braak_cer))

top_braak_cer <- topTable(efit_braak_cer, coef = "braak_num", n = 10, sort.by = "P") %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, P.Value, adj.P.Val, description)
top_braak_cer

# Compare cortex vs cerebellum Braak-associated genes
cat("Braak-associated in cortex but not cerebellum:\n")
top_braak_pc <- topTable(efit_braak_pc, coef = "braak_num", n = 10, sort.by = "P")
setdiff(top_braak_pc$external_gene_name, top_braak_cer$external_gene_name)

cat("Braak-associated in cerebellum but not cortex:\n")
setdiff(top_braak_cer$external_gene_name, top_braak_pc$external_gene_name)

cat("Shared Braak-associated genes across both regions:\n")
intersect(top_braak_pc$external_gene_name, top_braak_cer$external_gene_name)



# ============================================================
# CROSS-TISSUE INTERACTION: CORTEX vs CEREBELLUM
# (within-individual, all three groups: CTL, AD, DSAD)
# Using cerebellum as internal reference per individual
# Site confound less problematic here: within-person Ctx-Cer 
# comparison is site-internal for each group
# ============================================================

# Use full dataset (both tissues, all groups)


table(dge$samples$group_tissue)
dim(dge)

# ============================================================
# DESIGN: group_tissue (6 levels, no intercept)
# blocking on individual_id for paired Ctx/Cer per person
# ============================================================
design_int <- model.matrix(~0 + group_tissue, data = dge$samples)
colnames(design_int) <- levels(dge$samples$group_tissue)
# levels: AD_Cer, AD_PC, DSAD_Cer, DSAD_PC, CTL_Cer, CTL_PC

par(mfrow = c(1,2))
v_full <- voom(dge, design_int, plot = TRUE)
par(mfrow = c(1,1))

# Account for paired Cbl/Ctx from same individual
corfit <- duplicateCorrelation(v_full, design_int, 
                               block = dge$samples$individual_id)
corfit$consensus  # within-person correlation — expect positive ~0.3-0.6

vfit_full <- lmFit(v_full, design_int,
                   block = dge$samples$individual_id,
                   correlation = corfit$consensus)

# ============================================================
# CONTRASTS: difference-of-differences
# "Does the Ctx-Cer gap differ by disease group vs CTL?"
# ============================================================
contr_int <- makeContrasts(
  AD_region_effect   = (AD_PC - AD_Cer) - (CTL_PC - CTL_Cer),
  DSAD_region_effect = (DSAD_PC - DSAD_Cer) - (CTL_PC - CTL_Cer),
  levels = colnames(design_int)
)

vfit_full <- contrasts.fit(vfit_full, contrasts = contr_int)
efit_full <- eBayes(vfit_full)

summary(decideTests(efit_full))

# ============================================================
# TOP GENES PER CONTRAST
# ============================================================
cat("--- AD region effect (Ctx-Cer gap differs from CTL?) ---\n")
topTable(efit_full, coef = "AD_region_effect", n = 10, sort.by = "P") %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, P.Value, adj.P.Val, description)

cat("--- DSAD region effect (Ctx-Cer gap differs from CTL?) ---\n")
top_dsad_region <- topTable(efit_full, coef = "DSAD_region_effect", n = 20, sort.by = "P") %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, P.Value, adj.P.Val, description)
top_dsad_region

# ============================================================
# CHECK DSAD REGION EFFECT FOR AD-RELEVANT PATHWAYS
# autophagy/trafficking (from earlier KIF5A, TSC2, PEX6 finding)
# ============================================================
ad_pathway_genes <- c("KIF5A", "TSC2", "PEX6", "RAPGEF4", "MCTP2",
                      "TREM2", "AIF1", "C1QA", "C1QB", "TYROBP",
                      "VGF", "SCG2", "SCG3", "APP", "MAPT", "APOE")

dsad_results %>%
  filter(external_gene_name %in% ad_pathway_genes) %>%
  dplyr::select(external_gene_name, gene_biotype, logFC, adj.P.Val, description)

# ============================================================
# COMPARE AD vs DSAD REGION EFFECTS
# (which genes show disrupted regional pattern in BOTH vs CTL?)
# ============================================================
top_ad_region <- topTable(efit_full, coef = "AD_region_effect", 
                          n = 20, sort.by = "P")$external_gene_name
top_dsad_region_names <- top_dsad_region$external_gene_name

cat("Shared regional disruption in both AD and DSAD:\n")
intersect(top_ad_region, top_dsad_region_names)

cat("Regional disruption specific to DSAD only:\n")
setdiff(top_dsad_region_names, top_ad_region)

# ============================================================
# VISUALIZATIONS: DSAD REGION EFFECT
# (only analysis with substantial FDR-significant signal)
# ============================================================

# ---- 1. Mean-difference plot (static) ----
dt_full <- decideTests(efit_full)

plotMD(efit_full, column = "DSAD_region_effect", 
       status = dt_full[, "DSAD_region_effect"],
       main = "DSAD: Cortex-Cerebellum regional disruption",
       xlim = c(-10, 10))

# ---- 2. Heatmap: top DSAD region effect genes ----

# Get top 50 DSAD region effect genes by adjusted p-value
dsad_results <- topTable(efit_full, coef = "DSAD_region_effect", 
                            n = Inf, sort.by = "P")
top50_genes <- dsad_results$ensembl_gene_id[1:50]
top50_genes <- top50_genes[!is.na(top50_genes)]

# Get log-CPM for these genes across all samples
lcpm_full <- cpm(dge, log = TRUE)
i <- which(rownames(lcpm_full) %in% top50_genes)

# Color by group_tissue for column labels
group_tissue_labels <- dge$samples$group_tissue
mycol <- colorpanel(1000, "blue", "white", "red")

# Gene labels: use symbol where available, Ensembl ID otherwise
gene_labels <- ifelse(
  is.na(dge$genes$external_gene_name[i]) | 
    dge$genes$external_gene_name[i] == "",
  rownames(lcpm_full)[i],
  dge$genes$external_gene_name[i]
)

heatmap.2(lcpm_full[i, ],
          scale = "row",
          labRow = gene_labels,
          labCol = group_tissue_labels,
          col = mycol,
          trace = "none",
          density.info = "none",
          margin = c(10, 8),
          lhei = c(2, 10),
          dendrogram = "column",
          main = "Top 50 DSAD region-effect genes")

# ---- 3. Volcano plot: DSAD region effect ----

dsad_volcano <- dsad_results %>%
  mutate(
    sig = adj.P.Val < 0.05,
    label = ifelse(sig & abs(logFC) > 2, external_gene_name, "")
  )

ggplot(dsad_volcano, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("grey60", "red3"),
                     labels = c("Not significant", "FDR < 0.05")) +
  geom_text(aes(label = label), size = 2.5, hjust = -0.1,
            check_overlap = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "DSAD: Cortex vs Cerebellum regional disruption",
       x = "Log2 Fold Change (PC - Cer difference vs CTL)",
       y = "-log10(P-value)",
       color = "") +
  theme_bw() +
  theme(legend.position = "bottom")

dsad_results %>%
  filter(grepl("^S100", external_gene_name)) %>%
  dplyr::select(external_gene_name, logFC, adj.P.Val, description)



# ============================================================
# PATHWAY ENRICHMENT: DSAD REGION EFFECT
# Using clusterProfiler for GO and KEGG enrichment
# ============================================================


# ============================================================
# PREPARE GENE LISTS
# ============================================================


# Significant genes (FDR < 0.05): 883 down + 662 up
dsad_sig <- dsad_results %>% filter(adj.P.Val < 0.05)
cat("Total significant DSAD region effect genes:", nrow(dsad_sig), "\n")
cat("Up:", sum(dsad_sig$logFC > 0), "Down:", sum(dsad_sig$logFC < 0), "\n")

# For enrichment: need Entrez IDs (clusterProfiler works best with these)
# Map Ensembl → Entrez via org.Hs.eg.db
library(org.Hs.eg.db)

sig_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = dsad_sig$ensembl_gene_id[!is.na(dsad_sig$ensembl_gene_id)],
                                    columns = "ENTREZID",
                                    keytype = "ENSEMBL") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

cat("Genes with Entrez ID mapping:", nrow(sig_entrez), "\n")

# Universe: all tested genes (for background in enrichment)
universe_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = dsad_results$ensembl_gene_id[!is.na(dsad_results$ensembl_gene_id)],
                                         columns = "ENTREZID",
                                         keytype = "ENSEMBL") %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

cat("Universe size:", length(universe_entrez), "\n")


# ============================================================
# BUILD RANKED GENE LIST FOR GSEA
# ============================================================

# Derive from dsad_results — add rank score
dsad_ranked <- dsad_results %>%
  filter(!is.na(ensembl_gene_id)) %>%
  mutate(rank_score = sign(logFC) * -log10(P.Value))

# Map Ensembl → Entrez
ranked_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = dsad_ranked$ensembl_gene_id,
                                       columns = "ENTREZID",
                                       keytype = "ENSEMBL") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# Join rank scores with Entrez IDs
ranked_genes <- dsad_ranked %>%
  left_join(ranked_entrez, by = c("ensembl_gene_id" = "ENSEMBL")) %>%
  filter(!is.na(ENTREZID)) %>%
  arrange(desc(rank_score))

# Remove duplicate Entrez IDs — keep highest absolute rank score per ID
ranked_genes_dedup <- ranked_genes %>%
  group_by(ENTREZID) %>%
  slice_max(abs(rank_score), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(rank_score))

# Named numeric vector — required format for gseGO/gseKEGG
gene_rank_vec <- setNames(ranked_genes_dedup$rank_score, 
                          ranked_genes_dedup$ENTREZID)

# Verify
cat("Total genes in ranked list:", length(gene_rank_vec), "\n")
cat("Unique Entrez IDs:", length(unique(names(gene_rank_vec))), "\n")
head(gene_rank_vec)
tail(gene_rank_vec)

# Now run GSEA KEGG
gsea_kegg <- gseKEGG(
  geneList      = gene_rank_vec,
  organism      = "hsa",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE,
  seed          = TRUE
)

cat("GSEA KEGG significant pathways:", nrow(as.data.frame(gsea_kegg)), "\n")

as.data.frame(gsea_kegg) %>%
  dplyr::select(Description, NES, p.adjust, setSize) %>%
  arrange(p.adjust) %>%
  head(20)

# Also retry gseGO with deduplicated vector
gsea_go <- gseGO(
  geneList      = gene_rank_vec,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE,
  seed          = TRUE
)

cat("GSEA GO BP significant terms:", nrow(as.data.frame(gsea_go)), "\n")

as.data.frame(gsea_go) %>%
  dplyr::select(Description, NES, p.adjust, setSize) %>%
  arrange(p.adjust) %>%
  head(20)


# ============================================================
# VISUALIZE KEGG GSEA (17 significant pathways)
# ============================================================
dotplot(gsea_kegg, showCategory = 17,
        title = "DSAD region effect: GSEA KEGG pathways") +
  theme(axis.text.y = element_text(size = 8))

# GSEA enrichment plots for key pathways
gseaplot2(gsea_kegg,
          geneSetID = "hsa05010",
          title = "Alzheimer Disease pathway: DSAD region effect")

gseaplot2(gsea_kegg,
          geneSetID = "hsa04142",
          title = "Lysosome biogenesis: DSAD region effect")

gseaplot2(gsea_kegg,
          geneSetID = "hsa05208",
          title = "ROS pathway: DSAD region effect")

# ============================================================
# VISUALIZE GO GSEA (99 significant terms)
# ============================================================
# Top 20 by adjusted p-value
as.data.frame(gsea_go) %>%
  dplyr::select(Description, NES, p.adjust, setSize) %>%
  arrange(p.adjust) %>%
  head(20)

dotplot(gsea_go, showCategory = 20,
        title = "DSAD region effect: GSEA GO Biological Process") +
  theme(axis.text.y = element_text(size = 7))

# Check for AD-relevant GO terms specifically
as.data.frame(gsea_go) %>%
  filter(grepl("autophagy|lysosom|oxidative|mitochond|axon|vesicle|proteostasis|ubiquitin",
               Description, ignore.case = TRUE)) %>%
  dplyr::select(Description, NES, p.adjust, setSize) %>%
  arrange(p.adjust)





# Check which individuals have both tissues
tissue_per_individual <- meta %>%
  group_by(individual_id) %>%
  summarise(
    has_ctx = any(brain_part == "prefrontal_cortex"),
    has_cer = any(brain_part == "cerebellum"),
    n_samples = n()
  )

# Complete pairs (both tissues)
complete_pairs <- tissue_per_individual %>% filter(has_ctx & has_cer)
cat("Complete pairs (both tissues):", nrow(complete_pairs), "\n")

# Missing one tissue
cat("Cortex only (no cerebellum):\n")
tissue_per_individual %>% filter(has_ctx & !has_cer) %>% pull(individual_id)

cat("Cerebellum only (no cortex):\n")
tissue_per_individual %>% filter(has_cer & !has_ctx) %>% pull(individual_id)