# =====================================================================
# ECC-MR mechanistic follow-up: coloc + PheWAS for decoded pleiotropic loci
#
# Prerequisites: ECCMR_main_analysis.R has been run, so res_bmi etc. and
#   iop_combine / glaucoma_combine / bmi_combine / t2d_combine
#   are still in the environment (columns: chr, pos, rsid, ref, alt, beta, se, pval, maf, N)
#
# Dependencies: coloc, ieugwasr (PheWAS requires an OpenGWAS token, see Section 3),
#       ggplot2, reshape2 (optional)
# =====================================================================

## ------------------------------------------------------------------
## 1. Coloc: exposure-outcome colocalisation at each decoded locus
##    Question: is the pleiotropy "one causal variant acting on both traits" (high PP.H4,
##    genuine horizontal pleiotropy) or "two distinct variants bundled by LD" (high PP.H3, LD artefact)?
## ------------------------------------------------------------------
library(coloc)

# Extract one locus region from the full sumstats and align alleles (reuses the harmonise logic)
.region_pair <- function(exp_df, out_df, chr, pos_lo, pos_hi) {
  e <- exp_df[exp_df$chr == chr & exp_df$pos >= pos_lo & exp_df$pos <= pos_hi, ]
  o <- out_df[out_df$chr == chr & out_df$pos >= pos_lo & out_df$pos <= pos_hi, ]
  m <- merge(e, o, by = "rsid", suffixes = c("_X", "_Y"))
  same <- m$ref_X == m$ref_Y & m$alt_X == m$alt_Y
  flip <- m$ref_X == m$alt_Y & m$alt_X == m$ref_Y
  m <- m[same | flip, ]; flip <- flip[same | flip]
  m$beta_Y[flip] <- -m$beta_Y[flip]
  m$maf_Y[flip]  <- 1 - m$maf_Y[flip]
  pal <- paste0(m$ref_X, m$alt_X) %in% c("AT","TA","CG","GC")
  m <- m[!(pal & m$maf_X > 0.42 & m$maf_X < 0.58), ]
  m <- m[!duplicated(m$rsid), ]
  m
}

# Run coloc.abf for all decoded loci
#   res    : object returned by run_eccmr_pair()
#   s_out  : case proportion for a binary outcome (set NULL to skip, e.g. FinnGen H7_GLAUCOMA)
#   win    : window by which locus boundaries are extended
run_coloc_decoded <- function(res, exp_df, out_df, win = 250000,
                              outcome_type = c("cc", "quant"),
                              s_out = NULL) {
  outcome_type <- match.arg(outcome_type)
  loci <- unique(res$pleiotropic_snps$locus)
  cat("Colocalisation analysis:", length(loci), "decoded loci\n")

  rows <- lapply(loci, function(L) {
    sub <- res$pleiotropic_snps[res$pleiotropic_snps$locus == L, ]
    chr <- sub$chr_X[1]
    lo  <- min(sub$pos_X) - win; hi <- max(sub$pos_X) + win
    m <- .region_pair(exp_df, out_df, chr, lo, hi)
    if (nrow(m) < 50) {
      cat(sprintf("  locus %d (chr%d): only %d SNPs in region, skipped\n", L, chr, nrow(m)))
      return(NULL)
    }
    ## coloc requires MAF strictly in (0,1): filter missing/degenerate values before building datasets
    maf1 <- pmin(m$maf_X, 1 - m$maf_X)
    maf2 <- pmin(m$maf_Y, 1 - m$maf_Y)
    ok <- is.finite(maf1) & maf1 > 0 & maf1 < 1 &
           is.finite(maf2) & maf2 > 0 & maf2 < 1
    if (sum(!ok) > 0)
      message(sprintf("  locus %s: dropped %d SNPs with invalid MAF", L, sum(!ok)))
    m <- m[ok, ]; maf1 <- maf1[ok]; maf2 <- maf2[ok]
    d1 <- list(beta = m$beta_X, varbeta = m$se_X^2, N = m$N_X,
               snp = m$rsid, type = "quant", MAF = maf1)
    d2 <- list(beta = m$beta_Y, varbeta = m$se_Y^2, N = m$N_Y,
               snp = m$rsid, type = outcome_type, MAF = maf2)
    if (outcome_type == "cc" && !is.null(s_out)) d2$s <- s_out
    cl <- tryCatch(coloc.abf(d1, d2),
                   error = function(e) { message("  locus ", L, " failed: ",
                                                 conditionMessage(e)); NULL })
    if (is.null(cl)) return(NULL)
    pp <- cl$summary
    data.frame(
      locus = L, chr = chr, pos_lo = lo + win, pos_hi = hi - win,
      n_snps_region = nrow(m),
      n_decoded = nrow(sub),
      PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"],
      PP.H2 = pp["PP.H2.abf"], PP.H3 = pp["PP.H3.abf"],
      PP.H4 = pp["PP.H4.abf"],
      verdict = ifelse(pp["PP.H4.abf"] > 0.8, "shared causal variant (genuine pleiotropy)",
                ifelse(pp["PP.H3.abf"] > 0.8, "LD artefact (distinct variants)",
                       "inconclusive")),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

## ------------------------------------------------------------------
## 2. PheWAS: which other traits do the decoded SNPs hit -> infer the biological channel of pleiotropy
##
##   An OpenGWAS token is required:
##     1. log in with GitHub at https://api.opengwas.io/profile/ to generate a token
##     2. add one line to ~/.Renviron:  OPENGWAS_JWT=your_token
##   136 SNPs are queried one by one with a 1s rate limit; about 3-5 minutes in total
## ------------------------------------------------------------------
library(ieugwasr)

run_phewas_decoded <- function(res, pval = 5e-8, batch = c("ieu-b","ebi-a","finn-b"),
                               sleep_s = 1.0) {
  if (Sys.getenv("OPENGWAS_JWT") == "")
    stop("OPENGWAS_JWT not found; please configure the token in ~/.Renviron first")
  snps <- res$pleiotropic_snps$rsid
  cat("PheWAS query for", length(snps), "decoded SNPs ...\n")
  out <- lapply(snps, function(rs) {
    r <- tryCatch(phewas(variants = rs, pval = pval, batch = batch),
                  error = function(e) NULL)
    Sys.sleep(sleep_s)
    r
  })
  out <- do.call(rbind, out)
  cat("Retrieved", nrow(out), "SNP-trait associations (P <", pval, ")\n")
  out
}

## universal biological-domain classifier for PheWAS traits
classify_trait <- function(tr) {
  tr <- tolower(tr)
  ifelse(grepl("glucose|a1c|glycated|glycaemi|glycemi|hba1c|diabetes|insulin|medication|sugar", tr),
         "Glycaemic",
         ifelse(grepl("mass index|weight|body fat|waist|hip|metabolic rate|lean mass|height|mineral density|adiposity|obesity|overweight|anthropometric|fat percentage|trunk fat|arm fat|leg fat|whole body", tr),
                "Adiposity",
                ifelse(grepl("cholesterol|cholesteryl|\\bldl\\b|\\bhdl\\b|triglycer|lipid|lipoprotein|apolipoprotein|apob|apoa|statin|remnant|phospholipid", tr),
                       "Lipid",
                       ifelse(grepl("blood pressure|systolic|diastolic|hypertension|pulse|arterial pressure", tr),
                              "Blood pressure",
                              ifelse(grepl("platelet|white blood|red blood|lymphocyte|neutrophil|monocyte|eosinophil|basophil|haematocrit|hematocrit|cell count|corpuscular|haemoglobin|hemoglobin|reticulocyte", tr),
                                     "Haematologic",
                                     ifelse(grepl("coronary|heart|cardiovascular|angina|myocardial|vascular|stroke|artery|arterial disease|atherosclero|ischaemi|ischemi|atrial", tr),
                                            "Cardiovascular", "Other"))))))
}
## ------------------------------------------------------------------
## 3. PheWAS heatmap (supplementary-material quality; size adapts to content)
##    rows = traits (top by frequency), columns = decoded SNPs (rsid, grouped by locus),
##    colour = -log10(P); a PNG is always written; PDF uses the cairo device (failure tolerated)
##    note: the panoramic figure is >178 mm wide, better suited to the supplement; plot representative loci for the main text
## ------------------------------------------------------------------
## ------------------------------------------------------------------
## 2b. automatic nearest-gene annotation of loci (biomaRt / Ensembl, requires internet)
##     returns a locus_names vector, passed directly to plot_phewas_heatmap()
## ------------------------------------------------------------------
annotate_loci_genes <- function(res, known = c("453" = "TCF7L2", "638" = "FTO",
                                               "331" = "JAZF1")) {
  if (!requireNamespace("biomaRt", quietly = TRUE))
    stop("Please install biomaRt first: BiocManager::install('biomaRt')")
  loc <- res$pleiotropic_snps
  id_col  <- intersect(c("rsid", "snp", "SNP"), colnames(loc))[1]
  chr_col <- intersect(c("chr", "chr_X", "CHR"), colnames(loc))[1]
  pos_col <- intersect(c("pos", "pos_X", "POS"), colnames(loc))[1]
  loc <- loc[, c(id_col, chr_col, pos_col, "locus")]
  colnames(loc) <- c("snp", "chr", "pos", "locus")

  ## representative SNP position per locus: the one with the strongest signal
  rep <- do.call(rbind, lapply(split(loc, loc$locus), function(x) x[1, ]))
  mart <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                              mirror = "www")
  genes <- vapply(seq_len(nrow(rep)), function(i) {
    g <- tryCatch(
      biomaRt::getBM(attributes = "hgnc_symbol",
                     filters = c("chromosome_name", "start", "end"),
                     values = list(rep$chr[i], rep$pos[i] - 5e4, rep$pos[i] + 5e4),
                     mart = mart),
      error = function(e) data.frame(hgnc_symbol = character(0)))
    g <- g$hgnc_symbol[g$hgnc_symbol != ""]
    if (length(g) == 0) "" else g[1]
  }, character(1))
  out <- setNames(genes, as.character(rep$locus))
  ## manually verified names take precedence (standard usage in the literature)
  out[names(known)] <- known
  message("Locus gene annotation:\n", paste(sprintf("  L%s -> %s", names(out), out),
                                   collapse = "\n"))
  out
}

## nearest/overlapping genes for the 21 decoded BMI->T2D loci (GRCh37 coordinates, Ensembl annotation,
## manually verified; KCNJ11/GIPR/SLC2A2/PPARG/TCF7L2/FTO/GNPDA2/VEGFA/BCL2 are
## classic metabolic genes confirmed in the literature)
locus_genes_bmi <- c("108" = "COBLL1", "133" = "PPARG",  "138" = "RARB",
                     "167" = "ADCY5",  "174" = "ARHGEF26", "179" = "SLC2A2",
                     "202" = "GNPDA2", "212" = "UNC5C",  "258" = "JADE2",
                     "280" = "ZNF322", "288" = "LRFN2",  "291" = "VEGFA",
                     "314" = "ARG1",   "331" = "JAZF1",  "453" = "TCF7L2",
                     "457" = "ZRANB1", "470" = "KCNJ11", "490" = "SMCO4",
                     "638" = "FTO",    "691" = "BCL2",   "707" = "GIPR")

plot_phewas_heatmap <- function(ph, res, top_traits = 25,
                                out_prefix = "Fig_phewas_heatmap",
                                locus_names = locus_genes_bmi) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible())
  library(ggplot2)

  ## column names adapt automatically (chr / chr_X, pos / pos_X, rsid / snp)
  loc <- res$pleiotropic_snps
  id_col  <- intersect(c("rsid", "snp", "SNP"), colnames(loc))[1]
  chr_col <- intersect(c("chr", "chr_X", "CHR"), colnames(loc))[1]
  pos_col <- intersect(c("pos", "pos_X", "POS"), colnames(loc))[1]
  if (any(is.na(c(id_col, chr_col, pos_col))) || !"locus" %in% colnames(loc))
    stop("res$pleiotropic_snps lacks rsid/chr/pos/locus columns; actual columns: ",
         paste(colnames(loc), collapse = ", "))
  loc <- loc[, c(id_col, chr_col, pos_col, "locus")]
  colnames(loc) <- c("snp", "chr", "pos", "locus")

  ph$trait <- iconv(ph$trait, from = "UTF-8", to = "ASCII//TRANSLIT")  # strip special characters
  ## phewas results carry chr/position columns; drop them first to avoid suffixed names after merging
  ph_sub <- ph[, setdiff(colnames(ph), c("chr", "pos", "position", "locus"))]

  ## ---- trait deduplication: repeated definitions of the same biological trait across datasets are merged,
  ##      keeping the dataset with the largest sample size N ----
  n_before <- length(unique(ph_sub$trait))
  ## normalise: lowercase + strip source suffixes such as "(UKB data field 21001)"
  canon <- tolower(trimws(sub("\\s*\\(.*$", "", ph_sub$trait)))
  ## unify common synonymous spellings (capitalisation / spelling variants)
  canon <- gsub("^body mass index.*", "Body mass index", canon)
  canon <- gsub("^type 2 diabetes.*|^diabetes$", "Type 2 diabetes", canon)
  canon <- gsub("^(glycated )?h(ae)?emoglobin a1c.*", "HbA1c", canon)
  canon <- gsub("^glucose.*", "Glucose", canon)
  canon <- gsub("^sex hormone-binding globulin.*", "SHBG", canon)
  canon <- gsub("^hip circumference.*", "Hip circumference", canon)
  canon <- gsub("^waist circumference.*", "Waist circumference", canon)
  canon <- gsub("^weight$", "Weight", canon)
  canon <- gsub("^body fat percentage.*", "Body fat percentage", canon)
  canon <- gsub("^height$", "Height", canon)
  ph_sub$trait <- paste0(toupper(substr(canon, 1, 1)),
                         substr(canon, 2, nchar(canon)))   # capitalise the first letter for display
  ## keep only the largest-N dataset per normalised trait
  if (all(c("id", "n") %in% colnames(ph_sub))) {
    N_by <- tapply(ph_sub$n, list(ph_sub$trait, ph_sub$id),
                   function(x) max(x, na.rm = TRUE))
    best_id <- apply(N_by, 1, function(x) colnames(N_by)[which.max(x)])
    keep_row <- mapply(function(tr, id) isTRUE(best_id[tr] == id),
                       ph_sub$trait, ph_sub$id)
    ph_sub <- ph_sub[keep_row, ]
    message("Trait deduplication: ", n_before, " -> ", length(unique(ph_sub$trait)),
            " traits (largest-N dataset kept per trait)")
  }

  d <- merge(ph_sub, loc, by.x = "rsid", by.y = "snp")
  d <- d[order(d$chr, d$pos), ]

  cnt  <- table(d$trait)
  keep <- names(sort(cnt, decreasing = TRUE))[seq_len(min(top_traits, length(cnt)))]
  d    <- d[d$trait %in% keep, ]
  d$log10p <- pmin(-log10(d$p), 40)              # cap to prevent extreme values from swallowing the colour scale

  ## ---- Y axis: traits grouped by biological domain (universal classifier covering glycaemic/adiposity/lipid/blood-pressure/haematologic/cardiovascular) ----
  classify <- function(tr) {
    tr <- tolower(tr)
    ifelse(grepl("glucose|a1c|glycated|glycaemi|glycemi|hba1c|diabetes|insulin|medication|sugar", tr),
           "Glycaemic",
    ifelse(grepl("mass index|weight|body fat|waist|hip|metabolic rate|lean mass|height|mineral density|adiposity|obesity|overweight|anthropometric|fat percentage|trunk fat|arm fat|leg fat|whole body", tr),
           "Adiposity",
    ifelse(grepl("cholesterol|cholesteryl|\bldl\b|\bhdl\b|triglycer|lipid|lipoprotein|apolipoprotein|apob|apoa|statin|remnant|phospholipid", tr),
           "Lipid",
    ifelse(grepl("blood pressure|systolic|diastolic|hypertension|pulse|arterial pressure", tr),
           "Blood pressure",
    ifelse(grepl("platelet|white blood|red blood|lymphocyte|neutrophil|monocyte|eosinophil|basophil|haematocrit|hematocrit|cell count|corpuscular|haemoglobin|hemoglobin|reticulocyte", tr),
           "Haematologic",
    ifelse(grepl("coronary|heart|cardiovascular|angina|myocardial|vascular|stroke|artery|arterial disease|atherosclero|ischaemi|ischemi|atrial", tr),
           "Cardiovascular", "Other"))))))
  }
  grp_levels <- c("Glycaemic", "Adiposity", "Lipid", "Blood pressure",
                  "Haematologic", "Cardiovascular", "Other")
  d$group <- factor(classify(d$trait), levels = grp_levels)
  d$group <- droplevels(d$group)
  ## within each group, sort by number of associated SNPs (descending)
  ord_t <- names(sort(tapply(d$rsid, d$trait, function(x) length(unique(x))),
                      decreasing = TRUE))
  d$trait <- factor(d$trait, levels = rev(ord_t))

  ## ---- X axis: keep only rsids on the main panel; gene names go on a top annotation track (patchwork) ----
  d$rsid <- factor(d$rsid, levels = unique(d$rsid))
  loc_map <- tapply(d$locus, d$rsid, function(x) x[1])
  boundaries <- which(diff(as.integer(factor(loc_map[levels(d$rsid)],
                                             levels = unique(loc_map)))) != 0) + 0.5
  ## without a gene table (locus_names empty) the annotation track is left blank; no interruption
  has_genes <- length(locus_names) > 0
  if (has_genes) {
    gene_map <- setNames(locus_names[as.character(loc_map)], levels(d$rsid))
    gene_map[is.na(gene_map)] <- ""
    ## each named locus is labelled once, over its middle SNP
    for (g in unique(gene_map[gene_map != ""])) {
      idx <- which(gene_map == g)
      keep <- idx[ceiling(length(idx) / 2)]
      gene_map[setdiff(idx, keep)] <- ""
    }
  } else {
    gene_map <- setNames(rep("", length(levels(d$rsid))), levels(d$rsid))
  }

  ## adaptive canvas: widen each column when few SNPs; font size scales with canvas width (fs),
  ## so text does not become too small when a wide figure is shrunk to the page
  n_snp <- nlevels(d$rsid)
  per_snp <- if (n_snp <= 10) 10 else if (n_snp <= 20) 8 else if (n_snp <= 40) 5 else 3.4
  w_mm <- max(80 + n_snp * per_snp, 170)
  fs <- max(1, w_mm / 350)

  g <- ggplot(d, aes(rsid, trait, fill = log10p)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = boundaries, colour = "grey40", linewidth = 0.3) +
    scale_fill_gradientn(colours = c("grey95", "#F4A582", "#B2182B"),
                         values = scales::rescale(c(8, 20, 40)),
                         limits = c(8, 40), oob = scales::squish,
                         name = expression(-log[10](italic(P)))) +
    facet_grid(rows = vars(group), scales = "free_y", space = "free_y",
               switch = "y") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 8 * fs) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4.5 * fs),
          axis.text.y = element_text(size = 7 * fs),
          panel.grid   = element_blank(),
          panel.spacing = unit(2 * fs, "mm"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 8 * fs),
          legend.position = "right",
          plot.margin = margin(5, 5, 5, 5))

  ## top gene annotation track: known loci are labelled (red, bold), others blank
  g_anno <- ggplot(data.frame(rsid = factor(names(gene_map), levels = names(gene_map)),
                              gene = gene_map, y = 1),
                   aes(rsid, y, label = gene)) +
    geom_text(angle = 90, hjust = 0, size = 2.6 * fs, colour = "#B2182B",
              fontface = "bold") +
    geom_vline(xintercept = boundaries, colour = "grey40", linewidth = 0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 2.2))) +
    theme_void(base_size = 8 * fs) +
    theme(plot.margin = margin(5, 5, 0, 5))

  grp_cnt <- table(d$group)
  n_grp <- as.integer(grp_cnt)
  ## height: row height adapts and never exceeds the width (keep landscape/square, avoid overly tall figures)
  per_row <- if (sum(n_grp) <= 30) 5 else 4
  h_mm <- min(65 + sum(n_grp) * per_row, w_mm)

  if (has_genes && requireNamespace("patchwork", quietly = TRUE)) {
    g_all <- g_anno + g + patchwork::plot_layout(heights = c(1, 6), guides = "collect") &
      theme(legend.position = "right")
  } else {
    message("Hint: install patchwork for the top gene annotation track (install.packages('patchwork'))")
    g_all <- g
  }

  ggsave(paste0(out_prefix, ".png"), g_all, width = w_mm, height = h_mm,
         units = "mm", dpi = 300, limitsize = FALSE)
  tryCatch(
    ggsave(paste0(out_prefix, ".pdf"), g_all, width = w_mm, height = h_mm,
           units = "mm", device = cairo_pdf, limitsize = FALSE),
    error = function(e) message("PDF export failed (font issue); PNG is available: ",
                                conditionMessage(e))
  )
  message("Saved: ", out_prefix, ".png (", w_mm, " x ", h_mm, " mm)\n",
          "trait groups: ", paste(names(grp_cnt), grp_cnt, sep = "=", collapse = " / "))
  invisible(g_all)
}

## ==================================================================
## main workflow (BMI -> T2D as example; switch to res_iop + the corresponding data frames likewise)
## ==================================================================

## --- coloc ---
## s_out = outcome case proportion (fill in the actual proportion if T2D comes from UKB/DIAM; NULL if unsure)
coloc_bmi <- run_coloc_decoded(res_bmi, bmi_combine, t2d_combine,
                               outcome_type = "cc", s_out = NULL)
write.csv(coloc_bmi, "coloc_BMI_T2D_decoded_loci.csv", row.names = FALSE)
print(coloc_bmi[, c("locus","chr","n_decoded","PP.H3","PP.H4","verdict")])

## --- PheWAS ---
ph_bmi <- run_phewas_decoded(res_bmi, pval = 5e-8)
ph_bmi$domain <- classify_trait(ph_bmi$trait) 
write.csv(ph_bmi, "phewas_BMI_T2D_decoded_snps.csv", row.names = FALSE)

## --- locus gene annotation ---
## the built-in locus_genes_bmi is used by default (all 21 loci manually annotated; no internet needed)
plot_phewas_heatmap(ph_bmi, res_bmi, top_traits = 25,
                    out_prefix = "Fig_phewas_heatmap")
## to use biomaRt online annotation instead (requires access to Ensembl):
## locus_genes <- annotate_loci_genes(res_bmi)
## plot_phewas_heatmap(ph_bmi, res_bmi, top_traits = 25,
##                     locus_names = locus_genes,
##                     out_prefix = "Fig_phewas_heatmap")

## ------------------------------------------------------------------
## LDL -> CAD mechanistic analysis (Willer 2013 GLGC -> GCST003116)
## binary outcome; CARDIoGRAMplusC4D 2015 case proportion = 60801/184306, approx. 0.33
## ------------------------------------------------------------------
## --- coloc ---
coloc_ldl <- run_coloc_decoded(res_ldl, ldl_combine, cad_combine,
                               outcome_type = "cc", s_out = 0.33)
write.csv(coloc_ldl, "coloc_LDL_CAD_decoded_loci.csv", row.names = FALSE)
print(coloc_ldl[, c("locus","chr","n_decoded","PP.H3","PP.H4","verdict")])

## --- PheWAS ---
ph_ldl <- run_phewas_decoded(res_ldl, pval = 5e-8)
ph_ldl$domain <- classify_trait(ph_ldl$trait) 
write.csv(ph_ldl, "phewas_LDL_CAD_decoded_snps.csv", row.names = FALSE)

## --- heatmap ---
## gene annotation for the 5 decoded LDL loci (Ensembl GRCh37 offline annotation, 2026-09):
##   6=EVI5, 54=SH2B3 (pleiotropy hotspot), 19=FN1, 67=TOMM40/APOE, 26=CSNK1G3
locus_genes_ldl <- c("6" = "EVI5", "54" = "SH2B3", "19" = "FN1",
                     "67" = "TOMM40/APOE", "26" = "CSNK1G3")
plot_phewas_heatmap(ph_ldl, res_ldl, top_traits = 25,
                    locus_names = locus_genes_ldl,
                    out_prefix = "Fig_phewas_heatmap_LDL")

# =====================================================================
# Interpretation notes
# 1. coloc PP.H4 > 0.8: shared causal variant -> genuine horizontal pleiotropy (TCF7L2 expected);
#    PP.H3 > 0.8: distinct variants bundled by LD -> downgrade to "linked signal" in the manuscript.
# 2. In PheWAS, look for whether the traits hit by decoded SNPs form a "channel beyond the exposure"
#    (e.g. TCF7L2 -> fasting insulin/HOMA-B supports an insulin-secretion channel).
# =====================================================================




# ============================================================
# Coloc regional plots for decoded pleiotropic loci (LocusZoom-style two-track figure)
# Companion to the mechanism section of this script:
#   first run .region_pair() and run_coloc_decoded() to obtain
#   the coloc results table, then draw each locus with plot_coloc_region()
#
# Figure layout (two tracks sharing the genomic x axis):
#   top track:    -log10(P) of the exposure
#   bottom track: -log10(P) of the outcome, y axis reversed ("back to back")
#   red diamonds = pleiotropic SNPs decoded by ECC-MR
#   grey points  = other SNPs in the region
#   red dashed line = lead SNP (smallest exposure P among decoded SNPs), labelled in the top track
#   PP.H3 / PP.H4 of the locus are annotated top right; the title appears on the top track only
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## ------------------------------------------------------------
## annotate_loci_genes_rest(): lightweight alternative to biomaRt (needs only httr/jsonlite,
## uses the Ensembl REST API, unaffected by Bioconductor mirrors; GRCh37 via the grch37 subdomain)
##   res   : object returned by run_eccmr_pair()
##   known : manually verified c("locus" = "gene") entries; they take precedence over automatic annotation
##   returns: named vector, locus id -> gene symbol, passed directly to locus_genes
## ------------------------------------------------------------
annotate_loci_genes_rest <- function(res, known = c(), build = c("grch37", "grch38")) {
  build <- match.arg(build)
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE))
    stop("Please run install.packages(c('httr','jsonlite')) first")
  base <- if (build == "grch37") "https://grch37.rest.ensembl.org"
  else "https://rest.ensembl.org"
  
  loc <- res$pleiotropic_snps
  id_col  <- intersect(c("rsid", "snp", "SNP"), names(loc))[1]
  chr_col <- intersect(c("chr", "chr_X", "CHR"), names(loc))[1]
  pos_col <- intersect(c("pos", "pos_X", "bp", "POS"), names(loc))[1]
  loc <- loc[, c(id_col, chr_col, pos_col, "locus")]
  colnames(loc) <- c("snp", "chr", "pos", "locus")
  
  ## representative position per locus: the largest |alpha| if available, otherwise the first
  rep <- do.call(rbind, lapply(split(loc, loc$locus), function(x) x[1, ]))
  genes <- vapply(seq_len(nrow(rep)), function(i) {
    url <- sprintf("%s/overlap/region/human/%s:%d-%d?feature=gene;content-type=application/json",
                   base, rep$chr[i], rep$pos[i] - 50000, rep$pos[i] + 50000)
    g <- tryCatch({
      r <- httr::GET(url)
      if (httr::status_code(r) != 200) return("")
      js <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))
      pc <- js$external_name[js$biotype == "protein_coding" & js$external_name != ""]
      if (length(pc) == 0) "" else pc[1]
    }, error = function(e) "")
    Sys.sleep(0.4)   # polite delay for the REST rate limit
    g
  }, character(1))
  out <- setNames(genes, as.character(rep$locus))
  out[names(known)] <- known
  message("Locus gene annotation:\n", paste(sprintf("  L%s -> %s", names(out), out), collapse = "\n"))
  out
}

## ------------------------------------------------------------
## plot_coloc_region(): draw the two-track regional plot for one locus
##
##   res        : object returned by run_eccmr_pair() (with $pleiotropic_snps)
##   exp_df     : full exposure sumstats (same as used for coloc)
##   out_df     : full outcome sumstats
##   locus      : locus id (a value in res$pleiotropic_snps$locus)
##   win        : window beyond locus boundaries (default +/-250 kb, as in coloc)
##   coloc_row  : the row of the coloc results table for this locus (with PP.H3, PP.H4),
##                pass NULL to skip the annotation
##   locus_genes: named vector, locus id -> gene symbol, e.g. c("54" = "SH2B3");
##                can be generated with annotate_loci_genes_rest(),
##                NULL keeps the numeric locus id in the title
##   exp_name   : exposure label for axis titles, e.g. "BMI" / "LDL-C"
##   out_name   : outcome label, e.g. "T2D" / "CAD"
##   out_pdf    : output pdf path; NULL returns the plot object without saving
## ------------------------------------------------------------
plot_coloc_region <- function(res, exp_df, out_df, locus, win = 250000,
                              coloc_row = NULL, locus_genes = NULL,
                              exp_name = "Exposure", out_name = "Outcome",
                              out_pdf = NULL) {
  
  ## --- 1. decoded SNPs and region boundaries for this locus (column names adapt: rsid/SNP/snp) ---
  sub <- res$pleiotropic_snps[res$pleiotropic_snps$locus == locus, ]
  if (nrow(sub) == 0) stop("locus ", locus, " not found in res$pleiotropic_snps")
  pick <- function(df, cands) { n <- intersect(cands, names(df));
  if (length(n) == 0) stop("column not found: ", paste(cands, collapse = "/"),
                           "; actual columns: ", paste(names(df), collapse = ", ")); df[[n[1]]] }
  rs_col <- pick(sub, c("rsid", "SNP", "snp"))
  chr    <- pick(sub, c("chr_X", "chr"))[1]
  pos_v  <- pick(sub, c("pos_X", "pos", "bp"))
  lo  <- min(pos_v) - win
  hi  <- max(pos_v) + win
  decoded_rs <- rs_col
  if (!any(decoded_rs %in% exp_df$rsid))
    message("Warning: none of the decoded SNP rsids match the exposure data; check rsid version consistency")
  
  ## --- 2. regional data (reuses .region_pair above) ---
  m <- .region_pair(exp_df, out_df, chr, lo, hi)
  if (nrow(m) < 10) { message("locus ", locus, ": too few SNPs in region, skipped"); return(invisible(NULL)) }
  
  d <- data.frame(
    rsid    = m$rsid,
    pos     = m$pos_X,
    lp_X    = -log10(p_X <- 2 * pnorm(-abs(m$beta_X / m$se_X))),
    lp_Y    = -log10(2 * pnorm(-abs(m$beta_Y / m$se_Y))),
    decoded = m$rsid %in% decoded_rs
  )
  
  ## --- 3. PP.H3 / PP.H4 annotation text ---
  lab <- if (!is.null(coloc_row)) {
    sprintf("PP.H3 = %.3f\nPP.H4 = %.3f", coloc_row$PP.H3, coloc_row$PP.H4)
  } else ""
  
  gname <- if (!is.null(locus_genes)) locus_genes[as.character(locus)] else NA_character_
  ttl <- if (!is.na(gname) && nzchar(gname))
    sprintf("%s locus (chr%d: %.2f-%.2f Mb)", gname, chr, lo/1e6, hi/1e6)
  else
    sprintf("chr%d: %.2f-%.2f Mb (locus %d)", chr, lo/1e6, hi/1e6, locus)
  
  ## --- 3b. lead SNP: smallest exposure P among decoded SNPs; dashed vertical line in both tracks, labelled in the top track ---
  lead <- if (any(d$decoded)) d[d$decoded, ][which.max(d$lp_X[d$decoded]), ] else NULL
  
  ## --- 4. the two tracks ---
  mk_track <- function(yvar, ylab, flip = FALSE, show_title = FALSE,
                       label_lead = FALSE) {
    g <- ggplot(d, aes(pos, .data[[yvar]]))
    if (!is.null(lead))
      g <- g + geom_vline(xintercept = lead$pos, linetype = "dashed",
                          colour = "#C0392B", linewidth = 0.4, alpha = 0.8)
    g <- g +
      geom_point(colour = "grey70", size = 1.1) +
      geom_point(data = d[d$decoded, ], colour = "#C0392B",
                 shape = 23, fill = "#E74C3C", size = 2.6, stroke = 0.6) +
      annotate("label", x = Inf, y = Inf, label = lab,
               hjust = 1.05, vjust = 1.2, size = 3.2,
               colour = "grey20", fill = "white", alpha = 0.85) +
      labs(y = ylab, x = NULL, title = if (show_title) ttl else NULL) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 10),
            panel.grid.minor = element_blank())
    if (label_lead && !is.null(lead))
      g <- g + annotate("label", x = lead$pos, y = Inf, label = lead$rsid,
                        vjust = 1.2, hjust = -0.1, size = 3.2,
                        colour = "#C0392B", fill = "white", fontface = "bold")
    if (flip) g <- g + scale_y_reverse()   # reverse the bottom track -> back-to-back
    g
  }
  
  g_x <- mk_track("lp_X", bquote(-log[10](italic(P)) ~ .(exp_name)),
                  show_title = TRUE, label_lead = TRUE)
  g_y <- mk_track("lp_Y", bquote(-log[10](italic(P)) ~ .(out_name)), flip = TRUE) +
    scale_x_continuous(labels = function(v) sprintf("%.2f", v/1e6)) +
    labs(x = sprintf("Position on chr%d (Mb)", chr)) +
    theme(axis.text.x = element_text(size = 8))
  
  g_all <- g_x / g_y + patchwork::plot_layout(heights = c(1, 1))
  
  ## --- 5. save (fall back to grDevices::pdf if cairo_pdf fails) ---
  if (!is.null(out_pdf)) {
    ok <- tryCatch({ ggplot2::ggsave(out_pdf, g_all, width = 180, height = 130,
                                     units = "mm", device = cairo_pdf); TRUE },
                   error = function(e) FALSE)
    if (!ok) ggplot2::ggsave(out_pdf, g_all, width = 180, height = 130,
                             units = "mm", device = grDevices::pdf)
    message("Saved: ", out_pdf)
  }
  invisible(g_all)
}
BiocManager::install("biomaRt")

## ============================================================
## Example 1: LDL-C -> CAD
## ============================================================
##
## ## prerequisites: the mechanism script has been run, so the environment has:
## ##   res_ldl               (run_eccmr_pair result)
## ##   ldl_combine, cad_combine  (full sumstats)
## ##   .region_pair()           (function from this script)
## ##   annotate_loci_genes_rest() (gene annotation via Ensembl REST)
## coloc_ldl <- read.csv("coloc_LDL_CAD_decoded_loci.csv")
##
## ## locus id -> gene name, two options:
## ## (a) via Ensembl REST, no biomaRt needed (recommended):
 genes_ldl <- annotate_loci_genes_rest(res_ldl,
                known = c("54" = "SH2B3", "67" = "APOE",
                          "19" = "FN1",  "26" = "CSNK1G3"))
## ## (b) if biomaRt is installed, use annotate_loci_genes() from section 2b:
## ## genes_ldl <- annotate_loci_genes(res_ldl, known = c(... same as above ...))
## ## fully manual (offline): genes_ldl <- c("54"="SH2B3", "67"="APOE", ...)
##
## ## single locus (SH2B3, locus 54):
 plot_coloc_region(res_ldl, ldl_combine, cad_combine,
                   locus = 54,
                   coloc_row = coloc_ldl[coloc_ldl$locus == 54, ],
                   locus_genes = genes_ldl,
                   exp_name = "LDL-C", out_name = "CAD",
                   out_pdf = "Fig_coloc_LDL_CAD_locus54_SH2B3.pdf")

## ## batch: all decoded loci:
 for (L in coloc_ldl$locus) {
   plot_coloc_region(res_ldl, ldl_combine, cad_combine,
                     locus = L,
                     coloc_row = coloc_ldl[coloc_ldl$locus == L, ],
                     locus_genes = genes_ldl,
                     exp_name = "LDL-C", out_name = "CAD",
                     out_pdf = sprintf("Fig_coloc_LDL_CAD_locus%03d.pdf", L))
 }
##
 ## Example 2: BMI -> T2D (21 decoded loci)
 ## ============================================================
 ##
 ## ## prerequisites: the mechanism script has been run, so the environment has:
 ## ##   res_bmi               (run_eccmr_pair result)
 ## ##   bmi_combine, t2d_combine  (full sumstats)
 ## ##   .region_pair()           (function from this script)
 ## ##   locus_genes_bmi          (manually verified gene names for the 21 loci, from this script)
 ## coloc_bmi <- read.csv("coloc_BMI_T2D_decoded_loci.csv")
 ##
 ## ## use the manually verified table directly; no internet needed:
  genes_bmi <- locus_genes_bmi
 ##
 ## ## single locus (TCF7L2, locus 453):
  plot_coloc_region(res_bmi, bmi_combine, t2d_combine,
                    locus = 453,
                    coloc_row = coloc_bmi[coloc_bmi$locus == 453, ],
                    locus_genes = genes_bmi,
                    exp_name = "BMI", out_name = "T2D",
                    out_pdf = "Fig_coloc_BMI_T2D_locus453_TCF7L2.pdf")
 
 ## ## batch: all 21 decoded loci (gene names included in file names):
  for (L in coloc_bmi$locus) {
    g <- genes_bmi[as.character(L)]
    tag <- if (!is.na(g) && nzchar(g)) paste0(L, "_", g) else L
    plot_coloc_region(res_bmi, bmi_combine, t2d_combine,
                      locus = L,
                      coloc_row = coloc_bmi[coloc_bmi$locus == L, ],
                      locus_genes = genes_bmi,
                      exp_name = "BMI", out_name = "T2D",
                      out_pdf = sprintf("Fig_coloc_BMI_T2D_locus%s.pdf", tag))
  }
 ##
 ## ## only loci judged as genuine pleiotropy by coloc (PP.H4 > 0.8; candidates for main figures):
 for (L in coloc_bmi$locus[coloc_bmi$PP.H4 > 0.8]) {
    g <- genes_bmi[as.character(L)]
    tag <- if (!is.na(g) && nzchar(g)) paste0(L, "_", g) else L
    plot_coloc_region(res_bmi, bmi_combine, t2d_combine,
                      locus = L,
                      coloc_row = coloc_bmi[coloc_bmi$locus == L, ],
                      locus_genes = genes_bmi,
                      exp_name = "BMI", out_name = "T2D",
                      out_pdf = sprintf("Fig_coloc_BMI_T2D_locus%s.pdf", tag))
 }
  
  
