# =====================================================================
# ECC-MR real-data analysis script
#   Application 1: intraocular pressure (IOP, ukb-b-14146) -> glaucoma (finn-b-H7_GLAUCOMA)
#   Application 2: BMI (bmi_combine) -> type 2 diabetes (t2d_combine)
#
# Input: four data frames already loaded in the R session
#   iop_combine, glaucoma_combine, bmi_combine, t2d_combine
#   columns in each: chr, pos, rsid, ref, alt, beta, se, pval, maf, N
#   (convention: beta is the effect of the alt allele)

# LD reference panel: 1000 Genomes EUR (PLINK binary format); see comments in get_ld_matrix()
# =====================================================================

library(eccmr)

## ------------------------------------------------------------------
## 0. Parameters
## ------------------------------------------------------------------
P_GWS      <- 5e-8      # genome-wide significance threshold for instruments
MAF_MIN    <- 0.01      # exclude rare variants
LOCUS_WIN  <- 250000    # locus window: +/-250 kb
MAX_PER_LOCUS <- 20     # max SNPs kept per locus (by P value, to bound the LD matrix size)
PALINDROME <- c(0.42, 0.58)  # MAF range for excluding palindromic SNPs
N_BOOT     <- 500       # number of ECC-MR LD-block bootstrap replicates
R2_THRESH  <- 0.3       # LD block partitioning threshold

## ------------------------------------------------------------------
## 1. data quality control
## ------------------------------------------------------------------
qc_sumstats <- function(df, name) {
  need <- c("chr","pos","rsid","ref","alt","beta","se","pval","maf","N")
  miss <- setdiff(need, colnames(df))
  if (length(miss)) stop(name, " is missing columns: ", paste(miss, collapse = ", "))
  df <- df[!is.na(df$beta) & !is.na(df$se) & !is.na(df$pval) &
             !is.na(df$chr) & !is.na(df$pos) & !is.na(df$maf), ]
  df <- df[df$se > 0 & is.finite(df$se), ]
  df <- df[nchar(df$ref) == 1 & nchar(df$alt) == 1, ]   # keep biallelic SNPs only
  df <- df[!duplicated(df$rsid), ]
  cat(sprintf("[%s] %d variants, chr %s, median N = %g\n",
              name, nrow(df), paste(sort(unique(df$chr)), collapse = ","),
              median(df$N, na.rm = TRUE)))
  df
}

## ------------------------------------------------------------------
## 2. instrument selection: P < 5e-8, MAF filter, locus-window grouping
##    no LD clumping -- ECC-MR works directly with correlated instruments
## ------------------------------------------------------------------
select_instruments <- function(exp, p_thr = P_GWS, maf_min = MAF_MIN,
                               win = LOCUS_WIN, max_per_locus = MAX_PER_LOCUS,
                               f_min = 10) {
  hit <- exp[exp$pval < p_thr & exp$maf >= maf_min & exp$maf <= 1 - maf_min, ]
  ## instrument-strength filter: F = (beta/se)^2 > f_min
  ## (for SNPs with P < 5e-8, z > 5.48 implies F > 30, so this filter usually passes automatically,
  ##  but it is applied explicitly to match the Methods description and keep weak instruments out)
  if (f_min > 0 && all(c("beta", "se") %in% names(hit))) {
    f_stat <- (hit$beta / hit$se)^2
    n_weak <- sum(!is.finite(f_stat) | f_stat <= f_min)
    hit <- hit[is.finite(f_stat) & f_stat > f_min, ]
    hit$f_stat <- (hit$beta / hit$se)^2
    cat(sprintf("Instrument-strength filter: dropped %d weak instruments with F <= %.0f (%d remain, min F = %.1f)\n",
                f_min, n_weak, nrow(hit),
                if (nrow(hit)) min(hit$f_stat) else NA))
  }
  hit <- hit[order(hit$chr, hit$pos), ]
  # greedy locus definition: consecutive same-chromosome segments with gaps < win form one locus
  locus_id <- integer(nrow(hit)); lid <- 0
  for (i in seq_len(nrow(hit))) {
    if (i == 1 || hit$chr[i] != hit$chr[i-1] || hit$pos[i] - hit$pos[i-1] > win) {
      lid <- lid + 1
    }
    locus_id[i] <- lid
  }
  hit$locus <- locus_id
  # keep the top max_per_locus SNPs by P value per locus
  hit <- do.call(rbind, lapply(split(hit, hit$locus), function(x) {
    x <- x[order(x$pval), ]
    if (nrow(x) > max_per_locus) x <- x[seq_len(max_per_locus), ]
    x
  }))
  hit <- hit[order(hit$chr, hit$pos), ]
  cat(sprintf("Instruments: %d GWS SNPs across %d loci\n",
              nrow(hit), length(unique(hit$locus))))
  hit
}

## ------------------------------------------------------------------
## 3. Harmonise: rsID matching + allele alignment + palindromic-SNP removal
##    beta is assumed to correspond to the alt allele in both datasets
## ------------------------------------------------------------------
harmonise_pair <- function(exp_hit, out) {
  m <- merge(exp_hit, out, by = "rsid", suffixes = c("_X", "_Y"))
  cat(sprintf("After rsID matching: %d / %d SNPs\n", nrow(m), nrow(exp_hit)))

  same  <- m$ref_X == m$ref_Y & m$alt_X == m$alt_Y
  flip  <- m$ref_X == m$alt_Y & m$alt_X == m$ref_Y
  drop0 <- !(same | flip)
  cat(sprintf("Dropped (alleles unmatchable): %d\n", sum(drop0)))
  m <- m[same | flip, ]
  flip <- flip[same | flip]

  m$beta_Y[flip] <- -m$beta_Y[flip]
  m$maf_Y[flip]  <- 1 - m$maf_Y[flip]

  # palindromic SNPs (A/T, C/G) with MAF near 0.5: strand unreliable, drop
  pal <- (paste0(m$ref_X, m$alt_X) %in% c("AT","TA","CG","GC"))
  amb <- pal & ((m$maf_X > PALINDROME[1] & m$maf_X < PALINDROME[2]) |
                (m$maf_Y > PALINDROME[1] & m$maf_Y < PALINDROME[2]))
  cat(sprintf("Dropped (ambiguous palindromes): %d\n", sum(amb)))
  m <- m[!amb, ]

  # non-palindromic SNPs with large allele-frequency differences (possible strand mismatch), gently dropped
  badfreq <- abs(m$maf_X - m$maf_Y) > 0.25
  m <- m[!badfreq, ]
  cat(sprintf("Retained after harmonisation: %d SNPs\n", nrow(m)))
  m
}

## ------------------------------------------------------------------
## 4. LD matrix: local 1000G EUR panel (preferred) or the ieugwasr API (fallback)
##
##   local: PANEL_DIR points to the panel folder; the function searches it recursively
##     for .bed/.bim/.fam triplets (both a single panel and per-chromosome splits are supported,
##     e.g. 1000G.EUR.QC or 1000G.EUR.QC.1 ... .22); PLINK 1.9 and 2 both work
##   API: with use = "api", ieugwasr::ld_matrix() is used; a token is required:
##     1. log in with GitHub at https://api.opengwas.io/profile/ and generate a token
##     2. add to ~/.Renviron:  OPENGWAS_JWT=your_token
##     (the API caps the number of variants per query; query per locus and assemble
##      a block-diagonal matrix -- cross-locus LD is approximately zero and does not affect partition_ld_blocks)
## ------------------------------------------------------------------
.find_bfiles <- function(panel_dir) {
  beds <- list.files(panel_dir, pattern = "\\.bed$", recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  pref <- sub("\\.bed$", "", beds, ignore.case = TRUE)
  ok <- file.exists(paste0(pref, ".bim")) & file.exists(paste0(pref, ".fam"))
  pref[ok]
}

.get_ld_one <- function(rsids, bfile, plink) {
  out_prefix <- tempfile("eccmr_ld")
  rs_file <- paste0(out_prefix, "_rsids.txt")
  writeLines(unique(rsids), rs_file)
  v2 <- grepl("plink2|v2", tryCatch(
    system2(plink, "--version", stdout = TRUE, stderr = TRUE)[1],
    error = function(e) ""), ignore.case = TRUE)
  # since PLINK2 v2.0.0-a.7.4 (2026), --r is split into --r-phased / --r-unphased;
  # MR needs the r of genotype dosages, i.e. --r-unphased
  rflag <- if (v2) "--r-unphased" else "--r"
  args <- c("--bfile", bfile, "--extract", rs_file, rflag, "square",
            if (v2) "--write-snplist", "--out", out_prefix)
  if (system2(plink, args) != 0) stop("PLINK run failed: ", bfile)
  ids <- if (v2) {
    read.table(paste0(out_prefix, ".snplist"), stringsAsFactors = FALSE)[[1]]
  } else {
    # PLINK 1.9: --r square outputs the extracted variants in bim order
    bim <- read.table(paste0(bfile, ".bim"), stringsAsFactors = FALSE)
    bim[[2]][bim[[2]] %in% unique(rsids)]
  }
  ld_file <- paste0(out_prefix, ".ld")
  if (!file.exists(ld_file)) {
    # some PLINK2 builds use other extensions; search as a fallback
    cand <- list.files(dirname(out_prefix),
                       pattern = paste0("^", basename(out_prefix)),
                       full.names = TRUE)
    # newer PLINK2 outputs .unphased.vcor1 / .phased.vcor1
    cand <- cand[grepl("\\.(ld|vcor1?|cor|matrix)(\\.gz)?$", cand)]
    if (!length(cand)) stop("PLINK LD output not found: ", out_prefix)
    ld_file <- cand[1]
  }
  ld <- as.matrix(read.table(ld_file, check.names = FALSE))
  if (nrow(ld) != length(ids))
    stop("PLINK output does not match the number of rsIDs; try PLINK2 and check the panel")
  rownames(ld) <- colnames(ld) <- ids
  keep <- rsids[rsids %in% ids]
  ld[keep, keep, drop = FALSE]
}

get_ld_matrix <- function(rsids, locus = NULL, chr = NULL,
                          use = c("panel", "api"),
                          panel_dir = NULL, plink = "plink2",
                          pop = "EUR", batch = 300L) {
  use <- match.arg(use)
  out <- matrix(0, length(rsids), length(rsids),
                dimnames = list(rsids, rsids))

  if (use == "panel") {
    bfiles <- .find_bfiles(panel_dir)
    if (!length(bfiles)) stop("No .bed/.bim/.fam panel found in ", panel_dir)
    message("Panel found: ", paste(basename(bfiles), collapse = ", "))
    if (length(bfiles) > 1 && !is.null(chr)) {
      # per-chromosome split layout: match chromosomes by the bim file-name suffix
      for (cc in unique(chr)) {
        ix <- which(chr == cc)
        # the chromosome number is the last numeric segment of the file name (e.g. 1000G.EUR.QC.4 or xxx.chr4)
        cand <- bfiles[grepl(paste0("[._](chr)?", cc, "$"), bfiles,
                             ignore.case = TRUE)]
        bf <- if (length(cand)) cand[1] else bfiles[1]
        message(sprintf("chr%s uses panel %s (%d SNPs)", cc, basename(bf), length(ix)))
        m <- .get_ld_one(rsids[ix], bf, plink)
        out[rownames(m), colnames(m)] <- m
      }
    } else {
      m <- .get_ld_one(rsids, bfiles[1], plink)
      out[rownames(m), colnames(m)] <- m
    }
  } else {  # api
    if (!requireNamespace("ieugwasr", quietly = TRUE))
      stop("Please install first: install.packages('ieugwasr')")
    if (Sys.getenv("OPENGWAS_JWT") == "")
      stop("OPENGWAS_JWT not found; please configure the token in ~/.Renviron first")
    chunks <- if (is.null(locus)) list(rsids) else split(rsids, locus)
    for (i in seq_along(chunks)) {
      for (b in split(chunks[[i]], ceiling(seq_along(chunks[[i]]) / batch))) {
        message(sprintf("LD query: batch %d/%d (%d SNPs)", i, length(chunks), length(b)))
        m <- tryCatch(ieugwasr::ld_matrix(b, pop = pop, with_alleles = FALSE),
                      error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
        if (is.null(m)) next
        nm <- sub("_[ATCG]+_[ATCG]+$", "", rownames(m))
        ok <- nm %in% rsids
        m <- m[ok, ok, drop = FALSE]; nm <- nm[ok]
        out[nm, nm] <- m
      }
    }
  }
  missing <- rsids[rowSums(out != 0) == 0 & diag(out) == 0]
  diag(out) <- 1
  if (length(missing))
    message("Note: ", length(missing), " SNPs not found in the reference panel (treated as independent)")
  out
}
## ------------------------------------------------------------------
## 5. complete ECC-MR analysis for one exposure-outcome pair
## ------------------------------------------------------------------
run_eccmr_pair <- function(exp_df, out_df, pair_name,
                           ld_mat = NULL, ld_use = c("panel", "api"),
                           panel_dir = NULL, plink = "plink2",
                           pop = "EUR", seed = 1,
                           cores = max(1L, parallel::detectCores() - 1L)) {
  ld_use <- match.arg(ld_use)
  cat("\n================ ", pair_name, " ================\n")
  exp_df <- qc_sumstats(exp_df, paste0(pair_name, " exposure"))
  out_df <- qc_sumstats(out_df, paste0(pair_name, " outcome"))

  hit <- select_instruments(exp_df)
  m   <- harmonise_pair(hit, out_df)

  if (is.null(ld_mat)) {
    ld_mat <- get_ld_matrix(m$rsid, locus = m$locus, chr = m$chr_X,
                            use = ld_use, panel_dir = panel_dir,
                            plink = plink, pop = pop)
  }
  m <- m[m$rsid %in% rownames(ld_mat), ]
  ld_mat <- ld_mat[m$rsid, m$rsid, drop = FALSE]

  blocks <- partition_ld_blocks(ld_mat, r2_thresh = R2_THRESH)
  cat(sprintf("LD blocks: %d (r^2 > %.1f)\n", length(blocks), R2_THRESH))

  ## ECC-MR main analysis (multi-core bootstrap; results independent of cores)
  fit <- eccmr(beta_X = m$beta_X, beta_Y = m$beta_Y, se_Y = m$se_Y,
               blocks = blocks, ld_mat = ld_mat,
               snp_names = m$rsid, n_boot = N_BOOT,
               seed = seed, cores = cores, verbose = TRUE)

  ## the other four comparators use LD-pruned SNPs (one per LD block),
  ## which is their correct usage; ECC-MR always uses all correlated SNPs in the main analysis.
  ## (requires eccmr >= 0.2.2)
  cmp <- compare_mr_methods(beta_X = m$beta_X, beta_Y = m$beta_Y,
                            se_Y = m$se_Y, blocks = blocks, ld_mat = ld_mat,
                            methods = c("IVW", "MR-Egger",
                                        "Weighted median", "MR-Lasso"),
                            prune = TRUE, n_boot = N_BOOT,
                            seed = seed, verbose = TRUE)

  ## pleiotropic locus list
  ## v0.2.1: fit$pleiotropic_index = index vector; fit$pleiotropic = named alpha vector
  pidx <- if (!is.null(fit$pleiotropic_index)) fit$pleiotropic_index
          else if (is.list(fit$pleiotropic)) fit$pleiotropic$index
          else match(names(fit$pleiotropic), m$rsid)
  pidx <- pidx[!is.na(pidx)]
  if (length(pidx)) {
    pleio <- m[pidx, c("rsid","chr_X","pos_X","locus")]
    pleio$alpha_hat <- fit$alpha[pidx]
    pleio_locus <- aggregate(rsid ~ locus, pleio, length)
    colnames(pleio_locus)[2] <- "n_pleiotropic_snps"
  } else {
    pleio <- data.frame(rsid = character(0), chr_X = integer(0),
                        pos_X = integer(0), locus = integer(0),
                        alpha_hat = numeric(0))
    pleio_locus <- data.frame(locus = integer(0), n_pleiotropic_snps = integer(0))
  }

  list(name = pair_name, data = m, ld_mat = ld_mat, blocks = blocks,
       fit = fit, compare = cmp, pleiotropic_snps = pleio,
       pleiotropic_loci = pleio_locus)
}

## ------------------------------------------------------------------
## 6. summary table of results (for the paper table / forest plot)
## ------------------------------------------------------------------
extract_estimates <- function(res) {
  f <- res$fit
  out <- data.frame(
    pair = res$name, method = "ECC-MR",
    theta = f$theta, se = f$se,
    ci_lo = f$ci[1], ci_hi = f$ci[2],
    stringsAsFactors = FALSE)
  ## compare_mr_methods returns a data.frame directly (one row per method);
  ## older versions may return list(results = ...) -- handled for compatibility
  cmp <- if (is.data.frame(res$compare)) res$compare else res$compare$results
  if (!is.null(cmp) && nrow(cmp)) {
    r <- cmp
    keep <- r$method != "ECC-MR"
    ## ci_lo/ci_hi are computed in the package; use them preferentially
    ci_lo <- if ("ci_lo" %in% names(r)) r$ci_lo[keep] else r$theta[keep] - 1.96*r$se[keep]
    ci_hi <- if ("ci_hi" %in% names(r)) r$ci_hi[keep] else r$theta[keep] + 1.96*r$se[keep]
    out <- rbind(out, data.frame(pair = res$name, method = r$method[keep],
                                 theta = r$theta[keep], se = r$se[keep],
                                 ci_lo = ci_lo, ci_hi = ci_hi))
  }
  out
}

## ------------------------------------------------------------------
## 7. Nature Communications style forest plot
##    (Arial/Helvetica, theme_classic, faceted panel labels, restrained two-colour scheme)
## ------------------------------------------------------------------
plot_forest_nc <- function(est_df, xlab, out_pdf, ecc_col = "#B2182B") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible())
  library(ggplot2)
  lv <- c("ECC-MR","IVW","MR-Egger","Weighted median","MR-Lasso",
          "MR-PRESSO","MR-Clust","CAUSE")
  est_df$method <- factor(est_df$method,
                          levels = lv[lv %in% est_df$method])
  p <- ggplot(est_df, aes(y = method, x = theta)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60",
               linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18,
                   linewidth = 0.5, colour = "grey35") +
    geom_point(aes(colour = method == "ECC-MR"), size = 2.4, shape = 15) +
    scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = ecc_col),
                        guide = "none") +
    facet_wrap(~pair, ncol = 1, scales = "free_x") +
    labs(x = xlab, y = NULL) +
    theme_classic(base_size = 9, base_family = "Arial") +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 9, hjust = 0),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          axis.text = element_text(colour = "black"))
  ## fonts: use the generic R sans family (text is silently dropped if Arial is not registered)
  fam <- "sans"
  p <- p + theme_classic(base_size = 9, base_family = fam) %+replace%
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 9, hjust = 0),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          axis.text = element_text(colour = "black"))
  ok <- tryCatch({ ggsave(out_pdf, p, width = 89, height = 60, units = "mm",
                          device = cairo_pdf); TRUE },
                 error = function(e) { message("PDF failed, saving PNG instead: ",
                                               conditionMessage(e)); FALSE })
  if (!ok) ggsave(sub("\\.pdf$", ".png", out_pdf), p,
                  width = 89, height = 60, units = "mm", dpi = 300)
  p
}

## ==================================================================
## Main workflow
## ==================================================================

## local panel: point to your panel folder (searched recursively for .bed/.bim/.fam triplets;
## per-chromosome splits such as 1000G.EUR.QC.1 ... .22 are supported)
PANEL_DIR <- "path/to/1000G_EUR_Phase3_plink"  # set to your local 1000G EUR panel folder
PLINK     <- "plink2"  # or the full path to the plink2 executable
## to use the OpenGWAS API instead, change ld_use below to "api"

res_iop <- run_eccmr_pair(iop_combine, glaucoma_combine,
                          pair_name = "IOP -> Glaucoma",
                          ld_use = "panel", panel_dir = PANEL_DIR,cores=2,
                          plink = PLINK, seed = 101)

res_bmi <- run_eccmr_pair(bmi_combine, t2d_combine,
                          pair_name = "BMI -> T2D",
                          ld_use = "panel", panel_dir = PANEL_DIR,cores=2,
                          plink = PLINK, seed = 202)
res_ldl <- run_eccmr_pair(ldl_combine, cad_combine,
                          pair_name = "LDL -> CAD",
                          ld_use = "panel", panel_dir = PANEL_DIR,cores=2,
                          plink = PLINK, seed = 303)
## summary output
est_all <- do.call(rbind, list(extract_estimates(res_iop),
                               extract_estimates(res_bmi),
                               extract_estimates(res_ldl)))
write.csv(est_all, "eccmr_realdata_estimates.csv", row.names = FALSE)
write.csv(res_iop$pleiotropic_snps, "eccmr_pleiotropic_IOP_Glaucoma.csv",
          row.names = FALSE)
write.csv(res_bmi$pleiotropic_snps, "eccmr_pleiotropic_BMI_T2D.csv",
          row.names = FALSE)
write.csv(res_ldl$pleiotropic_snps, "eccmr_pleiotropic_LDL_CAD.csv",
          row.names = FALSE)
plot_forest_nc(est_all,
               xlab = "Causal estimate (log-OR per unit exposure)",
               out_pdf = "Fig_ECCMR_realdata_forest.pdf")

## print main results
summary(res_iop$fit)
summary(res_bmi$fit)
summary(res_ldl$fit)
# =====================================================================
# Notes
# 1. the LD panel is searched recursively under PANEL_DIR for .bed/.bim/.fam; both a single panel and
#    per-chromosome splits are supported (file names like ...EUR.QC.1 / ...EUR.QC.chr1 match automatically).
#    if panel SNPs are not named by rsID (e.g. chr:pos), rename them with PLINK --update-name first,
#    or use ld_use = "api" with OpenGWAS (requires an OPENGWAS_JWT token).
# 2. ukb-b-14146 and finn-b-H7_GLAUCOMA come from different cohorts (UKB vs FinnGen),
#    so sample overlap is negligible; if both samples of BMI->T2D include UKB, state it in the paper.
# 3. N_BOOT = 500 matches the paper; set 50 for initial debugging.
# =====================================================================

## ------------------------------------------------------------------
## 5c. extended comparison: MR-PRESSO and CAUSE
##
##   one-time installation:
##     install.packages("MRPRESSO")
##     install.packages("devtools")   # or remotes
##     devtools::install_github("jean997/cause")
##     devtools::install_github("MRCIEU/mrclust")
##
##   same rule as compare_mr_methods: these methods also assume independent SNPs,
##   so LD-pruned instruments are used (the SNP with the smallest exposure P per LD block).
## ------------------------------------------------------------------

## keep 1 SNP per LD block (consistent with prune = TRUE in compare)
.prune_one_per_locus <- function(m) {
  m <- m[order(m$pval_X), ]
  m[!duplicated(m$locus), ]
}

## ---- MR-PRESSO: global pleiotropy test + outlier-corrected causal estimate ----
run_mrpresso <- function(res, n_dist = 1000, sig = 0.05) {
  if (!requireNamespace("MRPRESSO", quietly = TRUE))
    stop("Please run install.packages('MRPRESSO') first")
  d <- .prune_one_per_locus(res$data)
  cat(sprintf("MR-PRESSO (%s): %d pruned SNPs\n", res$name, nrow(d)))
  pr <- MRPRESSO::mr_presso(BetaOutcome  = d$beta_Y, BetaExposure = d$beta_X,
                            SdOutcome    = d$se_Y,   SdExposure   = d$se_X,
                            OUTLIERtest  = TRUE, DISTORTIONtest   = TRUE,
                            data = d, NbDistribution = n_dist,
                            SignifThreshold = sig)
  main <- pr$`Main MR results`
  ## when outliers exist the second row is the corrected estimate; report it preferentially
  use <- if (nrow(main) >= 2 && !is.na(main$`Causal Estimate`[2])) 2L else 1L
  n_out <- tryCatch(
    sum(pr$`MR-PRESSO results`$`Outlier Test`$Pvalue < sig / nrow(d), na.rm = TRUE),
    error = function(e) NA_integer_)
  glob_p <- tryCatch(pr$`MR-PRESSO results`$`Global Test`$Pvalue,
                     error = function(e) NA_real_)
  out <- data.frame(pair = res$name, method = "MR-PRESSO",
                    theta = main$`Causal Estimate`[use],
                    se    = main$Sd[use],
                    ci_lo = main$`Causal Estimate`[use] - 1.96 * main$Sd[use],
                    ci_hi = main$`Causal Estimate`[use] + 1.96 * main$Sd[use],
                    n_outliers = n_out, global_p = glob_p,
                    stringsAsFactors = FALSE)
  cat(sprintf("  theta = %.4f (SE %.4f), outliers = %s, global pleiotropy P = %s\n",
              out$theta, out$se,
              ifelse(is.na(n_out), "NA", n_out), signif(glob_p, 3)))
  out
}

## ---- CAUSE: shared-factor vs causal model (needs genome-wide sumstats for nuisance parameters) ----
##   res    : object returned by run_eccmr_pair() (provides the pruned instruments)
##   exp_df / out_df : full data frames (bmi_combine, t2d_combine etc.),
##                     used to estimate CAUSE genome-wide nuisance parameters (~1e5 SNPs needed)
run_cause <- function(res, exp_df, out_df, n_nuisance = 100000, seed = 1) {
  if (!requireNamespace("cause", quietly = TRUE))
    stop("Please run devtools::install_github('jean997/cause') first")
  fmt <- function(df) data.frame(snp = df$rsid, beta = df$beta, se = df$se,
                                 A1 = df$alt, A2 = df$ref)
  cd <- cause::gwas_merge(fmt(exp_df), fmt(out_df),
                          snp_name_cols = "snp",
                          beta_hat_cols = c("beta", "beta"),
                          se_cols       = c("se", "se"),
                          A1_cols       = c("A1", "A1"),
                          A2_cols       = c("A2", "A2"))
  inst <- .prune_one_per_locus(res$data)$rsid
  set.seed(seed)
  nui <- sample(cd$snp, min(n_nuisance, length(cd$snp)))
  cat(sprintf("CAUSE (%s): nuisance parameters on %d SNPs, fit on %d instruments\n",
              res$name, length(nui), length(inst)))
  params <- cause::est_cause_params(cd, variants = nui)
  cfit   <- cause::cause(cd, param_ests = params, variants = inst)
  print(summary(cfit)$tab)   # full sharing/causal model results, for cross-checking

  tab <- summary(cfit)$tab
  cr  <- tab[grepl("causal", rownames(tab), ignore.case = TRUE), ]
  bcol <- grep("beta", colnames(tab), ignore.case = TRUE)
  ## CAUSE reports posterior quantiles; the median is the point estimate and IQR/1.35 an approximate SE
  med <- as.numeric(cr[bcol[grepl("0.5|med", colnames(tab)[bcol])[1]]])
  q25 <- suppressWarnings(as.numeric(cr[bcol[grepl("0.25|0.025", colnames(tab)[bcol])[1]]]))
  q75 <- suppressWarnings(as.numeric(cr[bcol[grepl("0.75|0.975", colnames(tab)[bcol])[1]]]))
  se  <- if (!is.na(q25) && !is.na(q75)) (q75 - q25) / 1.349 else NA_real_
  out <- data.frame(pair = res$name, method = "CAUSE",
                    theta = med, se = se,
                    ci_lo = ifelse(is.na(q25), med - 1.96 * se, q25),
                    ci_hi = ifelse(is.na(q75), med + 1.96 * se, q75),
                    n_outliers = NA_integer_, global_p = NA_real_,
                    stringsAsFactors = FALSE)
  cat(sprintf("  causal model beta = %.4f (approximate SE %.4f)\n", out$theta, out$se))
  attr(out, "cause_fit") <- cfit
  out
}

## ---- MR-Clust: pathway clustering (no single point estimate; the largest non-junk cluster is reported) ----
##   MR-Clust answers a different question from point estimators: it clusters instruments by
##   their outcome effects; each cluster is a candidate causal pathway and "junk/null" clusters are pleiotropic noise.
##   to enter the comparison table we report the centre of the largest non-junk cluster (as in the literature),
##   and also output the number of clusters for mechanistic interpretation (cross-check with the ECC-MR locus list).
run_mrclust <- function(res) {
  if (!requireNamespace("mrclust", quietly = TRUE))
    stop("Please run devtools::install_github('MRCIEU/mrclust') first")
  d <- .prune_one_per_locus(res$data)
  cat(sprintf("MR-Clust (%s): %d pruned SNPs\n", res$name, nrow(d)))
  fit <- mrclust::mr_clust_em(theta    = d$beta_Y / d$beta_X,
                              theta_se = d$se_Y / abs(d$beta_X),
                              bx = d$beta_X, by = d$beta_Y,
                              bxse = d$se_X, byse = d$se_Y)
  best <- fit$results$best            # optimal cluster assignment per SNP
  ## cluster centres: exclude junk/null clusters, take the largest weight
  cs <- fit$results$clusters
  if (is.null(cs)) cs <- best
  ## accommodate different output structures across versions: find the table with cluster mean/weight
  tab <- if (all(c("cluster_mean", "cluster_weight") %in% colnames(cs))) cs else NULL
  if (is.null(tab)) {
    message("MR-Clust output structure not recognised; inspect fit$results manually; skipping comparison-table merge")
    print(names(fit$results)); print(utils::head(cs))
    return(invisible(fit))
  }
  ok <- !grepl("junk|null", tab$cluster, ignore.case = TRUE) &
        is.finite(tab$cluster_mean)
  main <- tab[ok, ][which.max(tab$cluster_weight[ok]), ]
  n_cl <- sum(ok)
  out <- data.frame(pair = res$name, method = "MR-Clust",
                    theta = main$cluster_mean, se = NA_real_,
                    ci_lo = NA_real_, ci_hi = NA_real_,
                    n_outliers = NA_integer_, global_p = NA_real_,
                    stringsAsFactors = FALSE)
  cat(sprintf("  non-junk clusters: %d, main cluster theta = %.4f (weight %.2f)\n",
              n_cl, out$theta, main$cluster_weight))
  attr(out, "mrclust_fit") <- fit
  out
}

## ==================================================================
## optional: run the extended comparisons (uncomment after installing the three packages)
## ==================================================================
## press_iop <- run_mrpresso(res_iop)
## press_bmi <- run_mrpresso(res_bmi)
## press_ldl <- run_mrpresso(res_ldl)
## clust_iop <- run_mrclust(res_iop)
## clust_bmi <- run_mrclust(res_bmi)
## clust_ldl <- run_mrclust(res_ldl)
## cause_iop <- run_cause(res_iop, iop_combine, glaucoma_combine)
## cause_bmi <- run_cause(res_bmi, bmi_combine, t2d_combine)
## cause_ldl <- run_cause(res_ldl, ldl_combine, cad_combine)
## est_all <- rbind(est_all,
##                  press_iop[, 1:6], press_bmi[, 1:6], press_ldl[, 1:6],
##                  cause_iop[, 1:6], cause_bmi[, 1:6], cause_ldl[, 1:6])
## plot_forest_nc(est_all,
##                xlab = "Causal estimate (log-OR per unit exposure)",
##                out_pdf = "Fig_forest_ext.pdf")
