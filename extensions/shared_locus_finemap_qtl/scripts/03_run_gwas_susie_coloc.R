#!/usr/bin/env Rscript

set.seed(20260817)
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 03_run_gwas_susie_coloc.R <project_root>")
root <- args[[1]]
prepared <- file.path(root, "data", "prepared")
ld_dir <- file.path(root, "data", "ld")
out <- file.path(root, "results", "gwas")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

loci <- c("block154", "block464", "block1671")
traits <- c("Howard_DEP", "FinnGen_DEP", "deLange_IBD", "FinnGen_IBD")
pairs <- list(
  c("external_disease_pair", "Howard_DEP", "deLange_IBD"),
  c("finngen_disease_pair", "FinnGen_DEP", "FinnGen_IBD"),
  c("cross_howard_finngen_ibd", "Howard_DEP", "FinnGen_IBD"),
  c("cross_finngen_dep_delange", "FinnGen_DEP", "deLange_IBD"),
  c("depression_cross_cohort", "Howard_DEP", "FinnGen_DEP"),
  c("ibd_cross_cohort", "deLange_IBD", "FinnGen_IBD")
)
sample_sizes <- fread(file.path(prepared, "trait_sample_sizes.tsv"))
L <- 10L
MAXIT <- 2000L
RIDGES <- c(0, 1e-4)
P12 <- c(1e-5, 1e-6)

read_ld <- function(path, n) {
  con <- file(path, "rb")
  on.exit(close(con))
  x <- readBin(con, what = "numeric", n = n * n, size = 8L, endian = "little")
  if (length(x) != n * n) stop("LD binary size mismatch: ", path)
  matrix(x, nrow = n, ncol = n, byrow = TRUE)
}

top_snp <- function(fit, idx) {
  if (is.na(idx) || idx < 1 || idx > nrow(fit$alpha)) return(NA_character_)
  nms <- fit$variable_name
  if (is.null(nms)) nms <- colnames(fit$alpha)
  nms[[which.max(fit$alpha[idx, ])]]
}

make_dataset <- function(dat, trait_name, maf, R) {
  row <- sample_sizes[sample_sizes[["trait"]] == trait_name][1]
  list(
    beta = dat[[paste0(trait_name, "_beta")]],
    varbeta = dat[[paste0(trait_name, "_se")]]^2,
    snp = dat$rsid,
    position = dat$pos,
    N = row$N,
    type = "cc",
    s = row$cases / row$N,
    MAF = maf,
    LD = R
  )
}

empty_susie <- function(locus, pair_name, t1, t2, ridge, p12, status, msg = NA_character_) {
  data.table(locus=locus, pair_type=pair_name, trait1=t1, trait2=t2, ridge=ridge,
             p12=p12, nsnps=NA_integer_, PP.H0=NA_real_, PP.H1=NA_real_,
             PP.H2=NA_real_, PP.H3=NA_real_, PP.H4=NA_real_, idx1=NA_integer_,
             idx2=NA_integer_, trait1_top_snp=NA_character_, trait2_top_snp=NA_character_,
             status=status, message=msg)
}

all_fit_qc <- list()
all_pip <- list()
all_cs <- list()
all_susie_coloc <- list()
all_abf <- list()
all_overlap <- list()
all_ld_qc <- list()

for (locus in loci) {
  dat <- fread(file.path(prepared, paste0(locus, "_four_gwas_common.tsv.gz")))
  bim <- fread(file.path(ld_dir, paste0(locus, "_1000G_EUR.bim")), header=FALSE,
               col.names=c("chr", "rsid", "cm", "pos", "a1", "a2"))
  dat <- dat[match(bim$rsid, dat$rsid)]
  if (anyNA(dat$rsid) || !identical(dat$rsid, bim$rsid)) stop(locus, ": BIM/data mismatch")
  freq <- fread(file.path(ld_dir, paste0(locus, "_1000G_EUR.frq")))
  maf <- freq$MAF[match(bim$rsid, freq$SNP)]
  if (anyNA(maf) || any(maf <= 0 | maf > 0.5)) stop(locus, ": invalid MAF")
  R0 <- read_ld(file.path(ld_dir, paste0(locus, "_LD.ld.bin")), nrow(bim))
  R0 <- (R0 + t(R0)) / 2
  diag(R0) <- 1
  R0[!is.finite(R0)] <- 0
  dimnames(R0) <- list(bim$rsid, bim$rsid)
  eig_min <- min(eigen(R0, symmetric=TRUE, only.values=TRUE)$values)
  all_ld_qc[[length(all_ld_qc)+1L]] <- data.table(
    locus=locus, n_snps=nrow(bim), n_reference=489L,
    max_asymmetry=max(abs(R0-t(R0))), min_diagonal=min(diag(R0)),
    max_diagonal=max(diag(R0)), min_eigenvalue=eig_min
  )

  locus_fits <- list()
  for (ridge in RIDGES) {
    R <- (1-ridge)*R0 + ridge*diag(nrow(R0))
    dimnames(R) <- dimnames(R0)
    fits <- list()
    datasets <- list()
    for (trait_name in traits) {
      ds <- make_dataset(dat, trait_name, maf, R)
      datasets[[trait_name]] <- ds
      mismatch <- tryCatch(
        estimate_s_rss(z=ds$beta/sqrt(ds$varbeta), R=R, n=ds$N),
        error=function(e) NA_real_
      )
      fit <- tryCatch(
        runsusie(ds, L=L, coverage=0.95, maxit=MAXIT,
                 repeat_until_convergence=FALSE, estimate_residual_variance=FALSE),
        error=function(e) structure(list(error=conditionMessage(e)), class="failed_susie")
      )
      fits[[trait_name]] <- fit
      if (inherits(fit, "failed_susie")) {
        all_fit_qc[[length(all_fit_qc)+1L]] <- data.table(
          locus=locus, trait=trait_name, ridge=ridge, converged=FALSE,
          n_credible_sets=NA_integer_, max_pip=NA_real_, mismatch_s=mismatch,
          min_cs_purity=NA_real_, median_cs_size=NA_real_, status="fit_error",
          message=fit$error
        )
      } else {
        cs_count <- if (is.null(fit$sets$cs)) 0L else length(fit$sets$cs)
        purity <- fit$sets$purity
        all_fit_qc[[length(all_fit_qc)+1L]] <- data.table(
          locus=locus, trait=trait_name, ridge=ridge, converged=isTRUE(fit$converged),
          n_credible_sets=cs_count, max_pip=max(fit$pip), mismatch_s=mismatch,
          min_cs_purity=if (is.null(purity)||nrow(purity)==0) NA_real_ else min(purity$min.abs.corr),
          median_cs_size=if (cs_count==0) NA_real_ else median(lengths(fit$sets$cs)),
          status="ok", message=NA_character_
        )
        all_pip[[length(all_pip)+1L]] <- data.table(
          locus=locus, trait=trait_name, ridge=ridge, rsid=dat$rsid,
          chr=dat$chr, pos=dat$pos, pip=fit$pip,
          beta=dat[[paste0(trait_name,"_beta")]], se=dat[[paste0(trait_name,"_se")]]
        )
        if (cs_count > 0) {
          for (cs_idx in seq_along(fit$sets$cs)) {
            members <- fit$sets$cs[[cs_idx]]
            all_cs[[length(all_cs)+1L]] <- data.table(
              locus=locus, trait=trait_name, ridge=ridge, cs_id=cs_idx,
              rsid=dat$rsid[members], chr=dat$chr[members], pos=dat$pos[members],
              pip=fit$pip[members], alpha=fit$alpha[cs_idx,members],
              cs_size=length(members), cs_top_snp=top_snp(fit,cs_idx),
              cs_min_abs_corr=purity[cs_idx,"min.abs.corr"][[1]],
              cs_mean_abs_corr=purity[cs_idx,"mean.abs.corr"][[1]],
              cs_median_abs_corr=purity[cs_idx,"median.abs.corr"][[1]]
            )
          }
        }
      }
    }
    locus_fits[[as.character(ridge)]] <- fits
    saveRDS(fits, file.path(out, paste0(locus,"_susie_fits_ridge_",format(ridge,scientific=TRUE),".rds")))

    for (pair in pairs) {
      pair_name <- pair[[1]]; t1 <- pair[[2]]; t2 <- pair[[3]]
      f1 <- fits[[t1]]; f2 <- fits[[t2]]
      for (prior in P12) {
        if (inherits(f1,"failed_susie") || inherits(f2,"failed_susie")) {
          all_susie_coloc[[length(all_susie_coloc)+1L]] <- empty_susie(locus,pair_name,t1,t2,ridge,prior,"fit_error")
        } else {
          res <- tryCatch(coloc.susie(f1,f2,p1=1e-4,p2=1e-4,p12=prior), error=function(e)e)
          if (inherits(res,"error") || is.null(res$summary) || nrow(res$summary)==0) {
            msg <- if (inherits(res,"error")) conditionMessage(res) else "no_signal_pair"
            all_susie_coloc[[length(all_susie_coloc)+1L]] <- empty_susie(locus,pair_name,t1,t2,ridge,prior,"not_estimable",msg)
          } else {
            sm <- as.data.table(res$summary)
            setnames(sm,c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf"),
                     c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"),skip_absent=TRUE)
            sm[,`:=`(locus=locus,pair_type=pair_name,trait1=t1,trait2=t2,ridge=ridge,p12=prior,
                     trait1_top_snp=vapply(idx1,function(i)top_snp(f1,i),character(1)),
                     trait2_top_snp=vapply(idx2,function(i)top_snp(f2,i),character(1)),
                     status="ok",message=NA_character_)]
            all_susie_coloc[[length(all_susie_coloc)+1L]] <- sm[,.(locus,pair_type,trait1,trait2,ridge,p12,nsnps,
              PP.H0,PP.H1,PP.H2,PP.H3,PP.H4,idx1,idx2,trait1_top_snp,trait2_top_snp,status,message)]
          }
        }
        abf <- tryCatch(coloc.abf(datasets[[t1]],datasets[[t2]],p1=1e-4,p2=1e-4,p12=prior),error=function(e)e)
        if (inherits(abf,"error")) {
          all_abf[[length(all_abf)+1L]] <- data.table(locus=locus,pair_type=pair_name,trait1=t1,trait2=t2,
            ridge=ridge,p12=prior,status="error",message=conditionMessage(abf))
        } else {
          s <- abf$summary
          all_abf[[length(all_abf)+1L]] <- data.table(locus=locus,pair_type=pair_name,trait1=t1,trait2=t2,
            ridge=ridge,p12=prior,nsnps=unname(s[["nsnps"]]),PP.H0=unname(s[["PP.H0.abf"]]),
            PP.H1=unname(s[["PP.H1.abf"]]),PP.H2=unname(s[["PP.H2.abf"]]),
            PP.H3=unname(s[["PP.H3.abf"]]),PP.H4=unname(s[["PP.H4.abf"]]),status="ok",message=NA_character_)
        }
      }
    }
  }

  primary_fits <- locus_fits[["0"]]
  for (pair in list(c("depression_cross_cohort","Howard_DEP","FinnGen_DEP"),c("ibd_cross_cohort","deLange_IBD","FinnGen_IBD"))) {
    pair_name <- pair[[1]]; t1 <- pair[[2]]; t2 <- pair[[3]]
    f1 <- primary_fits[[t1]]; f2 <- primary_fits[[t2]]
    if (!inherits(f1,"failed_susie") && !inherits(f2,"failed_susie") && length(f1$sets$cs)>0 && length(f2$sets$cs)>0) {
      for (i in seq_along(f1$sets$cs)) for (j in seq_along(f2$sets$cs)) {
        m1 <- f1$sets$cs[[i]]; m2 <- f2$sets$cs[[j]]
        exact <- intersect(dat$rsid[m1],dat$rsid[m2])
        top1 <- which.max(f1$alpha[i,]); top2 <- which.max(f2$alpha[j,])
        r2 <- R0[top1,top2]^2
        all_overlap[[length(all_overlap)+1L]] <- data.table(
          locus=locus,pair_type=pair_name,trait1=t1,trait2=t2,cs1=i,cs2=j,
          top1=dat$rsid[top1],top2=dat$rsid[top2],top_variant_r2=r2,
          n_exact_overlap=length(exact),exact_overlap=paste(exact,collapse=";"),
          replicated=length(exact)>0 || r2>=0.8
        )
      }
    } else {
      all_overlap[[length(all_overlap)+1L]] <- data.table(locus=locus,pair_type=pair_name,trait1=t1,trait2=t2,
        cs1=NA_integer_,cs2=NA_integer_,top1=NA_character_,top2=NA_character_,top_variant_r2=NA_real_,
        n_exact_overlap=0L,exact_overlap="",replicated=FALSE)
    }
  }
}

fwrite(rbindlist(all_ld_qc,fill=TRUE),file.path(out,"ld_matrix_qc.tsv"),sep="\t")
fwrite(rbindlist(all_fit_qc,fill=TRUE),file.path(out,"susie_fit_qc.tsv"),sep="\t")
fwrite(rbindlist(all_pip,fill=TRUE),file.path(out,"variant_pip.tsv.gz"),sep="\t")
fwrite(rbindlist(all_cs,fill=TRUE),file.path(out,"credible_set_members.tsv.gz"),sep="\t")
fwrite(rbindlist(all_susie_coloc,fill=TRUE),file.path(out,"gwas_coloc_susie.tsv"),sep="\t")
fwrite(rbindlist(all_abf,fill=TRUE),file.path(out,"gwas_coloc_abf.tsv"),sep="\t")
fwrite(rbindlist(all_overlap,fill=TRUE),file.path(out,"cross_cohort_credible_set_overlap.tsv"),sep="\t")
