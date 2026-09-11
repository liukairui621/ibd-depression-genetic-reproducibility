#!/usr/bin/env Rscript

set.seed(20260817)
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(susieR)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: 06_run_qtl_coloc.R <project_root>")
root <- args[[1]]
out <- file.path(root,"results","qtl")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
long <- fread(file.path(root,"data","qtl_prepared","qtl_prepared_long.tsv.gz"))
manifest <- fread(file.path(out,"eligible_qtl_traits.tsv"))
samples <- fread(file.path(root,"data","prepared","trait_sample_sizes.tsv"))
gwas_traits <- c("Howard_DEP","FinnGen_DEP","deLange_IBD","FinnGen_IBD")
P12 <- c(1e-5,1e-6)

read_ld <- function(path,n) {
  con <- file(path,"rb"); on.exit(close(con))
  x <- readBin(con,what="numeric",n=n*n,size=8L,endian="little")
  matrix(x,nrow=n,ncol=n,byrow=TRUE)
}

abf_rows <- list(); susie_rows <- list(); fit_rows <- list()
for (locus in unique(manifest$locus)) {
  bim <- fread(file.path(root,"data","ld",paste0(locus,"_1000G_EUR.bim")),header=FALSE,
               col.names=c("chr","rsid","cm","pos","a1","a2"))
  gwas <- fread(file.path(root,"data","prepared",paste0(locus,"_four_gwas_common.tsv.gz")))
  gwas <- gwas[match(bim$rsid,gwas$rsid)]
  freq <- fread(file.path(root,"data","ld",paste0(locus,"_1000G_EUR.frq")))
  maf_ref <- freq$MAF[match(bim$rsid,freq$SNP)]
  R0 <- read_ld(file.path(root,"data","ld",paste0(locus,"_LD.ld.bin")),nrow(bim))
  R0 <- (R0+t(R0))/2; diag(R0)<-1; R0[!is.finite(R0)]<-0
  dimnames(R0)<-list(bim$rsid,bim$rsid)

  for (key in manifest[manifest[["locus"]] == locus, qtl_key]) {
    meta <- manifest[qtl_key==key][1]
    q <- long[qtl_key==key]
    idx <- match(q$rsid,bim$rsid)
    keep <- !is.na(idx)
    q <- q[keep]; idx <- idx[keep]
    ord <- order(idx); q <- q[ord]; idx <- idx[ord]
    if (nrow(q)<200) next
    R <- R0[idx,idx,drop=FALSE]
    maf <- maf_ref[idx]
    R_ridge <- 0.9999*R + 0.0001*diag(nrow(R))
    dimnames(R_ridge)<-dimnames(R)
    qds <- list(beta=q$aligned_beta,varbeta=q$aligned_se^2,snp=q$rsid,
                position=q$pos37,N=meta$sample_size,type="quant",MAF=maf,LD=R_ridge)

    for (trait_name in gwas_traits) {
      srow <- samples[trait==trait_name][1]
      gd <- gwas[idx]
      gds <- list(beta=gd[[paste0(trait_name,"_beta")]],varbeta=gd[[paste0(trait_name,"_se")]]^2,
                  snp=gd$rsid,position=gd$pos,N=srow$N,type="cc",s=srow$cases/srow$N,
                  MAF=maf,LD=R_ridge)
      primary_abf <- NULL
      for (prior in P12) {
        abf <- tryCatch(coloc.abf(gds,qds,p1=1e-4,p2=1e-4,p12=prior),error=function(e)e)
        if (inherits(abf,"error")) {
          abf_rows[[length(abf_rows)+1L]] <- data.table(qtl_key=key,locus=locus,gene_symbol=meta$gene_symbol,
            gene_id=meta$gene_id,dataset_label=meta$dataset_label,qtl_type=meta$qtl_type,gwas_trait=trait_name,
            p12=prior,status="error",message=conditionMessage(abf))
        } else {
          s <- abf$summary
          row <- data.table(qtl_key=key,locus=locus,gene_symbol=meta$gene_symbol,gene_id=meta$gene_id,
            dataset_label=meta$dataset_label,qtl_type=meta$qtl_type,gwas_trait=trait_name,p12=prior,
            nsnps=unname(s[["nsnps"]]),PP.H0=unname(s[["PP.H0.abf"]]),PP.H1=unname(s[["PP.H1.abf"]]),
            PP.H2=unname(s[["PP.H2.abf"]]),PP.H3=unname(s[["PP.H3.abf"]]),PP.H4=unname(s[["PP.H4.abf"]]),
            status="ok",message=NA_character_)
          abf_rows[[length(abf_rows)+1L]] <- row
          if (prior==1e-5) primary_abf <- row
        }
      }
      gate <- !is.null(primary_abf) && (primary_abf$PP.H4>=0.5 || primary_abf$PP.H3>=0.5)
      if (!gate) {
        susie_rows[[length(susie_rows)+1L]] <- data.table(qtl_key=key,locus=locus,gene_symbol=meta$gene_symbol,
          gene_id=meta$gene_id,dataset_label=meta$dataset_label,qtl_type=meta$qtl_type,gwas_trait=trait_name,
          p12=NA_real_,PP.H4=NA_real_,status="not_screened_by_abf_gate",message=NA_character_)
        next
      }
      qfit <- tryCatch(runsusie(qds,L=10,coverage=0.95,maxit=2000,repeat_until_convergence=FALSE,
                                estimate_residual_variance=FALSE),error=function(e)e)
      gfit <- tryCatch(runsusie(gds,L=10,coverage=0.95,maxit=2000,repeat_until_convergence=FALSE,
                                estimate_residual_variance=FALSE),error=function(e)e)
      fit_rows[[length(fit_rows)+1L]] <- data.table(qtl_key=key,locus=locus,gene_symbol=meta$gene_symbol,
        dataset_label=meta$dataset_label,gwas_trait=trait_name,n_snps=nrow(q),
        qtl_converged=!inherits(qfit,"error")&&isTRUE(qfit$converged),
        gwas_converged=!inherits(gfit,"error")&&isTRUE(gfit$converged),
        qtl_n_cs=if(inherits(qfit,"error")||is.null(qfit$sets$cs)) NA_integer_ else length(qfit$sets$cs),
        gwas_n_cs=if(inherits(gfit,"error")||is.null(gfit$sets$cs)) NA_integer_ else length(gfit$sets$cs))
      if (inherits(qfit,"error")||inherits(gfit,"error")) {
        susie_rows[[length(susie_rows)+1L]] <- data.table(qtl_key=key,locus=locus,gene_symbol=meta$gene_symbol,
          gene_id=meta$gene_id,dataset_label=meta$dataset_label,qtl_type=meta$qtl_type,gwas_trait=trait_name,
          p12=NA_real_,PP.H4=NA_real_,status="fit_error",message=paste(if(inherits(qfit,"error"))conditionMessage(qfit) else "",if(inherits(gfit,"error"))conditionMessage(gfit) else ""))
        next
      }
      for (prior in P12) {
        res <- tryCatch(coloc.susie(gfit,qfit,p1=1e-4,p2=1e-4,p12=prior),error=function(e)e)
        if (inherits(res,"error")||is.null(res$summary)||nrow(res$summary)==0) {
          susie_rows[[length(susie_rows)+1L]] <- data.table(qtl_key=key,locus=locus,gene_symbol=meta$gene_symbol,
            gene_id=meta$gene_id,dataset_label=meta$dataset_label,qtl_type=meta$qtl_type,gwas_trait=trait_name,
            p12=prior,PP.H4=NA_real_,status="not_estimable",message=if(inherits(res,"error"))conditionMessage(res) else "no_signal_pair")
        } else {
          sm <- as.data.table(res$summary)
          setnames(sm,c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf"),
                   c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"),skip_absent=TRUE)
          sm[,`:=`(qtl_key=key,locus=locus,gene_symbol=meta$gene_symbol,gene_id=meta$gene_id,
                   dataset_label=meta$dataset_label,qtl_type=meta$qtl_type,gwas_trait=trait_name,p12=prior,
                   status="ok",message=NA_character_)]
          susie_rows[[length(susie_rows)+1L]] <- sm[,.(qtl_key,locus,gene_symbol,gene_id,dataset_label,qtl_type,
            gwas_trait,p12,nsnps,PP.H0,PP.H1,PP.H2,PP.H3,PP.H4,idx1,idx2,status,message)]
        }
      }
    }
  }
}

fwrite(rbindlist(abf_rows,fill=TRUE),file.path(out,"qtl_coloc_abf.tsv"),sep="\t")
fwrite(rbindlist(susie_rows,fill=TRUE),file.path(out,"qtl_coloc_susie.tsv"),sep="\t")
fwrite(rbindlist(fit_rows,fill=TRUE),file.path(out,"qtl_susie_fit_qc.tsv"),sep="\t")
