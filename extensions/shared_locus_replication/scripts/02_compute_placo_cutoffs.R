#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

root <- "/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication"
source(file.path(root, "software/PLACO/PLACO_v0.2.0.R"))

stats <- read.delim(file.path(root, "results/raw/pair_harmonization_stats.tsv"), check.names = FALSE)
target <- 1e-4

p_for_product <- function(x, varz, method, product_sign) {
  z1 <- sqrt(x)
  z2 <- if (product_sign == "positive") sqrt(x) else -sqrt(x)
  if (method == "placo_plus") {
    stop("CorZ must be passed for PLACO+")
  }
  placo(c(z1, z2), VarZ = varz)$p.placo
}

p_for_product_plus <- function(x, varz, corz, product_sign) {
  z1 <- sqrt(x)
  z2 <- if (product_sign == "positive") sqrt(x) else -sqrt(x)
  placo.plus(c(z1, z2), VarZ = varz, CorZ = corz)$p.placo.plus
}

find_cutoff <- function(varz, corz, method, product_sign) {
  fn <- if (method == "placo_plus") {
    function(x) p_for_product_plus(x, varz, corz, product_sign) - target
  } else {
    function(x) p_for_product(x, varz, method, product_sign) - target
  }
  upper <- 20
  while (fn(upper) > 0 && upper < 200) upper <- upper * 1.5
  if (fn(upper) > 0) stop("Could not bracket product cutoff")
  uniroot(fn, interval = c(0.05, upper), tol = 1e-10)$root
}

out <- vector("list", nrow(stats))
for (i in seq_len(nrow(stats))) {
  varz <- c(stats$var_z_dep[i], stats$var_z_ibd[i])
  method <- stats$method[i]
  corz <- stats$cor_z[i]
  pos <- find_cutoff(varz, corz, method, "positive")
  neg <- find_cutoff(varz, corz, method, "negative")

  grid <- seq(0.5, max(pos, neg) + 10, length.out = 80)
  eval_grid <- function(product_sign) {
    if (method == "placo_plus") {
      vapply(grid, p_for_product_plus, numeric(1), varz = varz, corz = corz, product_sign = product_sign)
    } else {
      vapply(grid, p_for_product, numeric(1), varz = varz, method = method, product_sign = product_sign)
    }
  }
  ppos <- eval_grid("positive")
  pneg <- eval_grid("negative")
  if (any(diff(ppos) > 1e-10) || any(diff(pneg) > 1e-10)) {
    stop(sprintf("Non-monotone PLACO tail detected for %s", stats$pair_id[i]))
  }

  out[[i]] <- data.frame(
    pair_id = stats$pair_id[i], method = method,
    var_z_dep = varz[1], var_z_ibd = varz[2], cor_z = corz,
    target_p = target, positive_product_cutoff = pos,
    negative_product_abs_cutoff = neg, monotonicity_check = "PASS"
  )
}

out <- do.call(rbind, out)
write.table(out, file.path(root, "results/raw/placo_product_cutoffs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(root, "provenance/R_session_cutoffs.txt"))
