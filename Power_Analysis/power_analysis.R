## Power calculations for continuous trait GWAS (additive 1-df)
## Uses noncentral chi-square distribution (matches QUANTO's analytic approach)

power_from_r2 <- function(N, alpha, r2) {
  if (r2 <= 0 || r2 >= 1) stop("r2 must be in (0, 1)")
  lambda <- N * r2 / (1 - r2)         # noncentrality parameter
  crit <- qchisq(1 - alpha, df = 1)   # critical value
  1 - pchisq(crit, df = 1, ncp = lambda)
}

r2_required_for_power <- function(N, alpha, target_power = 0.80,
                                  lower = 1e-10, upper = 0.999999, tol = 1e-12) {
  f <- function(r2) power_from_r2(N, alpha, r2) - target_power

  # bracket check
  if (f(lower) > 0) return(lower)
  if (f(upper) < 0) stop("Target power not achievable with r2 < 1 under these settings.")

  uniroot(f, lower = lower, upper = upper, tol = tol)$root
}

beta_from_r2 <- function(r2, maf) {
  if (maf <= 0 || maf >= 1) stop("maf must be in (0, 1)")
  # For standardized phenotype: r2 ≈ 2*maf*(1-maf)*beta^2
  sqrt(r2 / (2 * maf * (1 - maf)))
}

## ---- Inputs ----
strata <- data.frame(
  stratum = c("Meta", "AFR", "EUR", "AFR_female", "AFR_male", "EUR_female", "EUR_male"),
  N       = c(1002,   573,   429,   335,         238,       152,         277)
)

alphas <- data.frame(
  threshold = c("5e-8", "1e-5"),
  alpha     = c(5e-8,   1e-5)
)

mafs <- c(0.05, 0.10, 0.20, 0.30, 0.50)
target_power <- 0.80

## ---- Compute table ----
make_rows <- function(stratum, N) {
  do.call(rbind, lapply(1:nrow(alphas), function(j) {
    alpha <- alphas$alpha[j]
    thr <- alphas$threshold[j]

    r2_req <- r2_required_for_power(N, alpha, target_power)

    row <- data.frame(
      stratum = stratum,
      N = N,
      threshold = thr,
      alpha = alpha,
      power_target = target_power,
      r2_required = r2_req
    )

    # Add betas at selected MAFs
    for (m in mafs) {
      nm <- paste0("beta_MAF_", format(m, nsmall = 2))
      row[[nm]] <- beta_from_r2(r2_req, m)
    }

    # Sanity check (should be ~0.80)
    row$power_check <- power_from_r2(N, alpha, r2_req)

    row
  }))
}

raw_tbl <- do.call(rbind, lapply(1:nrow(strata), function(i) {
  make_rows(strata$stratum[i], strata$N[i])
}))

## ---- Format for supplement (rounded, readable) ----
fmt_tbl <- raw_tbl
fmt_tbl$r2_required  <- round(fmt_tbl$r2_required, 4)   # R^2 shown as fraction (e.g., 0.0395)
fmt_tbl$power_check  <- round(fmt_tbl$power_check, 3)

beta_cols <- grep("^beta_MAF_", names(fmt_tbl), value = TRUE)
fmt_tbl[beta_cols] <- lapply(fmt_tbl[beta_cols], function(x) round(x, 3))

# Also provide percent variance explained for readability
fmt_tbl$r2_percent <- round(100 * raw_tbl$r2_required, 2)

# Reorder columns nicely
fmt_tbl <- fmt_tbl[, c("stratum","N","threshold","alpha","power_target",
                       "r2_required","r2_percent", beta_cols, "power_check")]

print(fmt_tbl, row.names = FALSE)

## ---- Write CSV for supplement ----
write.csv(fmt_tbl, file = "MBFR_GWAS_power_table.csv", row.names = FALSE)

message("Wrote: MBFR_GWAS_power_table.csv")
