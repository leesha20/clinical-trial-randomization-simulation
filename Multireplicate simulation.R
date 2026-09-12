###############################################################################
## BIMS/STAT 5304 Final Project - MULTI-REPLICATE SIMULATION

###############################################################################

cat("=== Multi-replicate simulation ===\n")
cat("Working directory:", getwd(), "\n\n")

# ---- Setup -----------------------------------------------------------------

needed <- c("ggplot2", "dplyr", "tidyr", "scales", "knitr")
for (p in needed) if (!p %in% rownames(installed.packages()))
  install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales); library(knitr)
})

out_dir <- file.path(getwd(), "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Constants and core functions (same DGP as single-run script) ----------

BETA0 <- -0.3; BETA_DRUG <- 1.0
BETA_AGE <- 3.0; BETA_SEV <- -3.0
BETA_BMI <- -0.05; BETA_BASE <- -0.03
N <- 200

simulate_population <- function(n = N) {
  data.frame(
    id           = 1:n,
    Age_Group    = sample(c("Under 50","50-65","Over 65"), n, replace = TRUE,
                          prob = c(0.30, 0.40, 0.30)),
    Disease_Sev  = sample(c("Mild","Moderate","Severe"), n, replace = TRUE,
                          prob = c(0.50, 0.30, 0.20)),
    BMI          = rnorm(n, 28, 4),
    Baseline_SBP = rnorm(n, 145, 10),
    Noise_Var1   = rnorm(n, 0, 1),
    Noise_Var2   = rpois(n, 5)
  )
}
set_factors <- function(df) {
  df$Age_Group   <- factor(df$Age_Group,   levels = c("Under 50","50-65","Over 65"))
  df$Disease_Sev <- factor(df$Disease_Sev, levels = c("Mild","Moderate","Severe"))
  df
}
generate_Y <- function(df, treatment) {
  bmi_c <- df$BMI - 28; sbp_c <- df$Baseline_SBP - 145
  eta <- BETA0 + BETA_DRUG * treatment +
    BETA_AGE * (df$Age_Group   == "Over 65") +
    BETA_SEV * (df$Disease_Sev == "Severe") +
    BETA_BMI * bmi_c + BETA_BASE * sbp_c
  rbinom(length(eta), 1, 1 / (1 + exp(-eta)))
}

design_simple <- function(df) {
  trt <- rbinom(nrow(df), 1, 0.5)
  df$Drug <- trt; df$Y <- generate_Y(df, trt); df
}
design_pocock <- function(df, p_bias = 0.75) {
  n <- nrow(df); trt <- integer(n); trt[1] <- rbinom(1, 1, 0.5)
  for (i in 2:n) {
    prev <- df[1:(i-1), ]
    a <- df$Age_Group[i]; s <- df$Disease_Sev[i]
    nT_a <- sum(trt[1:(i-1)] == 1 & prev$Age_Group   == a)
    nC_a <- sum(trt[1:(i-1)] == 0 & prev$Age_Group   == a)
    nT_s <- sum(trt[1:(i-1)] == 1 & prev$Disease_Sev == s)
    nC_s <- sum(trt[1:(i-1)] == 0 & prev$Disease_Sev == s)
    G_T <- abs((nT_a + 1) - nC_a) + abs((nT_s + 1) - nC_s)
    G_C <- abs(nT_a - (nC_a + 1)) + abs(nT_s - (nC_s + 1))
    trt[i] <- if (G_T < G_C)      rbinom(1, 1, p_bias)
    else if (G_T > G_C) rbinom(1, 1, 1 - p_bias)
    else                rbinom(1, 1, 0.5)
  }
  df$Drug <- trt; df$Y <- generate_Y(df, trt); df
}
design_rpw <- function(df, beta = 3) {
  n <- nrow(df); trt <- integer(n); Y <- integer(n); uT <- 1; uC <- 1
  for (i in 1:n) {
    a <- rbinom(1, 1, uT / (uT + uC)); trt[i] <- a
    Y[i] <- generate_Y(df[i, , drop = FALSE], a)
    if (Y[i] == 1) { if (a == 1) uT <- uT + beta else uC <- uC + beta
    } else         { if (a == 1) uC <- uC + beta else uT <- uT + beta }
  }
  df$Drug <- trt; df$Y <- Y; df
}

# ---- Helpers ---------------------------------------------------------------

# Marginal categorical imbalance: sum over (factor x level) of |n_T - n_C|
cat_imbalance <- function(d) {
  age <- table(d$Age_Group, d$Drug)
  sev <- table(d$Disease_Sev, d$Drug)
  sum(abs(age[, "1"] - age[, "0"])) + sum(abs(sev[, "1"] - sev[, "0"]))
}

# Fit naive + adjusted, return drug coefficient summaries
fit_models <- function(d) {
  out <- list(naive = NULL, adj = NULL)
  for (which in c("naive", "adj")) {
    f <- if (which == "naive") Y ~ Drug
    else Y ~ Drug + Age_Group + Disease_Sev
    fit <- tryCatch(glm(f, data = d, family = binomial),
                    error = function(e) NULL)
    if (is.null(fit) || !"Drug" %in% rownames(summary(fit)$coefficients)) {
      out[[which]] <- c(beta_hat = NA, SE = NA, CI_low = NA, CI_high = NA)
      next
    }
    s <- summary(fit)$coefficients
    est <- s["Drug", "Estimate"]; se <- s["Drug", "Std. Error"]
    out[[which]] <- c(beta_hat = est, SE = se,
                      CI_low = est - 1.96 * se, CI_high = est + 1.96 * se)
  }
  out
}

# Run one full replicate: same population, three designs, naive + adjusted
one_replicate <- function() {
  pop <- set_factors(simulate_population(N))
  res <- list()
  for (design_name in c("Simple", "Pocock", "RPW")) {
    d <- switch(design_name,
                Simple = design_simple(pop),
                Pocock = design_pocock(pop),
                RPW    = design_rpw(pop))
    fits <- fit_models(d)
    res[[design_name]] <- list(
      pct_drug   = mean(d$Drug == 1),
      n_drug     = sum(d$Drug == 1),
      imbalance  = cat_imbalance(d),
      naive      = fits$naive,
      adj        = fits$adj
    )
  }
  res
}

# ---- Run the replicates ----------------------------------------------------

NREP <- 1000
set.seed(5304)
cat("Running", NREP, "replicates... (this takes a couple of minutes)\n")

# Pre-allocate result containers
designs <- c("Simple", "Pocock", "RPW")
results <- list()
for (dn in designs) {
  results[[dn]] <- data.frame(
    rep = 1:NREP, pct_drug = NA, imbalance = NA,
    naive_beta = NA, naive_SE = NA, naive_low = NA, naive_high = NA,
    adj_beta = NA, adj_SE = NA, adj_low = NA, adj_high = NA
  )
}

t0 <- Sys.time()
for (i in 1:NREP) {
  if (i %% 100 == 0) cat("  rep", i, "/", NREP, "\n")
  r <- one_replicate()
  for (dn in designs) {
    results[[dn]][i, "pct_drug"]   <- r[[dn]]$pct_drug
    results[[dn]][i, "imbalance"]  <- r[[dn]]$imbalance
    results[[dn]][i, "naive_beta"] <- r[[dn]]$naive["beta_hat"]
    results[[dn]][i, "naive_SE"]   <- r[[dn]]$naive["SE"]
    results[[dn]][i, "naive_low"]  <- r[[dn]]$naive["CI_low"]
    results[[dn]][i, "naive_high"] <- r[[dn]]$naive["CI_high"]
    results[[dn]][i, "adj_beta"]   <- r[[dn]]$adj["beta_hat"]
    results[[dn]][i, "adj_SE"]     <- r[[dn]]$adj["SE"]
    results[[dn]][i, "adj_low"]    <- r[[dn]]$adj["CI_low"]
    results[[dn]][i, "adj_high"]   <- r[[dn]]$adj["CI_high"]
  }
}
cat("Done in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "seconds.\n\n")

# ---- Summarize across replicates -------------------------------------------

summarize_design <- function(dn) {
  df <- results[[dn]]
  data.frame(
    Design = dn,
    pct_drug_mean = mean(df$pct_drug),
    pct_drug_sd   = sd(df$pct_drug),
    imbalance_mean = mean(df$imbalance),
    naive_mean = mean(df$naive_beta, na.rm = TRUE),
    naive_emp_SE = sd(df$naive_beta, na.rm = TRUE),
    naive_bias = mean(df$naive_beta, na.rm = TRUE) - BETA_DRUG,
    naive_coverage = mean(df$naive_low <= BETA_DRUG & df$naive_high >= BETA_DRUG,
                          na.rm = TRUE),
    adj_mean = mean(df$adj_beta, na.rm = TRUE),
    adj_emp_SE = sd(df$adj_beta, na.rm = TRUE),
    adj_bias = mean(df$adj_beta, na.rm = TRUE) - BETA_DRUG,
    adj_coverage = mean(df$adj_low <= BETA_DRUG & df$adj_high >= BETA_DRUG,
                        na.rm = TRUE)
  )
}
summary_tab <- bind_rows(lapply(designs, summarize_design))

cat("=== Summary across", NREP, "replicates ===\n")
cat("(true beta_drug = 1.0)\n\n")
print(summary_tab, row.names = FALSE, digits = 3)

write.csv(summary_tab, file.path(out_dir, "table_multireplicate_summary.csv"),
          row.names = FALSE)
cat("\nSaved:", file.path(out_dir, "table_multireplicate_summary.csv"), "\n\n")

# ---- Plots -----------------------------------------------------------------

# Long-format estimate frame
est_long <- bind_rows(lapply(designs, function(dn) {
  df <- results[[dn]]
  bind_rows(
    data.frame(Design = dn, Model = "Naive",    beta = df$naive_beta),
    data.frame(Design = dn, Model = "Adjusted", beta = df$adj_beta)
  )
})) %>%
  mutate(Model = factor(Model, levels = c("Naive", "Adjusted")),
         Design = factor(Design, levels = c("Simple", "Pocock", "RPW")))

# Plot A: density of beta_hat by design x model
pA <- ggplot(est_long, aes(x = beta, fill = Model)) +
  geom_density(alpha = 0.55, color = NA) +
  geom_vline(xintercept = BETA_DRUG, linetype = "dashed", color = "red") +
  facet_wrap(~ Design, ncol = 1) +
  coord_cartesian(xlim = c(-0.5, 3.5)) +
  labs(title = paste0("Sampling distribution of beta_drug across ", NREP, " replicates"),
       subtitle = "Red dashed line = true value (1.0)",
       x = expression(hat(beta)[drug]), y = "Density") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "plot_mc_distributions.png"), pA,
       width = 7.5, height = 6, dpi = 200)
cat("Saved: plot_mc_distributions.png\n")

# Plot B: mean estimate +/- empirical SE (the "clean" forest plot)
fp_summary <- est_long %>%
  group_by(Design, Model) %>%
  summarise(mean_beta = mean(beta, na.rm = TRUE),
            emp_SE = sd(beta, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(low = mean_beta - 1.96 * emp_SE,
         high = mean_beta + 1.96 * emp_SE,
         label = paste(Design, Model, sep = " - "))

pB <- ggplot(fp_summary, aes(mean_beta, label, color = Model)) +
  geom_vline(xintercept = BETA_DRUG, linetype = "dashed", color = "red") +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y",
                width = 0.2, linewidth = 0.8) +
  labs(title = paste0("Mean beta_drug +/- 1.96 x empirical SE (", NREP, " reps)"),
       subtitle = "Red dashed line = true value (1.0). This is the *systematic* behavior of each design.",
       x = expression("Mean " ~ hat(beta)[drug]), y = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "plot_mc_forest.png"), pB,
       width = 8, height = 4.5, dpi = 200)
cat("Saved: plot_mc_forest.png\n")

# Plot C: distribution of % allocated to Drug
alloc_long <- bind_rows(lapply(designs, function(dn)
  data.frame(Design = dn, pct_drug = results[[dn]]$pct_drug))) %>%
  mutate(Design = factor(Design, levels = c("Simple", "Pocock", "RPW")))

pC <- ggplot(alloc_long, aes(pct_drug, fill = Design)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  facet_wrap(~ Design, ncol = 1, scales = "free_y") +
  scale_x_continuous(labels = percent, limits = c(0.3, 0.85)) +
  labs(title = paste0("Distribution of final % on Drug across ", NREP, " replicates"),
       subtitle = "Pocock concentrates near 50%; RPW shifts right (ethical allocation); Simple varies",
       x = "Final proportion on DrugX", y = "Frequency") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "plot_mc_allocation.png"), pC,
       width = 7.5, height = 6, dpi = 200)
cat("Saved: plot_mc_allocation.png\n")

# Plot D: distribution of categorical imbalance score
imb_long <- bind_rows(lapply(designs, function(dn)
  data.frame(Design = dn, imbalance = results[[dn]]$imbalance))) %>%
  mutate(Design = factor(Design, levels = c("Simple", "Pocock", "RPW")))

pD <- ggplot(imb_long, aes(imbalance, fill = Design)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ Design, ncol = 1, scales = "free_y") +
  labs(title = "Categorical imbalance score (lower = better balance)",
       subtitle = "Sum over (factor x level) of |n_Treatment - n_Control|",
       x = "Total imbalance across Age_Group + Disease_Sev levels",
       y = "Frequency") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "plot_mc_imbalance.png"), pD,
       width = 7.5, height = 6, dpi = 200)
cat("Saved: plot_mc_imbalance.png\n")

cat("\n=== Multi-replicate run complete ===\n")
cat("All outputs in:", out_dir, "\n")