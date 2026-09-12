###############################################################################
## BIMS/STAT 5304 - Bonus Graphs for Presentation
##
## 1. Power curves: P(reject H0: beta_drug=0) vs sample size, by design
## 2. Imbalance over time: cumulative |n_T - n_C| within covariate cells
## 3. Bias-variance heatmap: MSE decomposition across Design x Model

###############################################################################

cat("=== Bonus graphs for presentation ===\n")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})

out_dir <- file.path(getwd(), "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 0. Helper functions (in case you run this standalone) -----------------

if (!exists("BETA_DRUG")) {
  BETA0 <- -0.3; BETA_DRUG <- 1.0
  BETA_AGE <- 3.0; BETA_SEV <- -3.0
  BETA_BMI <- -0.05; BETA_BASE <- -0.03
}

if (!exists("simulate_population")) {
  simulate_population <- function(n) {
    data.frame(
      id = 1:n,
      Age_Group = sample(c("Under 50","50-65","Over 65"), n, replace = TRUE,
                         prob = c(0.30, 0.40, 0.30)),
      Disease_Sev = sample(c("Mild","Moderate","Severe"), n, replace = TRUE,
                           prob = c(0.50, 0.30, 0.20)),
      BMI = rnorm(n, 28, 4), Baseline_SBP = rnorm(n, 145, 10),
      Noise_Var1 = rnorm(n, 0, 1), Noise_Var2 = rpois(n, 5)
    )
  }
  set_factors <- function(df) {
    df$Age_Group   <- factor(df$Age_Group,   levels=c("Under 50","50-65","Over 65"))
    df$Disease_Sev <- factor(df$Disease_Sev, levels=c("Mild","Moderate","Severe"))
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
      nT_a <- sum(trt[1:(i-1)]==1 & prev$Age_Group  ==a)
      nC_a <- sum(trt[1:(i-1)]==0 & prev$Age_Group  ==a)
      nT_s <- sum(trt[1:(i-1)]==1 & prev$Disease_Sev==s)
      nC_s <- sum(trt[1:(i-1)]==0 & prev$Disease_Sev==s)
      G_T <- abs((nT_a+1)-nC_a) + abs((nT_s+1)-nC_s)
      G_C <- abs(nT_a-(nC_a+1)) + abs(nT_s-(nC_s+1))
      trt[i] <- if (G_T < G_C) rbinom(1,1,p_bias)
      else if (G_T > G_C) rbinom(1,1,1-p_bias)
      else rbinom(1,1,0.5)
    }
    df$Drug <- trt; df$Y <- generate_Y(df, trt); df
  }
  design_rpw <- function(df, beta = 3) {
    n <- nrow(df); trt <- integer(n); Y <- integer(n); uT <- 1; uC <- 1
    for (i in 1:n) {
      a <- rbinom(1, 1, uT/(uT+uC)); trt[i] <- a
      Y[i] <- generate_Y(df[i,,drop=FALSE], a)
      if (Y[i]==1) { if (a==1) uT <- uT+beta else uC <- uC+beta
      } else        { if (a==1) uC <- uC+beta else uT <- uT+beta }
    }
    df$Drug <- trt; df$Y <- Y; df
  }
}

# ---- 1. Power curves -------------------------------------------------------
#
# For each (design, sample size N), run B replicates, fit adjusted model,
# count fraction with p < 0.05 for Drug coefficient.

cat("\n[1/3] Building power curves...\n")
NSIZES <- c(50, 100, 150, 200, 300)
B_POWER <- 300   # replicates per (design, N) - bump to 500+ if you have time

run_power <- function(design_fn, n, B = B_POWER) {
  reject <- 0; total <- 0
  for (b in 1:B) {
    pop <- set_factors(simulate_population(n))
    d <- design_fn(pop)
    fit <- tryCatch(
      glm(Y ~ Drug + Age_Group + Disease_Sev, data = d, family = binomial),
      error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit) || !"Drug" %in% rownames(summary(fit)$coefficients)) next
    p <- summary(fit)$coefficients["Drug", "Pr(>|z|)"]
    if (!is.na(p)) { reject <- reject + (p < 0.05); total <- total + 1 }
  }
  if (total == 0) return(NA_real_)
  reject / total
}

set.seed(5304)
power_grid <- expand.grid(N = NSIZES,
                          Design = c("Simple", "Pocock", "RPW"),
                          stringsAsFactors = FALSE)
power_grid$Power <- NA_real_
for (i in seq_len(nrow(power_grid))) {
  fn <- switch(power_grid$Design[i],
               Simple = design_simple, Pocock = design_pocock, RPW = design_rpw)
  power_grid$Power[i] <- run_power(fn, power_grid$N[i])
  cat(sprintf("  %s, N=%d: power = %.3f\n",
              power_grid$Design[i], power_grid$N[i], power_grid$Power[i]))
}

p_power <- ggplot(power_grid, aes(N, Power, color = Design)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(title    = "Power to detect the drug effect by design",
       subtitle = "Probability that p < 0.05 for Drug, adjusted model. Dashed line = 80%.",
       x = "Sample size N", y = "Power") +
  theme_minimal(base_size = 13)
ggsave(file.path(out_dir, "plot_power_curves.png"), p_power,
       width = 7, height = 4.5, dpi = 200)
write.csv(power_grid, file.path(out_dir, "table_power_curves.csv"),
          row.names = FALSE)
cat("  saved: plot_power_curves.png + table_power_curves.csv\n")

# ---- 2. Imbalance over time ------------------------------------------------
#
# After each patient assignment, compute the cumulative |n_T - n_C|
# summed across all (Age_Group level, Disease_Sev level) cells.

cat("\n[2/3] Building imbalance-over-time plot...\n")

cumulative_imbalance <- function(d) {
  n <- nrow(d)
  imb <- numeric(n)
  for (i in 1:n) {
    sub <- d[1:i, ]
    age_imb <- with(sub,
                    sum(abs(table(Age_Group, factor(Drug, levels = c(0,1)))[, "1"] -
                              table(Age_Group, factor(Drug, levels = c(0,1)))[, "0"])))
    sev_imb <- with(sub,
                    sum(abs(table(Disease_Sev, factor(Drug, levels = c(0,1)))[, "1"] -
                              table(Disease_Sev, factor(Drug, levels = c(0,1)))[, "0"])))
    imb[i] <- age_imb + sev_imb
  }
  imb
}

# Average across replicates so it shows the systematic pattern, not one run
NREP_IMB <- 100
set.seed(5304)
imb_acc <- list(Simple = numeric(200), Pocock = numeric(200), RPW = numeric(200))
for (b in 1:NREP_IMB) {
  if (b %% 20 == 0) cat("  imbalance rep", b, "/", NREP_IMB, "\n")
  pop <- set_factors(simulate_population(200))
  imb_acc$Simple <- imb_acc$Simple + cumulative_imbalance(design_simple(pop))
  imb_acc$Pocock <- imb_acc$Pocock + cumulative_imbalance(design_pocock(pop))
  imb_acc$RPW    <- imb_acc$RPW    + cumulative_imbalance(design_rpw(pop))
}
imb_df <- bind_rows(lapply(names(imb_acc), function(nm) {
  data.frame(Design = nm, Patient = 1:200,
             MeanImbalance = imb_acc[[nm]] / NREP_IMB)
}))

p_imb <- ggplot(imb_df, aes(Patient, MeanImbalance, color = Design)) +
  geom_line(linewidth = 1) +
  labs(title    = "Average covariate imbalance during enrollment",
       subtitle = paste0("Mean across ", NREP_IMB, " replicates. Lower = better balance."),
       x = "Patient enrollment order",
       y = "Cumulative |n_Treatment - n_Control| across cells") +
  theme_minimal(base_size = 13)
ggsave(file.path(out_dir, "plot_imbalance_over_time.png"), p_imb,
       width = 7, height = 4.5, dpi = 200)
cat("  saved: plot_imbalance_over_time.png\n")

# ---- 3. Bias-variance heatmap (MSE decomposition) --------------------------
#
# Uses results from multireplicate_simulation.R if it has been run.
# If not, run a quick 500-rep simulation here to compute bias and variance.

cat("\n[3/3] Building bias-variance heatmap...\n")

if (!exists("results")) {
  cat("  multireplicate results not in memory; running 500 quick reps...\n")
  NREP_BV <- 500
  designs_v <- c("Simple", "Pocock", "RPW")
  fit_models <- function(d) {
    out <- list()
    for (m in c("naive","adj")) {
      f <- if (m == "naive") Y ~ Drug
      else Y ~ Drug + Age_Group + Disease_Sev
      fit <- tryCatch(glm(f, data = d, family = binomial),
                      error = function(e) NULL)
      if (is.null(fit) || !"Drug" %in% rownames(summary(fit)$coefficients)) {
        out[[m]] <- NA_real_; next
      }
      out[[m]] <- summary(fit)$coefficients["Drug","Estimate"]
    }
    out
  }
  results <- list()
  for (dn in designs_v) results[[dn]] <- data.frame(
    naive_beta = numeric(NREP_BV), adj_beta = numeric(NREP_BV))
  set.seed(5304)
  for (i in 1:NREP_BV) {
    if (i %% 100 == 0) cat("    bv rep", i, "/", NREP_BV, "\n")
    pop <- set_factors(simulate_population(200))
    for (dn in designs_v) {
      d <- switch(dn, Simple = design_simple(pop),
                  Pocock = design_pocock(pop),
                  RPW    = design_rpw(pop))
      fits <- fit_models(d)
      results[[dn]]$naive_beta[i] <- fits$naive
      results[[dn]]$adj_beta[i]   <- fits$adj
    }
  }
}

bv_summary <- bind_rows(lapply(c("Simple","Pocock","RPW"), function(dn) {
  df <- results[[dn]]
  data.frame(
    Design = dn,
    Model  = c("Naive", "Adjusted"),
    bias   = c(mean(df$naive_beta, na.rm = TRUE) - BETA_DRUG,
               mean(df$adj_beta,   na.rm = TRUE) - BETA_DRUG),
    variance = c(var(df$naive_beta, na.rm = TRUE),
                 var(df$adj_beta,   na.rm = TRUE))
  )
})) %>%
  mutate(bias_sq = bias^2, MSE = bias_sq + variance,
         Design = factor(Design, levels = c("Simple","Pocock","RPW")),
         Model  = factor(Model,  levels = c("Naive","Adjusted")))

cat("\n  Bias-variance summary:\n")
print(bv_summary, digits = 3, row.names = FALSE)
write.csv(bv_summary, file.path(out_dir, "table_bias_variance.csv"),
          row.names = FALSE)

p_bv <- ggplot(bv_summary, aes(Design, Model, fill = MSE)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("Bias\u00b2 = %.3f\nVar = %.3f\nMSE = %.3f",
                                bias_sq, variance, MSE)),
            color = "white", size = 4, lineheight = 1.05) +
  scale_fill_gradient(low = "#2c7fb8", high = "#d7301f",
                      name = "MSE", labels = scales::number_format(accuracy=0.01)) +
  labs(title = "MSE decomposition: bias\u00b2 + variance",
       subtitle = "Lower MSE (blue) is better. Adjusted models trade variance to remove bias.",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())
ggsave(file.path(out_dir, "plot_bias_variance_heatmap.png"), p_bv,
       width = 7.5, height = 4.5, dpi = 200)
cat("  saved: plot_bias_variance_heatmap.png + table_bias_variance.csv\n")

cat("\n=== All three bonus graphs complete ===\n")
cat("Files in", out_dir, ":\n")
print(list.files(out_dir, pattern = "plot_(power|imbalance|bias)|table_(power|bias)"))