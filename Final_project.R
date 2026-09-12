###############################################################################
## BIMS/STAT 5304 Final Project
## Comparing Simple, Pocock-Simon, and Play-the-Winner Randomization
## via simulation with a logistic-regression DGP
###############################################################################

# ---- 0. Setup --------------------------------------------------------------

set.seed(5304)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(knitr)
  library(MASS)     # stepAIC
  library(glmnet)   # LASSO
  library(broom)    # tidy model output
})

# True parameters (the DGP we are trying to recover)
BETA0     <- -0.3
BETA_DRUG <-  1.0
BETA_AGE  <-  3.0    # for "Over 65"
BETA_SEV  <- -3.0    # for "Severe"
BETA_BMI  <- -0.05
BETA_BASE <- -0.03

N <- 200             # trial size

# ---- 1. Task 1.1: Simulate the patient population --------------------------

simulate_population <- function(n = N) {
  data.frame(
    id           = 1:n,
    Age_Group    = sample(c("Under 50", "50-65", "Over 65"),
                          n, replace = TRUE, prob = c(0.30, 0.40, 0.30)),
    Disease_Sev  = sample(c("Mild", "Moderate", "Severe"),
                          n, replace = TRUE, prob = c(0.50, 0.30, 0.20)),
    BMI          = rnorm(n, mean = 28,  sd = 4),
    Baseline_SBP = rnorm(n, mean = 145, sd = 10),
    Noise_Var1   = rnorm(n, 0, 1),
    Noise_Var2   = rpois(n, lambda = 5)
  )
}

# Make factor levels explicit so reference categories are stable
set_factors <- function(df) {
  df$Age_Group   <- factor(df$Age_Group,
                           levels = c("Under 50", "50-65", "Over 65"))
  df$Disease_Sev <- factor(df$Disease_Sev,
                           levels = c("Mild", "Moderate", "Severe"))
  df
}

# Generate Y from the true DGP given covariates and treatment
generate_Y <- function(df, treatment) {
  bmi_c  <- df$BMI - 28
  sbp_c  <- df$Baseline_SBP - 145
  eta <- BETA0 +
    BETA_DRUG * treatment +
    BETA_AGE  * (df$Age_Group   == "Over 65") +
    BETA_SEV  * (df$Disease_Sev == "Severe")  +
    BETA_BMI  * bmi_c +
    BETA_BASE * sbp_c
  p <- 1 / (1 + exp(-eta))
  rbinom(length(p), size = 1, prob = p)
}

pop <- set_factors(simulate_population(N))


# ---- 2. Task 1.2: Design 1 - Simple Randomization --------------------------

design_simple <- function(df) {
  trt <- rbinom(nrow(df), size = 1, prob = 0.5)
  df$Drug <- trt
  df$Y    <- generate_Y(df, trt)
  df
}


# ---- 3. Task 1.3: Design 2 - Pocock & Simon Minimization -------------------


design_pocock <- function(df, p_bias = 0.75) {
  n   <- nrow(df)
  trt <- integer(n)
  
  # First patient: 50/50
  trt[1] <- rbinom(1, 1, 0.5)
  
  for (i in 2:n) {
    # counts so far in each arm by factor level
    prev <- df[1:(i - 1), ]
    a_lvl <- df$Age_Group[i]
    s_lvl <- df$Disease_Sev[i]
    
    nT_age <- sum(trt[1:(i - 1)] == 1 & prev$Age_Group   == a_lvl)
    nC_age <- sum(trt[1:(i - 1)] == 0 & prev$Age_Group   == a_lvl)
    nT_sev <- sum(trt[1:(i - 1)] == 1 & prev$Disease_Sev == s_lvl)
    nC_sev <- sum(trt[1:(i - 1)] == 0 & prev$Disease_Sev == s_lvl)
    
    # imbalance if we assign to T  vs  if we assign to C
    G_if_T <- abs((nT_age + 1) - nC_age) + abs((nT_sev + 1) - nC_sev)
    G_if_C <- abs(nT_age - (nC_age + 1)) + abs(nT_sev - (nC_sev + 1))
    
    if (G_if_T < G_if_C) {
      trt[i] <- rbinom(1, 1, p_bias)            # prefer T
    } else if (G_if_T > G_if_C) {
      trt[i] <- rbinom(1, 1, 1 - p_bias)        # prefer C
    } else {
      trt[i] <- rbinom(1, 1, 0.5)               # tie
    }
  }
  
  df$Drug <- trt
  df$Y    <- generate_Y(df, trt)
  df
}


# ---- 4. Task 1.4: Design 3 - Randomized Play-the-Winner --------------------


design_rpw <- function(df, beta = 3) {
  n   <- nrow(df)
  trt <- integer(n)
  Y   <- integer(n)
  urnT <- 1
  urnC <- 1
  
  for (i in 1:n) {
    pT <- urnT / (urnT + urnC)
    a  <- rbinom(1, 1, pT)        # 1 = T, 0 = C
    trt[i] <- a
    
    # observe outcome for THIS patient (one row)
    Y[i] <- generate_Y(df[i, , drop = FALSE], a)
    
    # update urn
    if (Y[i] == 1) {              # success -> add same color
      if (a == 1) urnT <- urnT + beta else urnC <- urnC + beta
    } else {                      # failure -> add opposite color
      if (a == 1) urnC <- urnC + beta else urnT <- urnT + beta
    }
  }
  
  df$Drug <- trt
  df$Y    <- Y
  df
}


# ---- 5. Run the three designs on the SAME population -----------------------

dat_simple <- design_simple(pop)
dat_pocock <- design_pocock(pop)
dat_rpw    <- design_rpw(pop)

datasets <- list(Simple = dat_simple,
                 Pocock = dat_pocock,
                 RPW    = dat_rpw)


# ---- 6. Task 1.5: Allocation summary + cumulative-proportion plot ----------

allocation_summary <- sapply(datasets, function(d) {
  c(N_Treatment = sum(d$Drug == 1),
    N_Control   = sum(d$Drug == 0),
    Pct_Drug    = round(100 * mean(d$Drug == 1), 1))
})
cat("\n=== Allocation summary ===\n")
print(allocation_summary)

alloc_df <- bind_rows(lapply(names(datasets), function(nm) {
  d <- datasets[[nm]]
  data.frame(
    Design = nm,
    Patient = seq_len(nrow(d)),
    CumPropDrug = cumsum(d$Drug) / seq_len(nrow(d))
  )
}))

p_alloc <- ggplot(alloc_df, aes(Patient, CumPropDrug, color = Design)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9) +
  labs(title    = "Cumulative treatment allocation by design",
       y        = "Cumulative proportion assigned to DrugX",
       x        = "Patient enrollment order") +
  theme_minimal(base_size = 12)
ggsave("plot_allocation.png", p_alloc, width = 7, height = 4.5, dpi = 200)


# ---- 7. Task 2: EDA --------------------------------------------------------

eda_balance <- function(d, label) {
  cat("\n--- Balance tables:", label, "---\n")
  cat("Age Group x Treatment:\n")
  print(table(d$Age_Group, d$Drug, dnn = c("Age", "Drug")))
  cat("Disease Severity x Treatment:\n")
  print(table(d$Disease_Sev, d$Drug, dnn = c("Severity", "Drug")))
  
  # Continuous covariate t-tests
  bmi_t <- t.test(BMI ~ Drug, data = d)
  sbp_t <- t.test(Baseline_SBP ~ Drug, data = d)
  cat(sprintf("BMI: mean(C)=%.2f mean(T)=%.2f  t=%.2f  p=%.3f\n",
              bmi_t$estimate[1], bmi_t$estimate[2],
              bmi_t$statistic, bmi_t$p.value))
  cat(sprintf("SBP: mean(C)=%.2f mean(T)=%.2f  t=%.2f  p=%.3f\n",
              sbp_t$estimate[1], sbp_t$estimate[2],
              sbp_t$statistic, sbp_t$p.value))
  
  # Observed success rates
  rates <- d %>%
    group_by(Drug) %>%
    summarise(n = n(), successes = sum(Y), rate = mean(Y), .groups = "drop")
  cat("Observed success rates:\n"); print(rates)
  
  invisible(list(bmi_t = bmi_t, sbp_t = sbp_t, rates = rates))
}

cat("\n========== TASK 2: EDA ==========\n")
eda_results <- lapply(names(datasets),
                      function(nm) eda_balance(datasets[[nm]], nm))
names(eda_results) <- names(datasets)

# Combined visualization: success rates by Drug across designs
rate_df <- bind_rows(lapply(names(datasets), function(nm) {
  d <- datasets[[nm]]
  d %>% group_by(Drug) %>%
    summarise(rate = mean(Y), n = n(), .groups = "drop") %>%
    mutate(Design = nm,
           Arm = ifelse(Drug == 1, "DrugX", "Control"))
}))

p_rates <- ggplot(rate_df, aes(Design, rate, fill = Arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * rate)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Observed success rate by arm and design",
       y = "Success rate", x = NULL) +
  theme_minimal(base_size = 12)
ggsave("plot_success_rates.png", p_rates, width = 7, height = 4.5, dpi = 200)


# ---- 8. Task 3: Modeling ---------------------------------------------------

extract_drug <- function(fit, label) {
  s   <- summary(fit)$coefficients
  est <- s["Drug", "Estimate"]
  se  <- s["Drug", "Std. Error"]
  p   <- s["Drug", "Pr(>|z|)"]
  ci  <- suppressMessages(confint(fit, "Drug"))
  data.frame(
    Model    = label,
    beta_hat = est,
    SE       = se,
    p_value  = p,
    CI_low   = ci[1],
    CI_high  = ci[2],
    OR       = exp(est),
    OR_low   = exp(ci[1]),
    OR_high  = exp(ci[2])
  )
}

cat("\n========== TASK 3.1: NAIVE MODEL (Y ~ Drug) ==========\n")
naive_results <- bind_rows(lapply(names(datasets), function(nm) {
  fit <- glm(Y ~ Drug, data = datasets[[nm]], family = binomial)
  cbind(Design = nm, extract_drug(fit, "Naive"))
}))
print(naive_results, row.names = FALSE)

cat("\n========== TASK 3.2: ADJUSTED MODEL ==========\n")
adj_results <- bind_rows(lapply(names(datasets), function(nm) {
  fit <- glm(Y ~ Drug + Age_Group + Disease_Sev,
             data = datasets[[nm]], family = binomial)
  cbind(Design = nm, extract_drug(fit, "Adjusted"))
}))
print(adj_results, row.names = FALSE)

# Combined coefficient table
coef_table <- bind_rows(naive_results, adj_results) %>%
  arrange(Design, Model)
write.csv(coef_table, "table_drug_coefs.csv", row.names = FALSE)

# Forest-style plot of beta_drug across designs / models
p_forest <- ggplot(coef_table,
                   aes(x = beta_hat, y = paste(Design, Model, sep = " - "),
                       color = Model)) +
  geom_vline(xintercept = BETA_DRUG, linetype = "dashed", color = "red") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2) +
  labs(title = "Drug coefficient estimates (95% CI). Red = true value (1.0)",
       x = expression(hat(beta)[drug]), y = NULL) +
  theme_minimal(base_size = 12)
ggsave("plot_forest.png", p_forest, width = 8, height = 4.5, dpi = 200)


# ---- 9. Task 3.3: Model selection on the Pocock dataset --------------------

cat("\n========== TASK 3.3: MODEL SELECTION (Pocock dataset) ==========\n")

full_fit <- glm(Y ~ Drug + Age_Group + Disease_Sev +
                  BMI + Baseline_SBP + Noise_Var1 + Noise_Var2,
                data = dat_pocock, family = binomial)

cat("\n-- Full model summary --\n")
print(summary(full_fit))

# Stepwise AIC
step_fit <- stepAIC(full_fit, direction = "both", trace = FALSE)
cat("\n-- Stepwise AIC selected model --\n")
print(formula(step_fit))
print(summary(step_fit)$coefficients)

# LASSO via glmnet
X <- model.matrix(Y ~ Drug + Age_Group + Disease_Sev +
                    BMI + Baseline_SBP + Noise_Var1 + Noise_Var2,
                  data = dat_pocock)[, -1]
y <- dat_pocock$Y

set.seed(5304)
cv_lasso <- cv.glmnet(X, y, family = "binomial", alpha = 1)

cat("\n-- LASSO coefficients at lambda.min --\n")
print(coef(cv_lasso, s = "lambda.min"))
cat("\n-- LASSO coefficients at lambda.1se --\n")
print(coef(cv_lasso, s = "lambda.1se"))

# Save coefficient paths plot
png("plot_lasso_path.png", width = 1400, height = 900, res = 200)
plot(cv_lasso)
dev.off()

# Compare beta_drug across selection methods
sel_summary <- data.frame(
  Method   = c("Full", "Stepwise AIC", "LASSO (lambda.min)", "LASSO (lambda.1se)"),
  beta_drug = c(coef(full_fit)["Drug"],
                coef(step_fit)["Drug"],
                as.numeric(coef(cv_lasso, s = "lambda.min")["Drug", ]),
                as.numeric(coef(cv_lasso, s = "lambda.1se")["Drug", ]))
)
cat("\n-- Drug coefficient across selection methods (Pocock dataset) --\n")
print(sel_summary, row.names = FALSE)
write.csv(sel_summary, "table_model_selection.csv", row.names = FALSE)


# ---- 10. Save datasets for the appendix ------------------------------------

write.csv(dat_simple, "dataset_simple.csv", row.names = FALSE)
write.csv(dat_pocock, "dataset_pocock.csv", row.names = FALSE)
write.csv(dat_rpw,    "dataset_rpw.csv",    row.names = FALSE)

cat("\nAll outputs written. Files: plot_allocation.png, plot_success_rates.png,\n",
    "plot_forest.png, plot_lasso_path.png, table_drug_coefs.csv,\n",
    "table_model_selection.csv, dataset_*.csv\n")