###############################################################################
## BIMS/STAT 5304 Final Project - PRESENTATION TABLES

## Outputs to output/:
##   - table_1_design_overview.csv      | the master comparison table
##   - table_2_anchor_probs.csv         | DGP anchor probabilities
##   - table_3_balance_pocock.csv       | balance tables for Pocock
##   - table_3_balance_simple.csv       | balance tables for Simple
##   - table_3_balance_rpw.csv          | balance tables for RPW
##   - table_4_continuous_balance.csv   | BMI / SBP balance + t-tests
##   - table_5_success_rates.csv        | observed success rates by arm
##   - table_6_drug_coefs.csv           | naive vs adjusted across designs
##   - table_7_model_selection.csv      | LASSO + stepwise on Pocock
##   - table_8_mc_summary.csv           | multi-replicate summary
##   - table_9_when_to_use.csv          | decision matrix for conclusion

###############################################################################

cat("=== Building presentation tables ===\n\n")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(knitr)
})

out_dir <- file.path(getwd(), "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Helper: pretty print a data frame to the console with knitr::kable
show <- function(df, caption) {
  cat("\n###", caption, "###\n")
  print(kable(df, format = "simple", digits = 3))
  cat("\n")
}

# ---- Sanity checks ---------------------------------------------------------

needed_objects <- c("dat_simple", "dat_pocock", "dat_rpw")
missing_objs <- needed_objects[!sapply(needed_objects, exists)]
if (length(missing_objs) > 0) {
  stop("Missing objects: ", paste(missing_objs, collapse = ", "),
       "\nRun final_project.R first, then source this script.")
}

have_mc <- exists("results") &&
  all(c("Simple","Pocock","RPW") %in% names(results))
if (!have_mc) {
  cat("NOTE: multireplicate results not in memory.\n")
  cat("      Skipping Table 8. Run multireplicate_simulation.R first to include it.\n\n")
}

datasets <- list(Simple = dat_simple, Pocock = dat_pocock, RPW = dat_rpw)

# ---- Table 1: Design overview (the master table) ---------------------------

t1 <- data.frame(
  Design = c("Simple", "Pocock", "RPW"),
  N_Drug = sapply(datasets, function(d) sum(d$Drug == 1)),
  N_Control = sapply(datasets, function(d) sum(d$Drug == 0)),
  Pct_on_Drug = sapply(datasets,
                       function(d) sprintf("%.1f%%", 100 * mean(d$Drug == 1))),
  Drug_Success_Rate = sapply(datasets,
                             function(d) sprintf("%.1f%%",
                                                 100 * mean(d$Y[d$Drug == 1]))),
  Control_Success_Rate = sapply(datasets,
                                function(d) sprintf("%.1f%%",
                                                    100 * mean(d$Y[d$Drug == 0]))),
  row.names = NULL
)
show(t1, "Table 1: Design overview")
write.csv(t1, file.path(out_dir, "table_1_design_overview.csv"), row.names = FALSE)

# ---- Table 2: Anchor probabilities from the DGP ----------------------------

# Reproduce the four anchor probabilities the assignment cites
inv_logit <- function(x) 1 / (1 + exp(-x))
B0 <- -0.3; BD <- 1.0; BA <- 3.0; BS <- -3.0
t2 <- data.frame(
  Patient_Profile = c("Younger, non-severe",
                      "Younger, non-severe",
                      "Over 65, non-severe",
                      "Younger, Severe"),
  Treatment = c("Placebo", "DrugX", "Placebo", "Placebo"),
  Linear_Predictor = c(B0,
                       B0 + BD,
                       B0 + BA,
                       B0 + BS),
  Pr_Success = sapply(c(B0, B0 + BD, B0 + BA, B0 + BS),
                      function(x) sprintf("%.1f%%", 100 * inv_logit(x)))
)
show(t2, "Table 2: DGP anchor probabilities")
write.csv(t2, file.path(out_dir, "table_2_anchor_probs.csv"), row.names = FALSE)

# ---- Table 3: Categorical balance tables -----------------------------------

balance_table <- function(d, factor_col, label) {
  tab <- table(d[[factor_col]], d$Drug)
  out <- as.data.frame.matrix(tab)
  names(out) <- c("Control", "Treatment")
  out$Total <- out$Control + out$Treatment
  out$Imbalance <- abs(out$Control - out$Treatment)
  out <- cbind(Level = rownames(out), out)
  rownames(out) <- NULL
  out
}

for (nm in names(datasets)) {
  d <- datasets[[nm]]
  age_tab <- balance_table(d, "Age_Group", nm)
  age_tab$Factor <- "Age Group"
  sev_tab <- balance_table(d, "Disease_Sev", nm)
  sev_tab$Factor <- "Disease Severity"
  combined <- rbind(age_tab, sev_tab)
  combined <- combined[, c("Factor", "Level", "Control", "Treatment", "Total",
                           "Imbalance")]
  show(combined, paste0("Table 3 (", nm, "): Categorical balance"))
  write.csv(combined,
            file.path(out_dir, paste0("table_3_balance_", tolower(nm), ".csv")),
            row.names = FALSE)
}

# ---- Table 4: Continuous covariate balance ---------------------------------

t4 <- bind_rows(lapply(names(datasets), function(nm) {
  d <- datasets[[nm]]
  bmi_t <- t.test(BMI ~ Drug, data = d)
  sbp_t <- t.test(Baseline_SBP ~ Drug, data = d)
  data.frame(
    Design = nm,
    Variable = c("BMI", "Baseline SBP"),
    Mean_Control = c(bmi_t$estimate[1], sbp_t$estimate[1]),
    Mean_Treatment = c(bmi_t$estimate[2], sbp_t$estimate[2]),
    Difference = c(diff(bmi_t$estimate), diff(sbp_t$estimate)),
    t_statistic = c(bmi_t$statistic, sbp_t$statistic),
    p_value = c(bmi_t$p.value, sbp_t$p.value)
  )
}))
show(t4, "Table 4: Continuous covariate balance (t-tests)")
write.csv(t4, file.path(out_dir, "table_4_continuous_balance.csv"),
          row.names = FALSE)

# ---- Table 5: Success rates by arm and design ------------------------------

t5 <- bind_rows(lapply(names(datasets), function(nm) {
  d <- datasets[[nm]]
  data.frame(
    Design = nm,
    Arm = c("Control", "DrugX"),
    N = c(sum(d$Drug == 0), sum(d$Drug == 1)),
    Successes = c(sum(d$Y[d$Drug == 0]), sum(d$Y[d$Drug == 1])),
    Failures = c(sum(d$Drug == 0) - sum(d$Y[d$Drug == 0]),
                 sum(d$Drug == 1) - sum(d$Y[d$Drug == 1])),
    Success_Rate = c(mean(d$Y[d$Drug == 0]), mean(d$Y[d$Drug == 1]))
  )
}))
show(t5, "Table 5: Observed success rates by arm")
write.csv(t5, file.path(out_dir, "table_5_success_rates.csv"), row.names = FALSE)

# ---- Table 6: Drug coefficients across naive and adjusted ------------------

extract_drug <- function(fit) {
  s <- summary(fit)$coefficients
  est <- s["Drug","Estimate"]; se <- s["Drug","Std. Error"]
  p <- s["Drug","Pr(>|z|)"]
  ci <- suppressMessages(confint(fit, "Drug"))
  data.frame(
    beta_hat = est, SE = se, p_value = p,
    CI_low = ci[1], CI_high = ci[2],
    OR = exp(est),
    OR_CI_low = exp(ci[1]), OR_CI_high = exp(ci[2])
  )
}

t6 <- bind_rows(lapply(names(datasets), function(nm) {
  d <- datasets[[nm]]
  naive_fit <- glm(Y ~ Drug, data = d, family = binomial)
  adj_fit   <- glm(Y ~ Drug + Age_Group + Disease_Sev, data = d, family = binomial)
  rbind(
    cbind(Design = nm, Model = "Naive",    extract_drug(naive_fit)),
    cbind(Design = nm, Model = "Adjusted", extract_drug(adj_fit))
  )
}))
t6 <- t6 %>% arrange(Design, Model)
show(t6, "Table 6: Drug coefficient estimates")
write.csv(t6, file.path(out_dir, "table_6_drug_coefs.csv"), row.names = FALSE)

# Pretty-printed version with formatted CI
t6_pretty <- t6 %>%
  mutate(Estimate_CI = sprintf("%.2f (%.2f, %.2f)", beta_hat, CI_low, CI_high),
         OR_CI = sprintf("%.2f (%.2f, %.2f)", OR, OR_CI_low, OR_CI_high),
         p = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))) %>%
  dplyr::select(Design, Model, `beta (95% CI)` = Estimate_CI, SE,
                `OR (95% CI)` = OR_CI, p)
show(t6_pretty, "Table 6 (pretty): Formatted for the report")
write.csv(t6_pretty, file.path(out_dir, "table_6_drug_coefs_pretty.csv"),
          row.names = FALSE)

# ---- Table 7: Model selection (Pocock dataset) -----------------------------

if (exists("sel_summary")) {
  show(sel_summary, "Table 7: Drug coefficient across selection methods")
  write.csv(sel_summary, file.path(out_dir, "table_7_model_selection.csv"),
            row.names = FALSE)
} else if (exists("step_fit") && exists("cv_lasso")) {
  full_fit <- glm(Y ~ Drug + Age_Group + Disease_Sev + BMI + Baseline_SBP +
                    Noise_Var1 + Noise_Var2,
                  data = dat_pocock, family = binomial)
  t7 <- data.frame(
    Method = c("Full", "Stepwise AIC", "LASSO (lambda.min)", "LASSO (lambda.1se)"),
    beta_drug = c(coef(full_fit)["Drug"], coef(step_fit)["Drug"],
                  as.numeric(coef(cv_lasso, s = "lambda.min")["Drug", ]),
                  as.numeric(coef(cv_lasso, s = "lambda.1se")["Drug", ]))
  )
  show(t7, "Table 7: Drug coefficient across selection methods")
  write.csv(t7, file.path(out_dir, "table_7_model_selection.csv"),
            row.names = FALSE)
} else {
  cat("Skipping Table 7 - run Task 3.3 from final_project.R first.\n\n")
}

# ---- Table 8: Multi-replicate summary --------------------------------------

if (have_mc) {
  TRUE_BETA <- 1.0
  t8 <- bind_rows(lapply(c("Simple","Pocock","RPW"), function(dn) {
    df <- results[[dn]]
    data.frame(
      Design = dn,
      Mean_pct_Drug = sprintf("%.1f%%", 100 * mean(df$pct_drug)),
      SD_pct_Drug = sprintf("%.3f", sd(df$pct_drug)),
      Naive_Mean = mean(df$naive_beta, na.rm = TRUE),
      Naive_Bias = mean(df$naive_beta, na.rm = TRUE) - TRUE_BETA,
      Naive_EmpSE = sd(df$naive_beta, na.rm = TRUE),
      Naive_Coverage = sprintf("%.1f%%",
                               100 * mean(df$naive_low <= TRUE_BETA &
                                            df$naive_high >= TRUE_BETA, na.rm = TRUE)),
      Adj_Mean = mean(df$adj_beta, na.rm = TRUE),
      Adj_Bias = mean(df$adj_beta, na.rm = TRUE) - TRUE_BETA,
      Adj_EmpSE = sd(df$adj_beta, na.rm = TRUE),
      Adj_Coverage = sprintf("%.1f%%",
                             100 * mean(df$adj_low <= TRUE_BETA &
                                          df$adj_high >= TRUE_BETA, na.rm = TRUE))
    )
  }))
  show(t8, "Table 8: Multi-replicate summary (true beta_drug = 1.0)")
  write.csv(t8, file.path(out_dir, "table_8_mc_summary.csv"), row.names = FALSE)
}

# ---- Table 9: When to use each design --------------------------------------

t9 <- data.frame(
  Situation = c("Regulatory Phase III, equal allocation expected",
                "Strong known prognostic factors (age, severity)",
                "Severe disease, large anticipated effect",
                "Pediatric or oncology trial",
                "Small pilot / Phase I",
                "When transparency is paramount",
                "Adaptive trial with ethical urgency"),
  Recommendation = c("Pocock-Simon",
                     "Pocock-Simon",
                     "RPW",
                     "RPW",
                     "Simple",
                     "Simple",
                     "RPW"),
  Rationale = c("Maximum precision; balances on prognostic factors",
                "Minimization removes random imbalance on key covariates",
                "Allocates more patients to better arm",
                "Ethical priority over statistical efficiency",
                "No assumptions, robust, easy to explain",
                "Pure randomization is least open to manipulation",
                "Real-time response-adaptive allocation favors winning arm")
)
show(t9, "Table 9: When to use each design")
write.csv(t9, file.path(out_dir, "table_9_when_to_use.csv"), row.names = FALSE)

# ---- Done ------------------------------------------------------------------

cat("=== All tables built ===\n")
cat("Files in", out_dir, ":\n")
print(list.files(out_dir, pattern = "^table_"))