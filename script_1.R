################################################################################
# Screen-detected high blood pressure: NHANES 2017-2018
#
# Research question:
# Do HbA1c, HDL cholesterol and hs-CRP improve identification of screen-detected 
# high BP beyond age, sex and BMI among adults without previously reported hypertension?
#
# Study design:
# Cross-sectional analysis of NHANES 2017-2018.
#
# Data source:
# Public-use NHANES examination, questionnaire and lab data.
################################################################################

# 1. Clear Environment and load any required packages ##########################

rm(list = ls())
if(!require('pacman'))install.packages('pacman')
pacman::p_load(haven, dplyr, ggplot2, survey, pROC, ranger, rsample, tidyr)

# 2. Download and import NHANES datasets required for analysis #################

files <- c("BPX_J", "BPQ_J", "BMX_J", "DEMO_J", "GHB_J", "HDL_J", "HSCRP_J")
base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/"
dir.create("data_raw", showWarnings = FALSE)

  # Only download file if it is not already locally present

for (file in files) {
  destination <- file.path("data_raw", paste0(file, ".XPT"))
  if (!file.exists(destination)) {
    download.file(paste0(base_url, file, ".XPT"), destfile = destination, 
      mode = "wb")}}

bpx <- read_xpt("data_raw/BPX_J.XPT") # BP examination
bpq <- read_xpt("data_raw/BPQ_J.XPT") # BP questionnaire
bmx <- read_xpt("data_raw/BMX_J.XPT") # Body measurements - BMI
demo <- read_xpt("data_raw/DEMO_J.XPT") # Demographics and other required info
ghb <- read_xpt("data_raw/GHB_J.XPT") # HbA1c
hdl <- read_xpt("data_raw/HDL_J.XPT") # HDL cholesterol
hscrp <- read_xpt("data_raw/HSCRP_J.XPT") # High sensitivity CRP

# 3. Define initial population for analysis ####################################

demo_processed <- demo %>%
  select(SEQN, RIDSTATR, RIAGENDR, RIDAGEYR,
    RIDRETH3, WTMEC2YR, SDMVPSU, SDMVSTRA)

demo_processed %>%
  summarise(n = n(), 
    missing_gender = sum(is.na(RIAGENDR)), 
    missing_age = sum(is.na(RIDAGEYR)), 
    missing_ethnicity = sum(is.na(RIDRETH3)), 
    missing_weight = sum(is.na(WTMEC2YR)))

  # Filter to get only MEC examined adults older than 20

demo_processed <- demo_processed %>%
  filter(RIDAGEYR>=20, RIDSTATR == 2, WTMEC2YR > 0)

# 4. Validate source datasets ##################################################

datasets <- list(bmx = bmx, 
                 bpq = bpq, 
                 bpx = bpx, 
                 demo_processed = demo_processed, 
                 ghb = ghb, 
                 hdl = hdl, 
                 hscrp = hscrp)

  # Ensure that all SEQN numbers are unique

sapply(datasets, function(x) {
  c(rows = nrow(x), 
    unique_SEQN = n_distinct(x$SEQN),
    duplicate_SEQN = sum(duplicated(x$SEQN)))
})

  # Ensure that variables required for analysis are all present before merging

variables_needed <- list(
  bmx = c("SEQN", "BMXBMI"), 
  bpq = c("SEQN", "BPQ020"),
  bpx = c("SEQN", "BPXSY1", "BPXDI1", "BPXSY2", "BPXDI2", "BPXSY3", "BPXDI3", 
          "BPXSY4", "BPXDI4"),
  demo_processed = c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "WTMEC2YR", "SDMVPSU", 
           "SDMVSTRA"),
  ghb = c("SEQN", "LBXGH"), 
  hdl = c("SEQN", "LBDHDD"), 
  hscrp = c("SEQN", "LBXHSCRP"))

for (name in names(variables_needed)) {
  missing_columns <- setdiff(variables_needed[[name]], 
                             names(datasets[[name]]))
  if (length(missing_columns) > 0) {
    stop("Missing columns in ", 
         name, ": ", 
         paste(missing_columns, collapse = ", ")
    )}}

# 6. Merge datasets together using SEQN number #################################
# Left join preserves missing values as NA 

analysis_raw <- demo_processed %>%
  left_join(bmx %>% select(SEQN, BMXBMI), by = "SEQN") %>%
  left_join(bpq %>% select(SEQN, BPQ020), by = "SEQN") %>%
  left_join(bpx %>% select(SEQN, BPXSY1, BPXDI1, BPXSY2, BPXDI2, BPXSY3, BPXDI3, 
    BPXSY4, BPXDI4), by = "SEQN") %>%
  left_join(ghb %>% select(SEQN, LBXGH), by = "SEQN") %>%
  left_join(hdl %>% select(SEQN, LBDHDD), by = "SEQN") %>%
  left_join(hscrp %>% select(SEQN, LBXHSCRP), by = "SEQN")

nrow(analysis_raw)

  # Check missing NA values

analysis_raw %>%
  summarise(n = n(), 
    missing_bmi = sum(is.na(BMXBMI)),
    missing_bp_question = sum(is.na(BPQ020)), 
    missing_gender = sum(is.na(RIAGENDR)), 
    missing_age = sum(is.na(RIDAGEYR)),
    missing_ethnicity = sum(is.na(RIDRETH3)),
    missing_weight = sum(is.na(WTMEC2YR)), 
    missing_hba1c = sum(is.na(LBXGH)),
    missing_hdl = sum(is.na(LBDHDD)), 
    missing_crp = sum(is.na(LBXHSCRP)))

# 7. Asess completeness of BP measurements #####################################

analysis_raw <- analysis_raw %>%
  mutate(
    valid_pair1 = !is.na(BPXSY1) & !is.na(BPXDI1),
    valid_pair2 = !is.na(BPXSY2) & !is.na(BPXDI2),
    valid_pair3 = !is.na(BPXSY3) & !is.na(BPXDI3),
    valid_pair4 = !is.na(BPXSY4) & !is.na(BPXDI4),
    num_bp_pairs = rowSums(across(c(
      valid_pair1, valid_pair2, valid_pair3, valid_pair4))),
    valid_bp = num_bp_pairs == 3)

table(analysis_raw$num_bp_pairs)



  # Calculate mean blood pressure using only complete paired readings

analysis_raw <- analysis_raw %>%
  mutate(mean_sbp = if_else(valid_bp,(
      coalesce(if_else(valid_pair1, BPXSY1, NA_real_), 0) +
      coalesce(if_else(valid_pair2, BPXSY2, NA_real_), 0) +
      coalesce(if_else(valid_pair3, BPXSY3, NA_real_), 0) +
      coalesce(if_else(valid_pair4, BPXSY4, NA_real_), 0)
      ) / 3, NA_real_),
    mean_dbp = if_else(valid_bp,(
      coalesce(if_else(valid_pair1, BPXDI1, NA_real_), 0) +
      coalesce(if_else(valid_pair2, BPXDI2, NA_real_), 0) +
      coalesce(if_else(valid_pair3, BPXDI3, NA_real_), 0) +
      coalesce(if_else(valid_pair4, BPXDI4, NA_real_), 0)
      ) / 3, NA_real_))

  # Define eligible population (1. no previous hypertension diagnosis, 
  # 2. three valid and paired BP measurements)

target_population <- analysis_raw %>%
  filter(BPQ020 == 2, valid_bp)
target_population %>%
  filter(mean_dbp == 0) %>%
  select(SEQN, RIDAGEYR, RIAGENDR, BPXSY1, BPXDI1, BPXSY2, BPXDI2, BPXSY3, 
         BPXDI3, BPXSY4, BPXDI4, mean_sbp, mean_dbp)

  # Identify missing biomarker or bmi participants but do not remove them

target_population <- target_population %>%
  mutate(complete_predictors = 
    !is.na(BMXBMI) & !is.na(LBXGH) &
    !is.na(LBDHDD) & !is.na(LBXHSCRP))

table(target_population$complete_predictors)

# 8. Define binary screening outcome using mean BP #############################
# High BP is defined here as SBP >= 140 mmHg or DBP >= 90 mmHG
# It is not diagnosed hypertension, rather it is screen-detected elevated BP

hist(target_population$mean_sbp, breaks = 30,
  main = "Distribution of Mean Systolic Blood Pressure",
  xlab = "Mean systolic blood pressure (mmHg)")

hist(target_population$mean_dbp, breaks = 30,
  main = "Distribution of Mean Diastolic Blood Pressure",
  xlab = "Mean diastolic blood pressure (mmHg)")

target_population <- target_population %>%
  mutate(high_bp = if_else(mean_sbp >= 140 | mean_dbp >= 90, 1L, 0L))

summary(target_population$mean_sbp)
summary(target_population$mean_dbp)
table(target_population$high_bp)


# 9. Add labels to binary outcomes #############################################

target_population <- target_population %>%
  mutate(high_bp_label = factor(high_bp, levels = c(0, 1), 
    labels = c("No high BP", "Screen-detected high BP")), 
    sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female"))
  )

# 10. Inspect biomarker distribution ###########################################

hist(target_population$LBXHSCRP, breaks = 40, 
     main = "Distribution of hs-CRP", xlab = "hs-CRP (mg/L)")
hist(target_population$LBXGH, breaks = 40, 
     main = "Distribution of HbA1c", xlab = "HbA1c %")
hist(target_population$LBDHDD, breaks = 40, 
     main = "Distribution of HDL", xlab = "HDL (mg/dL)")

# 11. Unweighted statistical analysis ##########################################
# Only describes analysed participants rather than entire population
  # Continuous variable stats

target_population %>%
  group_by(high_bp_label) %>%
  summarise(n = n(),
    age_mean = mean(RIDAGEYR, na.rm = TRUE), age_sd = sd(RIDAGEYR, na.rm = TRUE),
    bmi_mean = mean(BMXBMI, na.rm = TRUE), bmi_sd = sd(BMXBMI, na.rm = TRUE),
    hba1c_mean = mean(LBXGH, na.rm = TRUE), hba1c_sd = sd(LBXGH, na.rm = TRUE),
    hdl_mean = mean(LBDHDD, na.rm = TRUE), hdl_sd = sd(LBDHDD, na.rm = TRUE),
    
    # Median and IQR used due to significant skew
    
    crp_median = median(LBXHSCRP, na.rm = TRUE), crp_iqr = IQR(LBXHSCRP, na.rm = TRUE),
    sbp_mean = mean(mean_sbp), dbp_mean = mean(mean_dbp)
  )

  # Discrete variable stats

target_population %>%
  count(high_bp_label, sex) %>%
  group_by(high_bp_label) %>%
  mutate(percent = 100 * n / sum(n))

target_population %>%
  group_by(sex) %>%
  summarise(n = n(), high_bp_n = sum(high_bp), 
    high_bp_percent = 100 * mean(high_bp))

  # Age groups introduced to make figures easier to interpret

target_population <- target_population %>%
  mutate(age_group = cut(RIDAGEYR, 
      breaks = c(20, 30, 40, 50, 60, 70, Inf), right = FALSE,
      labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70+")))

age_results <- target_population %>%
  group_by(age_group) %>%
  summarise(n = n(), high_bp_n = sum(high_bp), 
  high_bp_percent = 100 * mean(high_bp))
age_results

  # Unweighted high BP frequency by age group

ggplot(age_results, aes(x = age_group, y = high_bp_percent)) +
  geom_col() +
  labs(title = "Screen-detected high blood pressure by age: unweighted sample",
    subtitle = "Adults without previously reported hypertension, NHANES 2017-2018",
    x = "Age group",
    y = "Participants with high blood pressure (%)") +
  theme_minimal()

# 12. Generate NHANES survey-design object #####################################
# Estimate prevalence to get survey weighted stats as NHANES is not a random
# sample - it uses unequal selection probabilities. Calculate based on full popn.
# Sampling weight gives an idea of how much an observed participant contributes
# to the target popn. estimate

survey_source <- analysis_raw %>%
  mutate(eligible_screening = !is.na(BPQ020) & BPQ020 == 2 & valid_bp,
    high_bp = if_else(valid_bp, as.integer(mean_sbp >= 140 | mean_dbp >= 90), NA_integer_),
    sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
    age_group = cut(RIDAGEYR, breaks = c(20, 30, 40, 50, 60, 70, Inf),
      right = FALSE, labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70+")))

full_design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
  nest = TRUE, data = survey_source)

nhanes_design <- subset(full_design, eligible_screening)

  # Estimate weighted proportion with elevated BP

weighted_prevalence <- svyciprop(~I(high_bp == 1), 
                        design = nhanes_design, method = "logit")
weighted_prevalence
confint(weighted_prevalence)

  # Include weighted prevalence for age group and sex

weighted_age <- svyby(~high_bp, ~age_group, nhanes_design, svymean, 
                      na.rm = TRUE, vartype = "ci")
weighted_age

weighted_sex <- svyby(~high_bp, ~sex, nhanes_design, svymean, 
                      na.rm = TRUE, vartype = "ci")
weighted_sex

# 13. Handle CRP skew ##########################################################
target_population %>%
  summarise(
    n = sum(!is.na(LBXHSCRP)),
    min = min(LBXHSCRP, na.rm = TRUE),
    q25 = quantile(LBXHSCRP, 0.25, na.rm = TRUE),
    median = median(LBXHSCRP, na.rm = TRUE),
    q75 = quantile(LBXHSCRP, 0.75, na.rm = TRUE),
    q95 = quantile(LBXHSCRP, 0.95, na.rm = TRUE),
    q99 = quantile(LBXHSCRP, 0.99, na.rm = TRUE),
    max = max(LBXHSCRP, na.rm = TRUE)
  )

  # Use log transformation to supress right skew 
  # log2 corresponds to one unit doubling hs-CRP

target_population <- target_population %>%
  mutate(log2_crp = log2(LBXHSCRP))
hist(target_population$log2_crp, breaks = 40, 
     main = "Distribution of log transformed hs-CRP", xlab = "log2 hs-CRP (mg/L)")

# 14. Generate model dataset ###################################################
# Create predictor variables on the survey source data.

survey_source <- survey_source %>%
  mutate(decade_age = RIDAGEYR / 10, log2_crp = if_else(
      !is.na(LBXHSCRP) & LBXHSCRP > 0, log2(LBXHSCRP), NA_real_),
    
    eligible_model = eligible_screening &
      !is.na(BMXBMI) & !is.na(LBXGH) & !is.na(LBDHDD) & !is.na(log2_crp))

model_data <- target_population %>%
  filter(complete_predictors) %>%
  mutate(decade_age = RIDAGEYR / 10)

  # Rebuild the full design after adding variables.

full_design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
  nest = TRUE, data = survey_source)

  # Generate model

model_design <- subset(full_design, eligible_model)

  # Model 1 - How well can basic characteristics identify screen-detected high 
  # BP using using a survey-weighted logistic-type model?

model_1 <- svyglm(high_bp ~ decade_age + sex + BMXBMI, design = model_design, 
                  family = quasibinomial())
summary(model_1)
odds_ratio_model_1 <- cbind(OR = exp(coef(model_1)), exp(confint(model_1)))
odds_ratio_model_1

  # Model 2 - Same baseline as model 1 but adding in biomarkers

model_2 <- svyglm(high_bp ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp, 
                  design = model_design, family = quasibinomial())
summary(model_2)
odds_ratio_model_2 <- cbind(OR = exp(coef(model_2)),exp(confint(model_2)))
odds_ratio_model_2
exp(10 * coef(model_2)["LBDHDD"])
exp(10 * confint(model_2)["LBDHDD", ])

weighted_prevalence
confint(weighted_prevalence)

summary(model_1)
summary(model_2)

# 15. Investigate unexpected high HDL confidence ###############################
# Model 2 identified a positive association between HDL and high BP
target_population %>%
  group_by(sex, high_bp_label) %>%
  summarise(n = sum(!is.na(LBDHDD)), mean_hdl = mean(LBDHDD, na.rm = TRUE),
    sd_hdl = sd(LBDHDD, na.rm = TRUE), median_hdl = median(LBDHDD, na.rm = TRUE),
    .groups = "drop")

target_population %>%
  group_by(age_group) %>%
  summarise(n = sum(!is.na(LBDHDD)), mean_hdl = mean(LBDHDD, na.rm = TRUE),
    sd_hdl = sd(LBDHDD, na.rm = TRUE), median_hdl = median(LBDHDD, na.rm = TRUE),
    .groups = "drop")

  # Exploratory analysis of potential non-linearity of HDL

model_2_hdl_nonlinear <- svyglm(high_bp ~ decade_age + sex + BMXBMI + LBXGH + 
  LBDHDD + I(LBDHDD^2) + log2_crp,
  design = model_design, family = quasibinomial())
summary(model_2_hdl_nonlinear)

  # Design-adjusted test of HDL^2 contribution - no strong evidence found

regTermTest(model_2_hdl_nonlinear,~ I(LBDHDD^2))

  # Investigate whether biomarkers collectively improve model after adjustment
  # for age, sex and BMI

regTermTest(model_2,~ LBXGH + LBDHDD + log2_crp)

# Generate prediction probabilities
model_data <- model_data %>%
  mutate(pred_model_1 = predict(model_1, type = "response"),
    pred_model_2 = predict(model_2, type = "response"))

summary(model_data$pred_model_1)
summary(model_data$pred_model_2)

# 16. AUROC analysis ###########################################################
# Compare how well each model can distinguish someone with screen-detected 
# high BP versus someone without it using. ROC curve gives trade-off 
# between sensitivity and specificity. DeLong test compares the AUROCs of the 
# two models but does not acount for NHANES design and is evaluated on 
# training data. 0.5 indicates chance ranking, 1 indicates perfect ranking

roc_model_1 <- roc(response = model_data$high_bp,predictor = model_data$pred_model_1)
roc_model_2 <- roc(response = model_data$high_bp,predictor = model_data$pred_model_2)
roc_test <- roc.test(roc_model_1, roc_model_2, paired = TRUE,method = "delong")

auc(roc_model_1)
auc(roc_model_2)
roc_test

# 17. Model callibration using Brier score #####################################
# Refelcts probability accuracy. Lower scores indicate better
# overall probability accuracy

brier_model_1 <- mean((model_data$pred_model_1 - model_data$high_bp)^2)
brier_model_2 <- mean((model_data$pred_model_2 - model_data$high_bp)^2)

brier_model_1
brier_model_2

# 18. Model calibration by predicted risk ######################################
# Participants into 10 equal sized groups by predicted probabilities. Compare
# mean predicted risk with observed frequency. Perfect calibration would lie on
# y=x. Summaries are unweighted

calibration_1 <- model_data %>%
  mutate(risk_decile = ntile(pred_model_1, 10)) %>%
  group_by(risk_decile) %>%
  summarise(n = n(), mean_predicted = mean(pred_model_1), 
            observed = mean(high_bp), .groups = "drop")
calibration_1

calibration_2 <- model_data %>%
  mutate(risk_decile = ntile(pred_model_2, 10)) %>%
  group_by(risk_decile) %>%
  summarise(n = n(), mean_predicted = mean(pred_model_2), 
            observed = mean(high_bp), .groups = "drop")
calibration_2

  # Combine calibration summaries for comparison.

calibration_combined <- bind_rows(
  calibration_1 %>% mutate(model = "Baseline logistic"),
  calibration_2 %>% mutate(model = "Biomarker logistic"))
ggplot(calibration_combined, aes(
    x = mean_predicted, y = observed,
    colour = model, group = model)) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.5) +
  geom_line() +
  coord_equal() +
  labs(title = "Apparent calibration of logistic models",
    subtitle = "NHANES 2017–2018; unweighted risk deciles",
    x = "Mean predicted probability",
    y = "Observed proportion with high BP",
    colour = "Model") +
  theme_minimal()

# 19. Fit original logistic models on the full dataset #########################
# Fit the original logistic models on the full modelling dataset.

original_model_1 <- glm(high_bp ~ decade_age + sex + BMXBMI,
  data = model_data, family = binomial())

original_model_2 <- glm(high_bp ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
  data = model_data, family = binomial())

  # Calculate apparent AUROC for each original model.
apparent_auc_1 <- as.numeric(auc(roc(model_data$high_bp,
  predict(original_model_1, newdata = model_data, type = "response"),
  levels = c(0, 1), direction = "<", quiet = TRUE)))

apparent_auc_2 <- as.numeric(auc(roc(model_data$high_bp,
  predict(original_model_2, newdata = model_data, type = "response"),
  levels = c(0, 1), direction = "<", quiet = TRUE)))

# 20. Bootstrap resampling #####################################################
# Help to estimate model optimism. GLM used as part of sample-level not 
# population-level prediction analysis
# 500 bootstraps for model 1

set.seed(2004)
BS <- 500
optimism_model_1 <- numeric(BS)

  # Draw BS sample

for (bs in 1:BS) {
  bootstrap_values <- sample(seq_len(nrow(model_data)), size = nrow(model_data),
    replace = TRUE)
  bootstrap_data <- model_data[bootstrap_values, ]
  
  # Fit model 1 in BS sample

  bootstrap_model <- glm(high_bp ~ decade_age + sex + BMXBMI,
    data = bootstrap_data, family = binomial())
  
  # Predict on bootstrap sample
  
  pred_boot <- predict(bootstrap_model, newdata = bootstrap_data, type = "response")
  
  # Predict on original sample
  
  pred_orig <- predict(bootstrap_model, newdata = model_data, type = "response")
  
  # Calculate AUROC
  
  auc_boot <- as.numeric(auc(roc(bootstrap_data$high_bp, pred_boot, quiet = TRUE)))
  auc_orig <- as.numeric(auc(roc(model_data$high_bp, pred_orig, quiet = TRUE)))
  
  optimism_model_1[bs] <- auc_boot - auc_orig
}

  # 500 bootstraps for model 2

set.seed(2004)
BS <- 500
optimism_model_2 <- numeric(BS)

  # Draw BS sample

for (bs in 1:BS) {
  bootstrap_values <- sample(seq_len(nrow(model_data)), size = nrow(model_data),
                             replace = TRUE)
  bootstrap_data <- model_data[bootstrap_values, ]
  
  # Fit model 2 in BS sample
  
  bootstrap_model <- glm(high_bp ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
                         data = bootstrap_data, family = binomial())
  
  # Predict on bootstrap sample
  
  pred_boot <- predict(bootstrap_model, newdata = bootstrap_data, type = "response")
  
  # Predict on original sample
  
  pred_orig <- predict(bootstrap_model, newdata = model_data, type = "response")
  
  # Calculate AUROC
  
  auc_boot <- as.numeric(auc(roc(bootstrap_data$high_bp, pred_boot, quiet = TRUE)))
  auc_orig <- as.numeric(auc(roc(model_data$high_bp, pred_orig, quiet = TRUE)))
  
  optimism_model_2[bs] <- auc_boot - auc_orig
}

  # Average optimism across 500 bootstrap samples

mean_optimism_1 <- mean(optimism_model_1)
mean_optimism_2 <- mean(optimism_model_2)

  # Correct original apparent AUROC for estimated optimism

corrected_auc_1 <- apparent_auc_1 - mean_optimism_1
corrected_auc_2 <- apparent_auc_2 - mean_optimism_2

bootstrap_summary <- data.frame(
  model = c("Baseline logistic", "Biomarker logistic"),
  apparent_auc = c(apparent_auc_1, apparent_auc_2),
  mean_optimism = c(mean_optimism_1, mean_optimism_2),
  corrected_auc = c(corrected_auc_1, corrected_auc_2)
)

bootstrap_summary


  # Convert high BP into a factor for ranger package
model_data <- model_data %>%
  mutate(high_bp_rf = factor(high_bp, levels = c(0, 1), labels = c("No", "Yes")))

# Fitting initial random forest without rebalancing
set.seed(2004)

# 21. Set up RF models #########################################################
# Want predicted probabilities (probability) out of it alongside how model
# prediction error changes when a predictors values are shuffled (importance).

rf_model <- ranger(high_bp_rf ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
  data = model_data, probability = TRUE, num.trees = 1000, importance = "permutation",
  seed = 2004)
rf_model

  # OOB represents the original observations not selected for a particular tree
  # Gives an internal performance diagnostic

rf_oob_prob <- rf_model$predictions[,"Yes"]
rf_oob_roc <- roc(response = model_data$high_bp, predictor = rf_oob_prob, quiet = TRUE)
auc(rf_oob_roc)

  # Interpret variable importance 

rf_importance <- sort(rf_model$variable.importance, decreasing = TRUE)
rf_importance

  # Vary mtry (how many predictors are considered at each tree split) to 
  # minimise Brier and maximise AUROC scores

mtry_results <- data.frame(mtry = 1:6, oob_auc = NA_real_, oob_brier = NA_real_)

for (m in 1:6) {
  rf_temp <- ranger(high_bp_rf ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
    data = model_data, probability = TRUE, num.trees = 1000, mtry = m, seed = 2004)
  oob_prob <- rf_temp$predictions[, "Yes"]
  mtry_results$oob_auc[m] <- as.numeric(
    auc(roc(model_data$high_bp, oob_prob, quiet = TRUE)))
  mtry_results$oob_brier[m] <- mean((oob_prob - model_data$high_bp)^2)
}
mtry_results

  # Vary minimum node size (smaller values allow for more detailed splits) 
  # whilst keeping mtry at 1

node_sizes <- c(5, 10, 20, 30, 50, 75, 100, 500)
node_results <- data.frame(min_node_size = node_sizes, oob_auc = NA_real_,
  oob_brier = NA_real_)

for (i in seq_along(node_sizes)) {
  rf_temp <- ranger(high_bp_rf ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
    data = model_data, probability = TRUE, num.trees = 1000, mtry = 1,
    min.node.size = node_sizes[i], seed = 2004)
  oob_prob <- rf_temp$predictions[, "Yes"]
  node_results$oob_auc[i] <- as.numeric(auc(roc(model_data$high_bp,
        oob_prob, quiet = TRUE)))

  node_results$oob_brier[i] <- mean((oob_prob - model_data$high_bp)^2)
}
node_results

  # Dial in node sizes (Optimum is 50)

node_sizes_specific <- c(50, 60, 70, 80, 90, 100, 120, 140)
node_results_specific <- data.frame(min_node_size = node_sizes_specific, oob_auc = NA_real_,
  oob_brier = NA_real_)

for (i in seq_along(node_sizes_specific)) {
  rf_temp <- ranger(high_bp_rf ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
    data = model_data, probability = TRUE, num.trees = 1000, mtry = 1,
    min.node.size = node_sizes_specific[i], seed = 2004)
  oob_prob <- rf_temp$predictions[, "Yes"]
  node_results_specific$oob_auc[i] <- as.numeric(auc(roc(model_data$high_bp,
    oob_prob, quiet = TRUE)))
  node_results_specific$oob_brier[i] <- mean((oob_prob - model_data$high_bp)^2)
}
node_results_specific

  # Use set mtry and node sizes to generate an exploratory RF estimate via
  # OOB tuning. Not a final validated model

set.seed(2004)

rf_exploratory <- ranger(high_bp_rf ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
  data = model_data, probability = TRUE, num.trees = 1000, mtry = 1,
  min.node.size = 50, importance = "permutation", seed = 2004)

rf_exploratory

rf_exploratory_prob <- rf_exploratory$predictions[, "Yes"]
rf_exploratory_roc <- roc(response = model_data$high_bp, predictor = rf_exploratory_prob,
  quiet = TRUE)

auc(rf_exploratory_roc)

# 22. Nested cross-validation ##################################################
# Compare AUROC and Brier score on the same participants
# Model 1 = Logistic regression: age + sex + BMI
# Model 2 = Logistic regression: Model 1 + biomarkers
# Model 3 = Tuned random forest using Model 2 predictors
# Outer fold estimates performance on untested participants
# Inner fold tunes RF parameters without touching outer test data

set.seed(2004)
outer_folds <- vfold_cv(model_data, v = 5, strata = high_bp_rf)

  # Define random forest hyperparameter grid based on exploratory OOB analysis

rf_grid <- expand.grid(mtry = c(1, 2, 3), min_node_size = c(75, 100, 125, 150, 175, 200))

  # Create table to store outer-fold results

nested_results <- data.frame(fold = 1:5,
  logistic_1_auc = NA_real_, logistic_2_auc = NA_real_, rf_auc = NA_real_,
  logistic_1_brier = NA_real_, logistic_2_brier = NA_real_, rf_brier = NA_real_,
  best_mtry = NA_integer_, best_node_size = NA_integer_)

  # Store individual predictions for calibration analysis where each participant
  # recives one prediction from each model when they are in the outer test fold

oof_predictions <- data.frame(SEQN = model_data$SEQN, high_bp = model_data$high_bp,
  fold = NA_integer_, logistic_1_prob = NA_real_, logistic_2_prob = NA_real_,
  rf_prob = NA_real_)

# Begin outer cross-validation

for (fold in 1:5) {
  cat("\nProcessing outer fold", fold, "of 5...\n")
  
  # Split into outer training and untouched test data
  
  outer_train <- analysis(outer_folds$splits[[fold]])
  outer_test <- assessment(outer_folds$splits[[fold]])

  # Logistic Model 1
  
  logistic_1 <- glm(high_bp ~ decade_age + sex + BMXBMI, data = outer_train,
    family = binomial())
  pred_logistic_1 <- predict(logistic_1, newdata = outer_test, type = "response")
  nested_results$logistic_1_auc[fold] <- as.numeric(
    auc(roc(response = outer_test$high_bp, predictor = pred_logistic_1, 
      levels = c(0, 1), direction = "<", quiet = TRUE)))
  nested_results$logistic_1_brier[fold] <- mean((pred_logistic_1 - outer_test$high_bp)^2)

  # Logistic Model 2
  
  logistic_2 <- glm(high_bp ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
    data = outer_train, family = binomial())
  pred_logistic_2 <- predict(logistic_2, newdata = outer_test, type = "response")
  nested_results$logistic_2_auc[fold] <- as.numeric(
    auc(roc(response = outer_test$high_bp, predictor = pred_logistic_2, 
      levels = c(0, 1), direction = "<", quiet = TRUE)))
  nested_results$logistic_2_brier[fold] <- mean((pred_logistic_2 - outer_test$high_bp)^2)

  # Create inner CV folds using only outer training data
  
  set.seed(2004 + fold)
  inner_folds <- vfold_cv(outer_train, v = 5, strata = high_bp_rf)

  # Create table for RF tuning results
  
  tuning_results <- rf_grid %>%
    mutate(mean_auc = NA_real_, sd_auc = NA_real_)

  # Test every RF hyperparameter combination
  
  for (g in 1:nrow(rf_grid)) {
    inner_auc <- numeric(5)
    
    # Run 5-fold inner CV for this hyperparameter combination
    
    for (i in 1:5) {
      inner_train <- analysis(inner_folds$splits[[i]])
      inner_validation <- assessment(inner_folds$splits[[i]])

      # Fit random forest to inner training data
      
      rf_inner <- ranger(high_bp_rf ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
        data = inner_train, probability = TRUE, num.trees = 1000, mtry = rf_grid$mtry[g],
        min.node.size = rf_grid$min_node_size[g], seed = 2004 + fold * 100 + g * 10 + i)
      
      # Predict held-out inner validation data
      
      validation_prob <- predict(rf_inner,data = inner_validation)$predictions[, "Yes"]
    
      # Calculate validation AUROC
      
      inner_auc[i] <- as.numeric(auc(roc(response = inner_validation$high_bp,
        predictor = validation_prob, quiet = TRUE)))
    }
    
    # Average performance across the five inner folds
    
    tuning_results$mean_auc[g] <- mean(inner_auc)
    tuning_results$sd_auc[g] <- sd(inner_auc)
  }
  
  # Select best RF hyperparameters based on highest mean AUROC, lower SD used 
  # in the case of ties
  
  best_parameters <- tuning_results %>%
    arrange(desc(mean_auc), sd_auc) %>%
    slice(1)
  
  best_mtry <- best_parameters$mtry
  best_node_size <- best_parameters$min_node_size
  
  # Store selected parameters
  
  nested_results$best_mtry[fold] <- best_mtry
  nested_results$best_node_size[fold] <- best_node_size
  
  # Fit tuned RF to ALL outer training data
  
  rf_final_fold <- ranger(high_bp_rf ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
    data = outer_train, probability = TRUE, num.trees = 1000,mtry = best_mtry,
    min.node.size = best_node_size, seed = 2004 + fold)

  # Predict the completely untouched outer test fold
  
  rf_test_prob <- predict(rf_final_fold,
    data = outer_test)$predictions[, "Yes"]
  
  nested_results$rf_auc[fold] <- as.numeric(auc(roc(
      response = outer_test$high_bp, predictor = rf_test_prob,
      levels = c(0, 1), direction = "<", quiet = TRUE)))
  
  nested_results$rf_brier[fold] <- mean(
    (rf_test_prob - outer_test$high_bp)^2)
  
  # Save held-out participant predictions. Match by SEQN rather than assuming 
  # that order of outer_test matches order of rows in the original model_data.
  
  test_rows <- match(outer_test$SEQN, oof_predictions$SEQN)
  oof_predictions$fold[test_rows] <- fold
  oof_predictions$logistic_1_prob[test_rows] <- pred_logistic_1
  oof_predictions$logistic_2_prob[test_rows] <- pred_logistic_2
  oof_predictions$rf_prob[test_rows] <- rf_test_prob
}
  
  # Check output completeness such that every participant should have one 
  # prediction from each of the three models

stopifnot(!anyDuplicated(oof_predictions$SEQN))
stopifnot(all(complete.cases(oof_predictions[, c(
    "fold", "logistic_1_prob", "logistic_2_prob", "rf_prob")])))

# 23. Performance summary ######################################################

nested_results

performance_summary <- nested_results %>%
  summarise(
    logistic_1_mean_auc = mean(logistic_1_auc), logistic_1_sd_auc = sd(logistic_1_auc),
    logistic_2_mean_auc = mean(logistic_2_auc), logistic_2_sd_auc = sd(logistic_2_auc),
    rf_mean_auc = mean(rf_auc), rf_sd_auc = sd(rf_auc),
    logistic_1_mean_brier = mean(logistic_1_brier), logistic_1_sd_brier = sd(logistic_1_brier),
    logistic_2_mean_brier = mean(logistic_2_brier), logistic_2_sd_brier = sd(logistic_2_brier),
    rf_mean_brier = mean(rf_brier), rf_sd_brier = sd(rf_brier))
performance_summary

# 24. Research figures #########################################################
# Figure 1

dir.create("figures", showWarnings = FALSE)
flow_counts <- data.frame(
  stage = c("NHANES 2017-2018 participants",
    "Adults aged 20+, MEC examined,\npositive examination weights",
    "No previously reported hypertension",
    "Three valid paired BP measurements",
    "Complete predictor data"),
  n = c(nrow(demo), nrow(demo_processed),
      sum(analysis_raw$BPQ020 == 2, na.rm = TRUE), nrow(target_population),
      nrow(model_data))) %>%
  mutate(y = rev(seq_len(n())),
    label = paste0(stage, "\n", format(n, big.mark = ",")))

flow_plot <- ggplot(flow_counts, aes(x = 1, y = y)) +
  geom_segment(aes(x = 1, xend = 1, y = y - 1),
    data = subset(flow_counts, y > 1), colour = "black",
    linewidth = 0.7) +
  geom_label(aes(label = label), fill = "white", linewidth = 0.3,
    size = 3.8, lineheight = 1.2) +
  coord_cartesian(xlim = c(0.5, 1.5), ylim = c(0.5, 5.5), clip = "off") +
  theme_void() +
  labs(title = "Study participant selection",
    subtitle = "NHANES 2017–2018"
  ) +
  theme(plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "black"))
flow_plot
ggsave(filename = "figures/figure1_participant_flow.png",
  plot = flow_plot, width = 8, height = 9, dpi = 300)

# Figure 2

age_plot_data <- weighted_age %>%
  as.data.frame() %>%
  mutate(prevalence = 100 * high_bp, lower_ci = 100 * ci_l, upper_ci = 100 * ci_u)

figure2 <- ggplot(age_plot_data, aes(x = age_group, y = prevalence)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.15, linewidth = 0.6, colour = "black") +
  geom_point(size = 3, colour = "black") +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 5),
    labels = function(x) paste0(x, "%")) +
  labs(title = "Screen-detected high blood pressure increases with age",
    subtitle = "NHANES 2017–2018 | Survey-weighted estimates",
    x = "Age group (years)",
    y = "Prevalence of screen-detected high BP",
    caption = "Points represent weighted prevalence; bars indicate 95% confidence intervals.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
figure2

ggsave("figures/figure2_weighted_prevalence_age.png",
  plot = figure2, width = 8, height = 5, dpi = 300)

# Figure 3

beta <- coef(model_2)
ci <- confint(model_2)
terms <- c("decade_age", "sexFemale", "BMXBMI",
  "LBXGH", "LBDHDD", "log2_crp")
multipliers <- c(1, 1, 1, 1, 10, 1)
or_plot_data <- data.frame(
  predictor = c("Age (per 10 years)", "Female (vs male)", "BMI (per kg/m²)",
    "HbA1c (per 1 percentage point)", "HDL (per 10 mg/dL)", "hs-CRP (per doubling)"),
  OR = exp(beta[terms] * multipliers),
  lower = exp(ci[terms, 1] * multipliers),
  upper = exp(ci[terms, 2] * multipliers))

or_plot_data$predictor <- factor(or_plot_data$predictor,
  levels = rev(or_plot_data$predictor))

figure3 <- ggplot(or_plot_data, aes(x = OR, y = predictor)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "blue") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
    height = 0.18, linewidth = 0.65, colour = "darkred") +
  geom_point(size = 3, colour = "darkred") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 1.5, 2, 3)) +
  labs(title = "Adjusted associations with screen-detected high BP",
    subtitle = "Survey-weighted multivariable logistic regression",
    x = "Adjusted odds ratio (log scale; 95% CI)",
    y = NULL,
    caption = paste(
      "NHANES 2017–2018.",
      "All predictors adjusted for one another.",
      "Associations do not establish causation.")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
figure3
ggsave("figures/figure3_adjusted_odds_ratios.png",
  plot = figure3, width = 9, height = 5.5, dpi = 300)

# Figure 4

performance_plot_data <- nested_results %>%
  select(fold, logistic_1_auc, logistic_2_auc, rf_auc,
    logistic_1_brier, logistic_2_brier, rf_brier) %>%
  pivot_longer(cols = -fold, names_to = "model_metric", values_to = "value") %>%
  mutate(
    model = case_when(
      grepl("^logistic_1", model_metric) ~ "Baseline logistic",
      grepl("^logistic_2", model_metric) ~ "Biomarker logistic",
      grepl("^rf_", model_metric) ~ "Random forest"),
    metric = if_else(
      grepl("_auc$", model_metric),
      "AUROC (higher is better)", "Brier score (lower is better)"),
    model = factor(
      model,
      levels = c("Baseline logistic", "Biomarker logistic", "Random forest")))

performance_means <- performance_plot_data %>%
  group_by(model, metric) %>%
  summarise(mean_value = mean(value), .groups = "drop")

figure4 <- ggplot(performance_plot_data,
  aes(x = model, y = value, colour = model)) +
  geom_point(position = position_jitter(width = 0.08, height = 0),
    size = 2.4, alpha = 0.65) +
  geom_point(data = performance_means, aes(x = model, y = mean_value),
    inherit.aes = FALSE, shape = 23, size = 4, fill = "white", colour = "black") +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_colour_manual(
    values = c("Baseline logistic" = "lightblue",
      "Biomarker logistic" = "royalblue", "Random forest" = "darkblue")) +
  labs(title = "Predictive performance across five outer CV folds",
    subtitle = "Baseline logistic regression, biomarker logistic regression and tuned random forest",
    x = NULL,y = "Performance metric",
    caption = paste(
      "Small points represent held-out fold performance;",
      "diamonds represent the mean across folds.",
      "Random-forest hyperparameters were tuned using inner CV.")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 20,hjust = 1),
    strip.text = element_text(face = "bold"))
figure4
ggsave(
  filename = "figures/figure4_model_comparison.png",
  plot = figure4, width = 11, height = 5.5, dpi = 300)

# Figure 5

calibration_oof <- oof_predictions %>%
  select(SEQN, high_bp, logistic_1_prob, logistic_2_prob, rf_prob) %>%
  pivot_longer(
    cols = c(logistic_1_prob, logistic_2_prob, rf_prob),
    names_to = "model", values_to = "predicted_probability") %>%
  mutate(
    model = recode(model, logistic_1_prob = "Baseline logistic",
      logistic_2_prob = "Biomarker logistic", rf_prob = "Random forest"))

calibration_oof_summary <- calibration_oof %>%
  group_by(model) %>%
  mutate(risk_decile = ntile(predicted_probability, 10)) %>%
  group_by(model, risk_decile) %>%
  summarise(n = n(), mean_predicted = mean(predicted_probability),
    observed = mean(high_bp), .groups = "drop")

figure5 <- ggplot(
  calibration_oof_summary,
  aes(x = mean_predicted, y = observed,
    colour = model, group = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
    colour = "grey", linewidth = 0.7) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c("Baseline logistic" = "lightblue",
      "Biomarker logistic" = "royalblue", "Random forest" = "darkblue")) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1),
    labels = function(x) paste0(round(x * 100), "%")) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1),
    labels = function(x) paste0(round(x * 100), "%")) +
  coord_equal() +
  labs(title = "Out-of-fold calibration of prediction models",
    subtitle = "Five-fold cross-validation | NHANES 2017–2018",
    x = "Mean predicted probability of high BP",
    y = "Observed proportion with high BP",
    colour = "Model",
    caption = paste(
      "Each point represents a model-specific predicted-risk decile.",
      "Dashed line indicates perfect calibration.",
      "Unweighted internal validation.")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), legend.position = "bottom",
    legend.title = element_blank())
figure5

ggsave(filename = "figures/figure5_oof_calibration.png",
  plot = figure5,width = 8, height = 6, dpi = 300)

