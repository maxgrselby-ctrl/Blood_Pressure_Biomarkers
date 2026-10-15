rm(list = ls())
if(!require('pacman'))install.packages('pacman')
pacman::p_load(haven, dplyr, ggplot2, survey, pROC)

# Import required files for data analysis
files <- c("BPX_J", "BPQ_J", "BMX_J", "DEMO_J", "GHB_J", "HDL_J", "HSCRP_J")
base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/"
dir.create("data_raw", showWarnings = FALSE)

for (file in files) {
  destination <- file.path("data_raw", paste0(file, ".XPT"))
  if (!file.exists(destination)) {
    download.file(paste0(base_url, file, ".XPT"), destfile = destination, 
      mode = "wb")}}
bpx <- read_xpt("data_raw/BPX_J.XPT")
bpq <- read_xpt("data_raw/BPQ_J.XPT")
bmx <- read_xpt("data_raw/BMX_J.XPT")
demo <- read_xpt("data_raw/DEMO_J.XPT")
ghb <- read_xpt("data_raw/GHB_J.XPT")
hdl <- read_xpt("data_raw/HDL_J.XPT")
hscrp <- read_xpt("data_raw/HSCRP_J.XPT")

# Select required demographics
demo_processed <- demo %>%
  select(SEQN, RIDSTATR, RIAGENDR, RIDAGEYR,
    RIDRETH3, WTMEC2YR, SDMVPSU, SDMVSTRA)
demo_processed %>%
  summarise(n = n(), missing_gender = sum(is.na(RIAGENDR)), 
    missing_age = sum(is.na(RIDAGEYR)), missing_ethnicity = sum(is.na(RIDRETH3)), 
    missing_weight = sum(is.na(WTMEC2YR)))

# Filter to get only MEC examined adults
demo_processed <- demo_processed %>%
  filter(RIDAGEYR>=20, RIDSTATR == 2, WTMEC2YR > 0)

datasets <- list(bmx = bmx, bpq = bpq, bpx = bpx, 
  demo_processed = demo_processed, ghb = ghb, hdl = hdl, hscrp = hscrp)

# Ensure that all SEQN numbers are unique
sapply(datasets, function(x) {
  c(rows = nrow(x), unique_SEQN = n_distinct(x$SEQN),
    duplicate_SEQN = sum(duplicated(x$SEQN)))
})

# Ensure that variables required for analysis are all present
variables_needed <- list(
  bmx = c("SEQN", "BMXBMI"), 
  bpq = c("SEQN", "BPQ020"),
  bpx = c("SEQN", "BPXSY1", "BPXDI1", "BPXSY2", "BPXDI2", "BPXSY3", "BPXDI3", 
          "BPXSY4", "BPXDI4"),
  demo = c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "WTMEC2YR", "SDMVPSU", 
           "SDMVSTRA"),
  ghb = c("SEQN", "LBXGH"), 
  hdl = c("SEQN", "LBDHDD"), 
  hscrp = c("SEQN", "LBXHSCRP"))

variables_needed$bpq %in% names(bpq)
variables_needed$bmx %in% names(bmx)
variables_needed$bpx %in% names(bpx)
variables_needed$demo %in% names(demo_processed)
variables_needed$ghb %in% names(ghb)
variables_needed$hdl %in% names(hdl)
variables_needed$hscrp %in% names(hscrp)

# Combine datasets by SEQN number
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
  summarise(n = n(), missing_bmi = sum(is.na(BMXBMI)),
    missing_bp_question = sum(is.na(BPQ020)), 
    missing_gender = sum(is.na(RIAGENDR)), missing_age = sum(is.na(RIDAGEYR)),
    missing_ethnicity = sum(is.na(RIDRETH3)),
    missing_weight = sum(is.na(WTMEC2YR)), missing_hba1c = sum(is.na(LBXGH)),
    missing_hdl = sum(is.na(LBDHDD)), missing_crp = sum(is.na(LBXHSCRP)))

# Identify participants with less than 3 s&d bp readings
analysis_raw <- analysis_raw %>%
  mutate(
    valid_pair1 = !is.na(BPXSY1) & !is.na(BPXDI1),
    valid_pair2 = !is.na(BPXSY2) & !is.na(BPXDI2),
    valid_pair3 = !is.na(BPXSY3) & !is.na(BPXDI3),
    valid_pair4 = !is.na(BPXSY4) & !is.na(BPXDI4),
    num_bp_pairs = rowSums(across(c(valid_pair1, valid_pair2, valid_pair3, 
      valid_pair4))),
    valid_bp = num_bp_pairs == 3)

table(analysis_raw$num_bp_pairs)

analysis_raw <- analysis_raw %>%
  mutate(mean_sbp = if_else(valid_bp, 
    rowMeans(across(c(BPXSY1, BPXSY2, BPXSY3, BPXSY4)), na.rm = TRUE), NA_real_),
    mean_dbp = if_else(valid_bp, 
    rowMeans(across(c(BPXDI1, BPXDI2, BPXDI3, BPXDI4)), na.rm = TRUE), NA_real_))

# Remove participants who cannot answer no to bp questionnaire or have less than 
# 3 s&d bp readings
table(analysis_raw$BPQ020, useNA = "ifany")
target_population <- analysis_raw %>%
  filter(BPQ020 == 2, valid_bp)

# Identify missing biomarker or bmi participants
target_population %>%
  summarise(n = n(), missing_bmi = sum(is.na(BMXBMI)),
    missing_hba1c = sum(is.na(LBXGH)), missing_hdl = sum(is.na(LBDHDD)),
    missing_crp = sum(is.na(LBXHSCRP)))

target_population <- target_population %>%
  mutate(complete_predictors = !is.na(BMXBMI) & !is.na(LBXGH) &
      !is.na(LBDHDD) & !is.na(LBXHSCRP))

table(target_population$complete_predictors)

# Analyse bp measurements
target_population %>%
  summarise(
    sbp1_min = min(BPXSY1, na.rm = TRUE), sbp1_max = max(BPXSY1, na.rm = TRUE), 
    sbp2_min = min(BPXSY2, na.rm = TRUE), sbp2_max = max(BPXSY2, na.rm = TRUE),
    sbp3_min = min(BPXSY3, na.rm = TRUE), sbp3_max = max(BPXSY3, na.rm = TRUE),
    sbp4_min = min(BPXSY4, na.rm = TRUE), sbp4_max = max(BPXSY4, na.rm = TRUE),
    
    dbp1_min = min(BPXDI1, na.rm = TRUE), dbp1_max = max(BPXDI1, na.rm = TRUE),
    dbp2_min = min(BPXDI2, na.rm = TRUE), dbp2_max = max(BPXDI2, na.rm = TRUE),
    dbp3_min = min(BPXDI3, na.rm = TRUE), dbp3_max = max(BPXDI3, na.rm = TRUE),
    dbp4_min = min(BPXDI4, na.rm = TRUE), dbp4_max = max(BPXDI4, na.rm = TRUE))

hist(target_population$mean_sbp, breaks = 30,
  main = "Distribution of Mean Systolic Blood Pressure",
  xlab = "Mean systolic blood pressure (mmHg)")

hist(target_population$mean_dbp, breaks = 30,
  main = "Distribution of Mean Diastolic Blood Pressure",
  xlab = "Mean diastolic blood pressure (mmHg)")

# BP data quality assessment
target_population <- target_population %>%
  mutate(high_bp = if_else(mean_sbp >= 140 | mean_dbp >= 90, 1L, 0L))

summary(target_population$mean_sbp)
summary(target_population$mean_dbp)
table(target_population$high_bp)

target_population %>%
  summarise(
    dbp1_zero = sum(BPXDI1 == 0, na.rm = TRUE),
    dbp2_zero = sum(BPXDI2 == 0, na.rm = TRUE),
    dbp3_zero = sum(BPXDI3 == 0, na.rm = TRUE),
    dbp4_zero = sum(BPXDI4 == 0, na.rm = TRUE),
    mean_dbp_zero = sum(mean_dbp == 0))

target_population %>%
  filter(mean_dbp == 0) %>%
  select(SEQN, RIDAGEYR, RIAGENDR, BPXSY1, BPXDI1, BPXSY2, BPXDI2, BPXSY3, 
    BPXDI3, BPXSY4, BPXDI4, mean_sbp, mean_dbp)

# Generate flow diagram
flow_counts <- tibble(
  stage = c("Total NHANES 2017-2018 participants",
    "Adults aged >=20 years and MEC examined",
    "No previously reported hypertension",
    "Three valid paired BP measurements",
    "Complete predictor data"),
  n = c(nrow(demo), nrow(demo_processed),
    analysis_raw %>% filter(BPQ020 == 2) %>% nrow(),
    analysis_raw %>% filter(BPQ020 == 2, valid_bp) %>% nrow(),
    analysis_raw %>% filter(BPQ020 == 2, valid_bp, !is.na(BMXBMI), !is.na(LBXGH),
      !is.na(LBDHDD), !is.na(LBXHSCRP)) %>% nrow()))

# Replace levels with labels
target_population <- target_population %>%
  mutate(high_bp_label = factor(high_bp, levels = c(0, 1), 
    labels = c("No high BP", "Screen-detected high BP")), 
    sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female"))
  )

# Inspect biomarker distribution
hist(target_population$LBXHSCRP, breaks = 40, 
     main = "Distribution of hs-CRP", xlab = "hs-CRP (mg/L)")
hist(target_population$LBXGH, breaks = 40, 
     main = "Distribution of HbA1c", xlab = "HbA1c %")
hist(target_population$LBDHDD, breaks = 40, 
     main = "Distribution of HDL", xlab = "HDL (mg/dL)")

# Continuous variable stats
target_population %>%
  group_by(high_bp_label) %>%
  summarise(n = n(),
    age_mean = mean(RIDAGEYR, na.rm = TRUE), age_sd = sd(RIDAGEYR, na.rm = TRUE),
    bmi_mean = mean(BMXBMI, na.rm = TRUE), bmi_sd = sd(BMXBMI, na.rm = TRUE),
    hba1c_mean = mean(LBXGH, na.rm = TRUE), hba1c_sd = sd(LBXGH, na.rm = TRUE),
    hdl_mean = mean(LBDHDD, na.rm = TRUE), hdl_sd = sd(LBDHDD, na.rm = TRUE),
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

target_population <- target_population %>%
  mutate(age_group = cut(RIDAGEYR, 
      breaks = c(20, 30, 40, 50, 60, 70, Inf), right = FALSE,
      labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70+")))

age_results <- target_population %>%
  group_by(age_group) %>%
  summarise(n = n(), high_bp_n = sum(high_bp), 
  high_bp_percent = 100 * mean(high_bp))
age_results

# Plot bp by age group
ggplot(age_results, aes(x = age_group, y = high_bp_percent)) +
  geom_col() +
  labs(title = "Screen-detected high blood pressure by age",
    subtitle = "Adults without previously reported hypertension, NHANES 2017-2018",
    x = "Age group",
    y = "Participants with high blood pressure (%)") +
  theme_minimal()

# Generate NHANES survey-design object and estimate prevalence
nhanes_design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, 
                  weights = ~WTMEC2YR, nest = TRUE, data = target_population)
nhanes_design

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

# Handle CRP skew through log transformation
target_population <- target_population %>%
  mutate(log2_crp = log2(LBXHSCRP))
hist(target_population$log2_crp, breaks = 40, 
     main = "Distribution of log transformed hs-CRP", xlab = "log2 hs-CRP (mg/L)")

# Generate model dataset
model_data <- target_population %>%
  filter(complete_predictors) %>%
  mutate(decade_age = RIDAGEYR/10)

model_design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, 
                          weights = ~WTMEC2YR, decade_age = RIDAGEYR/10, 
                          nest = TRUE, data = model_data)
model_design

# Model 1 - How well can basic characteristics identify screen-detected high BP?
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

# Investigate high HDL confidence
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

# Investigate potential non-linearity of HDL
model_2_hdl_nonlinear <- svyglm(high_bp ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + I(LBDHDD^2) + log2_crp,
  design = model_design, family = quasibinomial())
summary(model_2_hdl_nonlinear)
regTermTest(model_2_hdl_nonlinear,~ I(LBDHDD^2))

# Investigate whether biomarkers collectively improve model
regTermTest(model_2,~ LBXGH + LBDHDD + log2_crp)

# Generate prediction probabilities
model_data <- model_data %>%
  mutate(pred_model_1 = predict(model_1, type = "response"),
    pred_model_2 = predict(model_2, type = "response"))

summary(model_data$pred_model_1)
summary(model_data$pred_model_2)

# Compare how well each model can rank someone with screen-detected high BP 
# above someone without it
roc_model_1 <- roc(response = model_data$high_bp,predictor = model_data$pred_model_1)
roc_model_2 <- roc(response = model_data$high_bp,predictor = model_data$pred_model_2)
roc_test <- roc.test(roc_model_1, roc_model_2, paired = TRUE,method = "delong")

auc(roc_model_1)
auc(roc_model_2)
roc_test

# Model callibration
brier_model_1 <- mean((model_data$pred_model_1 - model_data$high_bp)^2)
brier_model_2 <- mean((model_data$pred_model_2 - model_data$high_bp)^2)

brier_model_1
brier_model_2

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

# Visualising apparent model callibration
ggplot(calibration_1, aes(x = mean_predicted, y = observed)) +
  geom_point(size = 3, colour = "red") +
  geom_line(colour = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "black") +
  coord_equal() +
  labs(title = "Calibration of baseline model", subtitle = "Age, sex and BMI",
    x = "Mean predicted probability", 
    y = "Observed proportion with screen-detected high BP") +
  theme_minimal()

ggplot(calibration_2, aes(x = mean_predicted, y = observed)) +
  geom_point(size = 3, colour = "blue") +
  geom_line(colour = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "black") +
  coord_equal() +
  labs(title = "Calibration of enhanced model", 
       subtitle = "Age, sex, BMI and biomarkers",
       x = "Mean predicted probability", 
       y = "Observed proportion with screen-detected high BP") +
  theme_minimal()

ggplot() +
  geom_line(data = calibration_1, aes(x = mean_predicted, 
                                      y = observed, colour = "Model 1")) +
  geom_line(data = calibration_2, aes(x = mean_predicted, 
                                      y = observed, colour = "Model 2")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "black") +
  coord_equal() +
  labs(title = "Calibration of models", x = "Mean predicted probability", 
       y = "Observed proportion with screen-detected high BP", colour = "Calibration") +
  scale_colour_manual(values = c("Model 1" = "red", "Model 2" = "blue")) +
  theme_minimal()

# Set up random bootstrap
set.seed(2004)
seq_len(nrow(model_data))
bootstrap_values <- sample(seq_len(nrow(model_data)), size = nrow(model_data),
  replace = TRUE)
bootstrap_data <- model_data[bootstrap_values, ]

# Fit model 1 to bootstrap
bootstrap_model_1 <- glm(high_bp ~ decade_age + sex + BMXBMI,
                         data = bootstrap_data, family = binomial())
summary(bootstrap_model_1)

# Fit model 2 to bootstrap
bootstrap_model_2 <- glm(high_bp ~ decade_age + sex + BMXBMI + LBXGH + LBDHDD + log2_crp,
                         data = bootstrap_data, family = binomial())
summary(bootstrap_model_2)

# Calculate area under curve for BS model 1&2
model_data$pred_logistic_1 <- predict(bootstrap_model_1, newdata = model_data,
  type = "response")
roc_prediction_1 <- roc(response = model_data$high_bp,
  predictor = model_data$pred_logistic_1, quiet = TRUE)
auc(roc_prediction_1)

model_data$pred_logistic_2 <- predict(bootstrap_model_2, newdata = model_data,
  type = "response")
roc_prediction_2 <- roc(response = model_data$high_bp,
  predictor = model_data$pred_logistic_2, quiet = TRUE)
auc(roc_prediction_2)










