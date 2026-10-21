################################################################################
# External validation: NHANES 2017-2018 models tested on NHANES 2015-2016
#
# Research question:
# Do HbA1c, HDL cholesterol and hs-CRP improve prediction of screen-detected
# high BP in an independent NHANES survey cycle?
#
# Study design:
# Fit logistic models using NHANES 2017-2018, then assess their unchanged
# predictions in NHANES 2015-2016. Predictive metrics are unweighted, matching
# the original project's internal cross-validation analysis.
################################################################################

# 1. Clear Environment and load any required packages ##########################

rm(list = ls())
if (!require('pacman')) install.packages('pacman')
pacman::p_load(haven, dplyr, ggplot2, pROC, tidyr)

# 2. Download and import NHANES datasets required for analysis #################

files_2017 <- c('BPX_J', 'BPQ_J', 'BMX_J', 'DEMO_J', 'GHB_J', 'HDL_J', 'HSCRP_J')
files_2015 <- c('BPX_I', 'BPQ_I', 'BMX_I', 'DEMO_I', 'GHB_I', 'HDL_I', 'HSCRP_I')

url_2017 <- 'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/'
url_2015 <- 'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/'
dir.create('data_raw', showWarnings = FALSE)
dir.create('figures', showWarnings = FALSE)

for (file in files_2017) {
  destination <- file.path('data_raw', paste0(file, '.XPT'))
  if (!file.exists(destination)) {
    download.file(paste0(url_2017, file, '.XPT'), destination, mode = 'wb')
  }
}

for (file in files_2015) {
  destination <- file.path('data_raw', paste0(file, '.XPT'))
  if (!file.exists(destination)) {
    download.file(paste0(url_2015, file, '.XPT'), destination, mode = 'wb')
  }
}

bpx_2017 <- read_xpt('data_raw/BPX_J.XPT')
bpq_2017 <- read_xpt('data_raw/BPQ_J.XPT')
bmx_2017 <- read_xpt('data_raw/BMX_J.XPT')
demo_2017 <- read_xpt('data_raw/DEMO_J.XPT')
ghb_2017 <- read_xpt('data_raw/GHB_J.XPT')
hdl_2017 <- read_xpt('data_raw/HDL_J.XPT')
hscrp_2017 <- read_xpt('data_raw/HSCRP_J.XPT')

bpx_2015 <- read_xpt('data_raw/BPX_I.XPT')
bpq_2015 <- read_xpt('data_raw/BPQ_I.XPT')
bmx_2015 <- read_xpt('data_raw/BMX_I.XPT')
demo_2015 <- read_xpt('data_raw/DEMO_I.XPT')
ghb_2015 <- read_xpt('data_raw/GHB_I.XPT')
hdl_2015 <- read_xpt('data_raw/HDL_I.XPT')
hscrp_2015 <- read_xpt('data_raw/HSCRP_I.XPT')

# 3. Define initial populations for analysis ###################################

demo_2017_processed <- demo_2017 %>%
  select(SEQN, RIDSTATR, RIAGENDR, RIDAGEYR, WTMEC2YR) %>%
  filter(RIDAGEYR >= 20, RIDSTATR == 2, WTMEC2YR > 0)

demo_2015_processed <- demo_2015 %>%
  select(SEQN, RIDSTATR, RIAGENDR, RIDAGEYR, WTMEC2YR) %>%
  filter(RIDAGEYR >= 20, RIDSTATR == 2, WTMEC2YR > 0)

# 4. Validate source datasets ##################################################

variables_needed <- list(
  bmx = c('SEQN', 'BMXBMI'),
  bpq = c('SEQN', 'BPQ020'),
  bpx = c('SEQN', paste0('BPXSY', 1:4), paste0('BPXDI', 1:4)),
  demo = c('SEQN', 'RIDAGEYR', 'RIAGENDR', 'WTMEC2YR'),
  ghb = c('SEQN', 'LBXGH'),
  hdl = c('SEQN', 'LBDHDD'),
  hscrp = c('SEQN', 'LBXHSCRP'))

check_datasets <- function(datasets, cycle) {
  for (name in names(variables_needed)) {
    missing_columns <- setdiff(variables_needed[[name]], names(datasets[[name]]))
    if (length(missing_columns) > 0) {
      stop(cycle, ': missing columns in ', name, ': ',
           paste(missing_columns, collapse = ', '))
    }
    if (anyDuplicated(datasets[[name]]$SEQN)) {
      stop(cycle, ': duplicate SEQN values in ', name)
    }
  }
}

datasets_2017 <- list(bmx = bmx_2017, bpq = bpq_2017, bpx = bpx_2017,
                      demo = demo_2017_processed, ghb = ghb_2017, 
                      hdl = hdl_2017, hscrp = hscrp_2017)
datasets_2015 <- list(bmx = bmx_2015, bpq = bpq_2015, bpx = bpx_2015,
                      demo = demo_2015_processed, ghb = ghb_2015, 
                      hdl = hdl_2015, hscrp = hscrp_2015)

check_datasets(datasets_2017, '2017-2018')
check_datasets(datasets_2015, '2015-2016')

# 5. Merge datasets together using SEQN number #################################

merge_cycle <- function(datasets) {
  datasets$demo %>%
    left_join(datasets$bmx %>% select(SEQN, BMXBMI), by = 'SEQN') %>%
    left_join(datasets$bpq %>% select(SEQN, BPQ020), by = 'SEQN') %>%
    left_join(datasets$bpx %>% select(SEQN, all_of(variables_needed$bpx[-1])),
              by = 'SEQN') %>%
    left_join(datasets$ghb %>% select(SEQN, LBXGH), by = 'SEQN') %>%
    left_join(datasets$hdl %>% select(SEQN, LBDHDD), by = 'SEQN') %>%
    left_join(datasets$hscrp %>% select(SEQN, LBXHSCRP), by = 'SEQN')
}

analysis_2017 <- merge_cycle(datasets_2017)
analysis_2015 <- merge_cycle(datasets_2015)

# 6. Assess completeness of BP measurements ####################################

prepare_bp <- function(analysis_raw) {
  analysis_raw <- analysis_raw %>%
    mutate(valid_pair1 = !is.na(BPXSY1) & !is.na(BPXDI1),
           valid_pair2 = !is.na(BPXSY2) & !is.na(BPXDI2),
           valid_pair3 = !is.na(BPXSY3) & !is.na(BPXDI3),
           valid_pair4 = !is.na(BPXSY4) & !is.na(BPXDI4),
           num_bp_pairs = valid_pair1 + valid_pair2 + valid_pair3 + valid_pair4,
           valid_bp = num_bp_pairs == 3)
  
  analysis_raw %>%
    mutate(mean_sbp = if_else(valid_bp, (
      coalesce(if_else(valid_pair1, BPXSY1, NA_real_), 0) +
        coalesce(if_else(valid_pair2, BPXSY2, NA_real_), 0) +
        coalesce(if_else(valid_pair3, BPXSY3, NA_real_), 0) +
        coalesce(if_else(valid_pair4, BPXSY4, NA_real_), 0)) / 3, NA_real_),
      mean_dbp = if_else(valid_bp, (
        coalesce(if_else(valid_pair1, BPXDI1, NA_real_), 0) +
          coalesce(if_else(valid_pair2, BPXDI2, NA_real_), 0) +
          coalesce(if_else(valid_pair3, BPXDI3, NA_real_), 0) +
          coalesce(if_else(valid_pair4, BPXDI4, NA_real_), 0)) / 3, NA_real_))
}

analysis_2017 <- prepare_bp(analysis_2017)
analysis_2015 <- prepare_bp(analysis_2015)

# 7. Define eligible population and screening outcome ##########################

target_2017 <- analysis_2017 %>% filter(BPQ020 == 2, valid_bp)
target_2015 <- analysis_2015 %>% filter(BPQ020 == 2, valid_bp)

target_2017 <- target_2017 %>%
  mutate(high_bp = if_else(mean_sbp >= 140 | mean_dbp >= 90, 1L, 0L),
         sex = factor(RIAGENDR, levels = c(1, 2), labels = c('Male', 'Female')))

target_2015 <- target_2015 %>%
  mutate(high_bp = if_else(mean_sbp >= 140 | mean_dbp >= 90, 1L, 0L),
         sex = factor(RIAGENDR, levels = c(1, 2), labels = c('Male', 'Female')))

# 8. Handle hs-CRP comparability between survey cycles #########################
  # CDC bridging equation: 2015-2016 hs-CRP <= 23 mg/L is converted to the
  # approximate 2017-2018 measurement scale before taking log.
  # Values above 23 mg/L are outside the equation's stated range.
target_2017 <- target_2017 %>% mutate(crp_for_model = LBXHSCRP)
target_2015 <- target_2015 %>%
  mutate(crp_outside_bridge = !is.na(LBXHSCRP) & LBXHSCRP > 23,
         crp_for_model = if_else(!is.na(LBXHSCRP) & LBXHSCRP <= 23,
                                 0.8695 * LBXHSCRP + 0.2954, LBXHSCRP))

# 9. Generate complete-case model datasets #####################################

model_2017 <- target_2017 %>%
  mutate(decade_age = RIDAGEYR / 10,
         log2_crp = if_else(!is.na(crp_for_model) & crp_for_model > 0,
                            log2(crp_for_model), NA_real_)) %>%
  filter(!is.na(decade_age), !is.na(sex), !is.na(BMXBMI),
         !is.na(LBXGH), !is.na(LBDHDD), !is.na(log2_crp))

model_2015_all <- target_2015 %>%
  mutate(decade_age = RIDAGEYR / 10,
         log2_crp = if_else(!is.na(crp_for_model) & crp_for_model > 0,
                            log2(crp_for_model), NA_real_)) %>%
  filter(!is.na(decade_age), !is.na(sex), !is.na(BMXBMI),
         !is.na(LBXGH), !is.na(LBDHDD), !is.na(log2_crp))

model_2015 <- model_2015_all %>% filter(!crp_outside_bridge)

cat('Development sample:', nrow(model_2017), '\n')
cat('Validation sample:', nrow(model_2015), '\n')
cat('2015-2016 excluded above hs-CRP bridge range:',
    sum(model_2015_all$crp_outside_bridge), '\n')
stopifnot(nrow(model_2017) > 0, nrow(model_2015) > 0,
          length(unique(model_2015$high_bp)) == 2)

# 10. Fit original logistic models using ONLY 2017-2018 #########################

original_model_1 <- glm(high_bp ~ decade_age + sex + BMXBMI,
                        data = model_2017, family = binomial())
original_model_2 <- glm(high_bp ~ decade_age + sex + BMXBMI +
                          LBXGH + LBDHDD + log2_crp,
                        data = model_2017, family = binomial())

model_2017 <- model_2017 %>%
  mutate(pred_model_1 = predict(original_model_1, newdata = model_2017,
                                type = 'response'),
         pred_model_2 = predict(original_model_2, newdata = model_2017,
                                type = 'response'))
stopifnot(!anyNA(model_2017$pred_model_1), !anyNA(model_2017$pred_model_2))

# 11. Generate external predictions ############################################

model_2015 <- model_2015 %>%
  mutate(pred_model_1 = predict(original_model_1, newdata = model_2015,
                                type = 'response'),
         pred_model_2 = predict(original_model_2, newdata = model_2015,
                                type = 'response'))
stopifnot(!anyNA(model_2015$pred_model_1), !anyNA(model_2015$pred_model_2))

# 12. Calculate external AUROC and Brier scores #################################

external_roc_1 <- roc(model_2015$high_bp, model_2015$pred_model_1,
                      levels = c(0, 1), direction = '<', quiet = TRUE)
external_roc_2 <- roc(model_2015$high_bp, model_2015$pred_model_2,
                      levels = c(0, 1), direction = '<', quiet = TRUE)

external_auc_1 <- as.numeric(auc(external_roc_1))
external_auc_2 <- as.numeric(auc(external_roc_2))
external_brier_1 <- mean((model_2015$pred_model_1 - model_2015$high_bp)^2)
external_brier_2 <- mean((model_2015$pred_model_2 - model_2015$high_bp)^2)

# 13. Assess calibration intercept and slope ###################################

calibration_stats <- function(observed, predicted) {
  logit_pred <- qlogis(pmin(pmax(predicted, 1e-8), 1 - 1e-8))
  intercept <- coef(glm(observed ~ 1, offset = logit_pred,
                        family = binomial()))[1]
  slope <- coef(glm(observed ~ logit_pred, family = binomial()))[2]
  c(intercept = unname(intercept), slope = unname(slope))
}

calibration_1 <- calibration_stats(model_2015$high_bp, model_2015$pred_model_1)
calibration_2 <- calibration_stats(model_2015$high_bp, model_2015$pred_model_2)

# 14. Summarise external validation results ####################################

external_results <- data.frame(
  model = c('Baseline logistic', 'Biomarker logistic'),
  n = nrow(model_2015),
  events = sum(model_2015$high_bp),
  auroc = c(external_auc_1, external_auc_2),
  brier = c(external_brier_1, external_brier_2),
  calibration_intercept = c(calibration_1['intercept'], calibration_2['intercept']),
  calibration_slope = c(calibration_1['slope'], calibration_2['slope']))

print(external_results)
cat('External AUROC improvement (biomarker - baseline):',
    external_auc_2 - external_auc_1, '\n')

# 15. Merge 2015 and 2017 datasets #############################################

model_2015 <- model_2015 %>%
  mutate(original_df = "2015")
model_2017 <- model_2017 %>%
  mutate(original_df = "2017") %>%
  mutate(crp_outside_bridge = "FILLER")

merged_data <- rbind(model_2015, model_2017)

# 15. Plot external calibration #################################################

calibration_plot_data <- merged_data %>%
  select(high_bp, pred_model_1, pred_model_2, original_df) %>%
  pivot_longer(cols = c(pred_model_1, pred_model_2),
               names_to = 'model', values_to = 'predicted_probability') %>%
  mutate(model = recode(model,
                        pred_model_1 = 'Baseline logistic',
                        pred_model_2 = 'Biomarker logistic')) %>%
  group_by(model, original_df) %>%
  mutate(risk_decile = ntile(predicted_probability, 10)) %>%
  group_by(model, original_df, risk_decile) %>%
  summarise(n = n(), mean_predicted = mean(predicted_probability),
            observed = mean(high_bp), .groups = 'drop') %>%
  mutate(model_cohort = paste(model, "-", original_df))



figure6 <- ggplot(
  calibration_plot_data,
  aes(x = mean_predicted, y = observed,
      colour = model_cohort, group = model_cohort)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey", linewidth = 0.7) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c("Baseline logistic - 2017" = "lightblue", "Biomarker logistic - 2017" = "royalblue",
               "Baseline logistic - 2015" = "red", "Biomarker logistic - 2015" = "darkred")) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1),
                     labels = function(x) paste0(round(x * 100), "%")) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1),
                     labels = function(x) paste0(round(x * 100), "%")) +
  coord_equal() +
  labs(title = 'External calibration of logistic models',
       subtitle = 'NHANES 2015-2016 | Models developed in 2017-2018',
       x = 'Mean predicted probability of high BP',
       y = 'Observed proportion with high BP', colour = 'Model and Cohort',
       caption = 'Unweighted validation; hs-CRP harmonised within bridging range.') +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.title = element_blank())
figure6

ggsave(filename = "figures/figure6_external_calibration.png",
       plot = figure6,width = 8, height = 6, dpi = 300)


