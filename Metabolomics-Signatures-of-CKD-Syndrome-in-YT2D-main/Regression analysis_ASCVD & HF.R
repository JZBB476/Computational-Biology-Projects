rm(list = ls())
library(tidyverse)
library(broom)
library(haven)
library(AHAprevent)
library(purrr)
setwd("/Users/mac/Desktop/2026 S1/SHBC2026")
raw_data <- read_dta('./nightingale_clean.dta')
raw_data2 <- read_dta('./MODYYT2D060126_Olink_v5_PREVENT_Nightingale.dta')
ascvd_clean <- raw_data2 %>% filter(!is.na(ascvd_std)) %>% 
  select(study_id, aha_prevent_hba1c_10yr_cvd, aha_prevent_hba1c_10yr_ascvd, ascvd_std, age_as_of_recruitment, gender, diabetes_duration, hba1c_value, gfr)
hf_clean <- raw_data2 %>% filter(!is.na(hf_std)) %>% 
  select(study_id, aha_prevent_hba1c_10yr_cvd, aha_prevent_hba1c_10yr_hf, hf_std, age_as_of_recruitment, gender, diabetes_duration, hba1c_value, gfr)

# 1. Clean data: check for non-numeric data and convert "TAG" to NA
vars_to_clean <- c("pyruvate", "bohbutyrate", "creatinine", "gln")
clean_data <- raw_data %>% mutate(across(all_of(vars_to_clean), ~ ifelse(.x == "TAG", NA, .x)))
sum(is.na(clean_data$pyruvate)) #number of NA

# 2. Log-transform and standardize variables (from total_c to s_hdl_tg_pct)
# Adjust the selection range to match your column order
std_data <- clean_data %>% mutate(across(total_c : s_hdl_tg_pct, as.numeric)) %>% 
  mutate(across(total_c:s_hdl_tg_pct, list(ln = ~ log(.x), std = ~ as.vector(scale(log(.x)))), #standardized by subtracting mean and divide by std
                .names = "{.fn}_{.col}"))

# 3. Identify standardized metabolite variables
metabolites <- names(std_data)[startsWith(names(std_data), "std_")]

# 4. Define covariates
covariates <- c("age_as_of_recruitment", "gender", "diabetes_duration", "hba1c_value", "gfr")

# 5. Merge
final_ascvd <- merge(ascvd_clean, std_data, by = "study_id") %>% as_tibble()
final_hf <- merge(hf_clean, std_data, by = "study_id") %>% as_tibble()

# 6. Run linear regressions for each metabolite
##ASCVD
ascvd_result <- map_dfr(metabolites, function(metabolites){
  formula <- as.formula(paste0("ascvd_std ~", metabolites, 
                    "+ age_as_of_recruitment + gender + diabetes_duration + hba1c_value + gfr"))
  
  model_ascvd  <- tryCatch(lm(formula, data = final_ascvd), error = function(e) NULL) #tryCatch catches errors, runs function(e) NULL (where e is the error object), and returns NULL 
  if (is.null(model_ascvd) || nobs(model_ascvd) == 0) {
    return(tibble(metabolite = metabolites, beta = NA_real_, se = NA_real_, pval = NA_real_))
  }
  tidy(model_ascvd) %>% filter(term == metabolites) %>% transmute(metabolite = metabolites, beta = estimate, se = std.error, pval = p.value)
})
ascvd_result <- ascvd_result %>% arrange(desc(beta)) %>% mutate(p.fdr = p.adjust(pval, method = "BH")) %>% 
  mutate(significant = as.integer(p.fdr < 0.05))
ascvd_result %>% summarise(n_failed = sum(beta = is.na(beta)), n_total = n())
#write.csv(ascvd_result, file = "ascvd_result_r.csv")

##HF
hf_result <- map_dfr(metabolites, function(metabolites){
  formula2 <- as.formula(paste0("hf_std ~", metabolites, " + ", paste(covariates, collapse = " + ")))
  
  model_hf <- tryCatch(lm(formula2, data = final_hf), error = function(e){NULL})
  if (is.null(model_hf) || nobs(model_hf) == 0){
    return(tibble(metabolite = metabolites, beta = NA_real_, se = NA_real_, pval = NA_real_))
  }
  tidy(model_hf) %>% filter(term == metabolites) %>% transmute(metabolite = metabolites, beta = estimate, se = std.error, pval = p.value) 
})
hf_result <- hf_result %>% arrange(desc(beta)) %>% mutate(p.fdr = p.adjust(pval, method = "BH")) %>% 
  mutate(significant = as.integer(p.fdr < 0.05))
hf_result %>% summarise(n_failed = sum(beta = is.na(beta)), n_total = n())
#write.csv(hf_result, file = "hf_result_r.csv")

# 7. Summarize significant results
num_sig_ascvd <- sum(ascvd_result$significant == 1, na.rm = TRUE)
cat("Significant metabolites for ASCVD (adjusted):", num_sig_ascvd)

num_sig_hf <- sum(hf_result$significant == 1, na.rm = TRUE)
cat("Significant metabolites for ASCVD (adjusted):", num_sig_hf)

# Print clean table of significant metabolites sorted by FDR p-value
## ascvd_met_results_fdr %>% filter(significant == 1) %>% arrange(fdr_pval) %>% select(metabolite, beta, fdr_pval) %>% print(row.names = FALSE)

# Print clean table of significant metabolites sorted by FDR p-value
## hf_met_results_fdr %>% filter(significant == 1) %>% arrange(fdr_pval) %>% select(metabolite, beta, fdr_pval) %>% print(row.names = FALSE)