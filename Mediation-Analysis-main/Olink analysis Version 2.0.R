rm(list = ls())
set.seed(123)
library(stringr)
library(tidyr)
library(purrr)
library(tidyverse)
library(readxl)
library(arrow)
library(OlinkAnalyze)
library(mediation)
library(broom)
library(BiocManager)

# 1. Read and process data =====================================================
setwd("/Volumes/BURAKI/Olink - Nightingale")
raw_redcap <- read_xlsx('./Data/mody_yt2d_230726-olink.xlsx') %>% filter(blood_status == 1) %>% filter(gad_result_inhouse == 0) %>% 
  filter(pathogenic_identified == 0) %>% mutate(diabetes_duration = age_as_of_recruitment - age_of_onset_years) 

raw_redcap <- raw_redcap %>% mutate(across(c(ethnicity, gender, ace_inhibitor, sglt2_inhibitor, arb), as.factor)) %>%
  mutate(across(c(urine_acr_ug_mg, egfr), as.numeric))%>%
  mutate(log_uACR = log(urine_acr_ug_mg))

raw_olink_ori <- read_parquet('/Volumes/BURAKI/Olink - Nightingale/Data/CT_Explore_384_CAM_INF_528_samples_NPX_2025-12-15.parquet')
raw_olink_check <- OlinkAnalyze::check_npx(df = raw_olink_ori, preferred_names = c(quant = "PCNormalizedNPX"))
olink_clean <- OlinkAnalyze::clean_npx(df = raw_olink_ori, check_log = raw_olink_check, verbose = TRUE)
olink <- olink_clean %>% rename_with(function(x){
  x %>% str_replace_all("([a-z0-9])([A-Z])", "\\1_\\2") %>% str_to_lower()
}) %>% rename(study_id = sample_id) #For all column names, add "_" in between words and change all to lowercase

#inspecting olink data
olink %>% summarise(n_participants = n_distinct(study_id), missing_value = sum(is.na(pcnormalized_npx)), 
                    n_panel = n_distinct(panel), missing_npx = sum(is.na(pcnormalized_npx)), n_olink = n_distinct(olink_id)) 
olink %>% count(sample_type, assay_type, sample_qc, assay_qc, sort = TRUE)

#merge olink with raw_redcap
new_olink <- olink %>% mutate(panel_assay = paste(panel, assay, sep = "_")) %>% 
  pivot_wider(id_cols = study_id, names_from = panel_assay, values_from = pcnormalized_npx)
raw_redcap <- raw_redcap %>% dplyr::select(study_id, homa2ir, log_uACR, age_as_of_recruitment, gender, ethnicity, 
                                    diabetes_duration, bmi, systolic_bp, hba1c_value, ace_inhibitor, sglt2_inhibitor, 
                                    egfr, arb)
analysis_data <- merge(raw_redcap, new_olink, by = "study_id") %>% as_tibble %>% drop_na() #remove all the NAs (317 samples left)

# 2. Analyses ==================================================================
#Identification of proteomic pathways linking insulin resistance to renal function in young-onset type 2 diabetes
#Level 1 of mediation analysis, screening proteins individually

#2.1 Relationship between HOMA2-IR and uACR (Total effect)======================
crude_model <- lm(log_uACR ~ homa2ir, data = analysis_data)
summary(crude_model)
adjusted_model <- lm(log_uACR ~ homa2ir + age_as_of_recruitment + gender + 
                       ethnicity + diabetes_duration + bmi + systolic_bp + hba1c_value + 
                       ace_inhibitor + sglt2_inhibitor + egfr + arb, data = analysis_data)
summary(adjusted_model)
total_effect <- tidy(adjusted_model) %>% filter(term == "homa2ir")
#write.csv(total_effect, file = "Total effect between homa2ir and log_uACR.csv")

#2.2 Relationship between HOMA2-IR and proteins (Path A)========================
protein_names <- analysis_data %>% dplyr::select(starts_with(c("C", "I"))) %>% names()
pathA_results <- map_dfr(protein_names, function(protein){
  formulaA <- as.formula(paste0("`", protein, "` ~ homa2ir + age_as_of_recruitment + gender + 
                                ethnicity + diabetes_duration + bmi + systolic_bp + hba1c_value + 
                                ace_inhibitor + sglt2_inhibitor + egfr + arb")) 
  #backticks used to escape special characters, spaces, or non-standard naming conventions in R variable names
  modelA <- lm(formulaA, data = analysis_data)
  tidy(modelA) %>% filter(term == "homa2ir") %>% mutate(protein = protein) #convert to structured data frame
})
pathA_results <- pathA_results %>% mutate(p_fdr = p.adjust(p.value, method = "BH")) %>% filter(p_fdr < 0.05) #Benjamini-Hochberg FDR
print(pathA_results) #30 significant proteins
#write.csv(pathA_results, file = "Relationship between homa2ir and proteins.csv")

#2.3 Relationship between proteins and uACR (Path B)============================
pathB_results <- map_dfr(protein_names, function(protein){
  formulaB <- paste0("log_uACR ~ `", protein, "` + age_as_of_recruitment + gender + 
                                ethnicity + diabetes_duration + bmi + systolic_bp + hba1c_value + 
                                ace_inhibitor + sglt2_inhibitor + egfr + arb")
  modelB <- lm(formulaB, data = analysis_data)
  tidy(modelB) %>% filter(term == protein) %>% mutate (protein = protein)
})
pathB_results <- pathB_results %>% mutate(p_fdr = p.adjust(p.value, method = "BH")) %>% filter(p_fdr < 0.05)
print(pathB_results) #188 significant proteins
#write.csv(pathB_results, file = "Relationship between proteins and log_uACR.csv")

#2.4 Extract mediator proteins from path A and path B===========================
candidate_mediators <- pathA_results %>% 
  dplyr::select(protein, beta_A = estimate, se_A = std.error, p_A = p.value, fdr_A = p_fdr) %>%
  inner_join(pathB_results %>% dplyr::select(protein, beta_B = estimate, se_B = std.error, p_B = p.value, fdr_B = p_fdr), 
             by = "protein")
candidate_mediators <- candidate_mediators %>% filter(fdr_A < 0.05, fdr_B < 0.05)
print(candidate_mediators)
#write.csv(candidate_mediators, file = "Candidate mediator proteins.csv")

#2.5 Relationship between HOMA2-IR and uACR (Direct effect)=====================
##Adjusted for the 17 candidate mediator proteins
candidate_mediator_protein <- paste0("`", candidate_mediators$protein, "`", collapse = " + ")
direct_effect <- map_dfr(candidate_mediator_protein, function(protein){
  formula_direct_eff <- paste0("log_uACR ~ homa2ir + ", candidate_mediator_protein, " + age_as_of_recruitment + gender + 
                                ethnicity + diabetes_duration + bmi + systolic_bp + hba1c_value + 
                                ace_inhibitor + sglt2_inhibitor + egfr + arb")
  model_direct_eff <- lm(formula_direct_eff, data = analysis_data)
  tidy(model_direct_eff) %>% filter(term == "homa2ir") %>% mutate(protein = protein)
})
direct_effect <- direct_effect %>% mutate(p.value = p.adjust(p.value, method = "BH"))
print(direct_effect)
#write.csv(direct_effect, file = "Direct effect between homa2ir and log_uACR (individual protein).csv")

#2.6 Path A and path B adjusted for proteins other than the candidate mediator proteins
protein_adjustedA_name <- setdiff(pathA_results$protein, candidate_mediators$protein)
protein_adjustedA <- paste0("`", protein_adjustedA_name, "`", collapse = " + ")
pathA_adjusted_results <- map_dfr(setdiff(protein_names, protein_adjustedA_name), function(protein){
  formula_adjustedA <- paste0("`", protein, "` ~ homa2ir + age_as_of_recruitment + gender + 
                                ethnicity + diabetes_duration + bmi + systolic_bp + hba1c_value + 
                                ace_inhibitor + sglt2_inhibitor + egfr + arb + ", protein_adjustedA)
  model_adjustedA <- lm(formula_adjustedA, data = analysis_data)
  tidy(model_adjustedA) %>% filter(term == "homa2ir") %>% mutate(protein = protein)
})
pathA_adjusted_results <- pathA_adjusted_results %>% mutate(p_fdr = p.adjust(p.value, method = "BH")) %>% 
  filter(p_fdr < 0.05)
print(pathA_adjusted_results) #6 signifcant proteins left

protein_adjustedB_name <- setdiff(pathB_results$protein, candidate_mediators$protein)
protein_adjustedB <- paste0("`", protein_adjustedB_name, "`", collapse = " + ")
pathB_adjusted_results <- map_dfr(setdiff(protein_names, protein_adjustedB_name), function(protein){
  formula_adjustedB <- paste0("log_uACR ~ `", protein, "` + age_as_of_recruitment + gender + 
                                ethnicity + diabetes_duration + bmi + systolic_bp + hba1c_value + 
                                ace_inhibitor + sglt2_inhibitor + egfr + arb + ", protein_adjustedB)
  model_adjustedB <- lm(formula_adjustedB, data = analysis_data)
  tidy(model_adjustedB) %>% filter(term == protein) %>% mutate(protein = protein)
})
pathB_adjusted_results <- pathB_adjusted_results %>% mutate(p_fdr = p.adjust(p.value, method = "BH")) %>%
  filter(p_fdr < 0.05)
print(pathB_adjusted_results) #0 signifcant proteins left

#2.7 Mediation proteins and related biological pathway analysis ==========================
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
mediator_protein_genes <- sub(".*_", "", candidate_mediators$protein)
mediator_entrezid <- lapply(mediator_protein_genes, function(x) {
  converted <- bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  unique(converted$ENTREZID)
})
mediator_entrezid <- unique(unlist(mediator_entrezid)) #converts to vector of gene
print(mediator_entrezid)

#Gene Ontology
#Biological Process (BP): Describes the biological goals that a gene product contributes to.
#Molecular Function (MF): Specifies the biochemical activity of the gene product.
#Cellular Component (CC): Indicates where in the cell the gene product functions.
mediator_GO_results <- enrichGO(gene = mediator_entrezid, OrgDb = org.Hs.eg.db,
                                keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05)
dotplot(mediator_GO_results) + ggtitle("Mediator Protein GO - Molecular Function Result") + 
  theme(plot.title = element_text(size = 20))


#2.8 Using STRING database to determine the interaction between mediation proteins ============
# 1. Initialize STRING database interface (v12 is current)
string_db <- STRINGdb$new(version = "12.0", 
                          species = 9606,            # Human NCBI Taxonomy ID
                          score_threshold = 400,     # Medium confidence cutoff (400 = 0.4, 700 = high confidence)
                          input_directory = "")

# 2. Example target gene vector (e.g., differentially expressed genes)
my_genes <- data.frame(gene_name = mediator_protein_genes)

# 3. Map gene symbols to STRING IDs
mapped_genes <- string_db$map(my_genes, "gene_name", removeUnmappedRows = TRUE)

# Basic network plot directly from STRINGdb
string_db$plot_network(mapped_genes$STRING_id)




# 3. High-dimensional Mediation Analysis (HIMA) ==================================
library(HIMA)

#3.1 create data frame for running HIMA
#Part 1 whole olink proteomic data
HIMA_analysis_data <- analysis_data %>% 
  dplyr::select(study_id, homa2ir, log_uACR, age_as_of_recruitment, gender, 
                ethnicity, diabetes_duration, bmi, systolic_bp, hba1c_value, 
                ace_inhibitor, sglt2_inhibitor, egfr, arb)

valid_study_id <- HIMA_analysis_data %>% pull(study_id)
new_olink_df <- semi_join(new_olink, HIMA_analysis_data, by = "study_id") %>% as.data.frame()
rownames(new_olink_df) <- new_olink_df$study_id
new_olink_df$study_id <- NULL

#convert all factor to binary variables
HIMA_analysis_data <- HIMA_analysis_data %>% mutate(gender = ifelse(gender == 1, 1, 0), #male from 1 to 1, female from 2 to 0
                                                    ace_inhibitor = ifelse(ace_inhibitor == 0, 0, 1),
                                                    sglt2_inhibitor = ifelse(sglt2_inhibitor == 0, 0, 1), 
                                                    arb = ifelse(arb == 0, 0, 1))

#dummy expansion on ethnicity due to more than 3 variables 
#(A dummy variable is a numerical variable used in statistical modeling to represent categorical/qualitative data)
ethnicity_dummies <- model.matrix(~ ethnicity, data = analysis_data)[, -1, drop = FALSE] #k-1 = 3 dummy variables are needed here
colnames(ethnicity_dummies) <- gsub("ethnicity", "ethnicity_", colnames(ethnicity_dummies))
HIMA_analysis_data <- HIMA_analysis_data %>% dplyr::select(-ethnicity) %>% bind_cols(as.data.frame(ethnicity_dummies))




hima_fit <- hima(formula = log_uACR ~ homa2ir + age_as_of_recruitment + gender + ethnicity_2 + ethnicity_3 + ethnicity_4 + diabetes_duration + 
                   bmi + systolic_bp + hba1c_value + ace_inhibitor + sglt2_inhibitor + egfr + arb, 
                 data.pheno = HIMA_analysis_data, data.M = new_olink_df, 
                 mediator.type = "gaussian", penalty = "MCP", scale = TRUE, sigcut = 0.05, verbose = TRUE)
summary(hima_fit)


#Part 2 Use only the 17 significant protein mediators identified
new_olink_df2 <- semi_join(new_olink, HIMA_analysis_data, by = "study_id") %>% 
  dplyr::select(study_id, candidate_mediators$protein) %>% as.data.frame()
rownames(new_olink_df2) <- new_olink_df2$study_id
new_olink_df2$study_id <- NULL

hima_fit2 <- hima(formula = log_uACR ~ homa2ir + age_as_of_recruitment + gender + ethnicity_2 + ethnicity_3 + ethnicity_4 + diabetes_duration + 
                   bmi + systolic_bp + hba1c_value + ace_inhibitor + sglt2_inhibitor + egfr + arb, 
                 data.pheno = HIMA_analysis_data, data.M = new_olink_df2, 
                 mediator.type = "gaussian", penalty = "MCP", scale = TRUE, sigcut = 0.05, verbose = TRUE)
summary(hima_fit2)
