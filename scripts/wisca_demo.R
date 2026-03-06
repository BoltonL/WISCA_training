library(tidyverse)
library(AMR)
library(lubridate)

##PREPROCESSING
data_dir <- "./Data"
data_name <- "wisca_clean_age_onset.csv"
data_raw <- read.csv(file = file.path(data_dir,data_name), na.strings = c("NA","NULL",".",""))

#let's see how the variables look?
str(data_raw)

#Refactor variables and ensure correct data types
data_in <- data_raw %>%
  mutate(blood_culture_date = ymd(blood_culture_date),
         hospital_name = as.factor(hospital_name),
         age_at_onset = as.numeric(age_at_onset),
         sex = as.factor(sex),
         birthweight_g = as.numeric(birthweight_g),
         infection_type = as.factor(infection_type)) 


antibiotic_in_demo <- c("amp+gent","meropenem","piptaz+amikacin","meropenem+vancomycin") #what antibiotics are you wanting to include in the wisca
exclude_path_demo <- c("Coagulase-Negative Staphylococcus", "(Unknown Name)","Candida Spp.") #any pathogen exclusions as str_to_title
deduped <- "no" #is your dataset deduplicated for first isolates
susceptible_I_in <- "S" #how would you like to code intermediate susceptibility
include_pathogens <- 10

source("./Code/wisca_preprocess_amr.R")  #script to run preprocessing

data_wisca_demo <- wisca_preprocess(x = data_in,
                                   antibiotic_in = antibiotic_in_demo, 
                                   pathogen_in = include_pathogens,
                                   exclude = exclude_path_demo,
                                   susceptible_I = susceptible_I_in,
                                   isolate_first = deduped)

antibiotic_wisca_demo <- data_wisca_demo[[1]]
data_wisca_in_demo <- data_wisca_demo[[2]] 




##RUNNING WISCA
wisca_cover_demo <- data_wisca_in_demo %>%
  group_by(INFECTION_TYPE) %>%
  wisca(#x = data_wisca_in, 
    antimicrobials = antibiotic_wisca_demo,
    col_mo = "mo",
    ab_transform = NULL)

params_wisca_hai <- retrieve_wisca_parameters(wisca_cover_demo)
wisca_cover_demo_plot <- wisca_cover_demo %>%
  pivot_longer(cols = "AMP+GEN":"TZP+AMK", names_to = "keyantimicrobials", values_to = "coverage" ) %>%
  mutate(coverage = str_remove_all(coverage,"\\%")) %>%
  separate(coverage,c("estimate","ci"), sep = "\\(" ) %>%
  mutate(ci = str_remove_all(ci, "\\)")) %>%
  separate(ci,c("l_ci","u_ci"), sep = "-") %>%
  mutate(hospital = str_to_title(INFECTION_TYPE)) %>%
  mutate(hospital = str_replace_all(INFECTION_TYPE,"_"," ")) %>%
  mutate(keyantimicrobials = as.factor(keyantimicrobials),
         hospital = as.factor(hospital)) %>%
  mutate(across(.cols = "estimate":"u_ci", .fns = as.numeric))

ggplot(wisca_cover_demo_plot, aes(x = keyantimicrobials)) + 
  geom_col(aes(y = estimate, fill = INFECTION_TYPE), position = "dodge") +
  geom_errorbar(aes(ymin = l_ci, ymax = u_ci, group = hospital), color = "black",position = position_dodge(width = 0.9), 
                width = 0.1, linewidth = 1) +
  theme_bw() +
  theme(text = element_text(family = "sans", size = 14),
        axis.text = element_text(size = 14)) + 
  labs(title = paste0("Average coverage for top ",include_pathogens," of circulating pathogens"),
       x = "Antibiotic regimens",
       y = "Coverage (95% credible intervals) (%)",
       fill = "Syndrome") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_brewer(palette = 3)
