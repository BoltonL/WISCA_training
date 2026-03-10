# From Data to Decision: Applying WISCA in Sepsis Using R
## Authors: Larisse Bolton, Matthijs Barends, Aislinn Cook  

Welcome! This repository contains all materials for the live training.  This training introduces the **Weighted Incidence Syndromic Combination Antibiogram (WISCA)** approach and demonstrates how to prepare data, run the WISCA function, and interpret outputs in R.

## What You Will Learn

By the end of this webinar, you will be able to:

- Explain what WISCA is and why it improves on traditional antibiograms  
- Prepare sepsis microbiology data for WISCA analysis  
- Run the custom WISCA R function  
- Interpret coverage probabilities and confidence intervals  
- Apply findings to empiric antibiotic decision-making  

## Set Up the Environment

There are various options, assuming you're in RStudio:

1. Click File > New Project > Version Control > Git > use https://github.com/BoltonL/WISCA_training as the Repository URL

2. Or, install `usethis` and run its `use_course()` function to set everything up:

   ```r
   if (!"usethis" %in% rownames(installed.packages())) install.packages("usethis")
   usethis::use_course("https://github.com/BoltonL/WISCA_training/archive/refs/heads/main.zip")
   ```

## Repository Contents

- `data/` → Synthetic sepsis dataset
- `scripts/` → Step-by-step analysis scripts, containing:
  - WISCA live coding file
  - WISCA teaching file with assignments
  - WISCA solution file with teaching solutions
- `output/` → Files generated from the code: plots and Excel files
- `slides/` → Webinar slides and recording (will be available after the webinar)

## Disclaimer

This training dataset is simulated.  
WISCA outputs should be validated with local microbiology and clinical expertise before informing real-world decisions.


We look forward to learning with you!
