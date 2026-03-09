
###############################################################################
# WISCA Training Script - Teaching Version with Exercises
# Weighted Incidence Syndromic Combination Antibiogram
#
# This version is designed for teaching and self-study.
# It contains:
#   1. step-by-step explanations
#   2. extra interpretation notes
#   3. small exercises for participants to try on their own
#
# The data used here are fictitious and suitable for training purposes.
#
# Larisse Bolton, Aislinn Cook, Matthijs Berends
###############################################################################

# Load required packages ------------------------------------------------------

# tidyverse:
#   used for importing, cleaning, counting, mutating, and plotting data
#
# AMR:
#   provides domain-specific tools for antimicrobial resistance analysis,
#   including WISCA
library(tidyverse)
library(AMR)

# Helpful RStudio tip:
# Use Code > Show Document Outline
# or Ctrl/Cmd + Shift + O
# to navigate this script more easily.


# Environment setup -----------------------------------------------------------

# If you have not yet set up your R environment, see the repository README:
# https://github.com/BoltonL/WISCA_training
#
# That README explains how to:
#   - install required packages
#   - download the repository
#   - run the training files locally


# Import data -----------------------------------------------------------------

# Read the fictitious cleaned dataset from the repository folder.
# This dataset contains:
#   - patient and hospital context
#   - infection syndrome information
#   - pathogen names
#   - AST results for several antimicrobials
clean <- read_csv("data/wisca_dummy_clean.csv")

# First look at the imported data structure
glimpse(clean)

# Open spreadsheet-style viewer in RStudio
View(clean)

# Why do this?
# Before any modelling, always check:
#   - whether columns were imported correctly
#   - whether character/numeric columns make sense
#   - whether the antimicrobial result columns look like S/I/R data


# Prepare data ----------------------------------------------------------------

# The AMR package uses special classes for several kinds of data.
#
# Two important ones here are:
#   - 'sir' for susceptibility interpretations (S, I, R)
#   - 'mo'  for microorganisms
#
# Converting columns to these classes allows AMR functions to understand them
# correctly and apply the package's microbiological reference knowledge.

clean <- clean |>
  mutate(
    # Convert all AST interpretation columns to the 'sir' class.
    #
    # These columns are assumed to run from 'ampicillin' to 'vancomycin'
    # in the dataset.
    across(
      c(ampicillin:vancomycin),
      as.sir
    ),

    # Convert pathogen names to AMR microorganism codes/class.
    #
    # This allows the package to standardise organism names and use them in
    # downstream AMR analyses.
    mo = as.mo(pathogen_name)
  )

# Review the transformed data
glimpse(clean)
View(clean)

# A good check after this step:
#   - antimicrobial columns should now show class <sir>
#   - the new 'mo' column should be present


# Familiarise yourself with the dataset --------------------------------------

# WISCA estimates empirical treatment coverage, so it is important to first
# understand the composition of the data: which hospitals, which syndromes,
# which patients, and which pathogens are represented.

# 1. How many isolates come from each hospital?
clean |>
  count(hospital_name) |>
  mutate(p = n / sum(n))

# 2. How many isolates are in each infection syndrome?
clean |>
  count(infection_type) |>
  mutate(p = n / sum(n))

# 3. What does the age distribution look like?
hist(clean$age_at_onset)

# 4. Which pathogens are most common?
clean |>
  count(mo, sort = TRUE) |>
  mutate(p = n / sum(n))

# Interpretation note:
# WISCA is weighted by the observed pathogen distribution. So if one organism
# is common in your data, it can strongly influence empiric coverage estimates.


# EXERCISE 1 ------------------------------------------------------------------
# Explore the training data yourself.
#
# Task:
# Create a count table for sex, ward, or another descriptive column if present.
#
# Example pattern:
# clean |> count(your_column_here) |> mutate(p = n / sum(n))
#
# Questions to consider:
#   - Is the dataset balanced across hospitals or syndromes?
#   - Could that affect interpretation of WISCA results?



# WISCA basics ----------------------------------------------------------------

# WISCA stands for:
# Weighted Incidence Syndromic Combination Antibiogram
#
# In practical terms, it asks:
# "Given the pathogens observed in this syndrome, what is the probability
# that a chosen empiric treatment regimen will cover the infection?"
#
# This differs from a conventional antibiogram:
#   - a conventional antibiogram often reports susceptibility per organism
#   - WISCA focuses on empiric regimen performance across the syndrome


# Different ways to specify antimicrobial regimens -----------------------------

# Example 1: single-agent regimens
clean |>
  wisca(
    antimicrobials = c(
      "ampicillin",
      "piperacillin_tazobactam",
      "meropenem"
    )
  )

# Example 2: include explicit combinations
clean |>
  wisca(
    antimicrobials = c(
      "ampicillin",
      "ampicillin + meropenem",
      "vancomycin + meropenem"
    )
  )

# Example 3: use an antimicrobial selector
# 'aminopenicillins()' is an AMR selector that refers to that drug group
clean |>
  wisca(
    antimicrobials = aminopenicillins()
  )

# Example 4: add one extra drug to all members of a selector
clean |>
  wisca(
    antimicrobials = aminopenicillins() + "gentamicin"
  )

# Example 5: automatically create several related combinations
clean |>
  wisca(
    antimicrobials = aminopenicillins() + c("", "amikacin", "gentamicin")
  )

# Teaching note:
# This is useful because in real-world empiric treatment work you often want to
# compare:
#   - monotherapy
#   - one or more combination regimens
#   - groups of related antimicrobials


# EXERCISE 2 ------------------------------------------------------------------
# Try defining your own set of regimens.
#
# Task:
# Run wisca() with a regimen list of your choice.
#
# Ideas:
#   - compare two monotherapies
#   - compare one monotherapy with one combination
#   - compare a selector-based group with an added aminoglycoside
#
# Example starting point:
# clean |>
#   wisca(
#     antimicrobials = c("ampicillin", "ampicillin + gentamicin")
#   )
#
# Questions:
#   - Does adding a second drug appear to improve coverage?
#   - Is the increase clinically large or only modest?



# Run a simple WISCA model ----------------------------------------------------

# We now save the result to an object.
#
# simulations:
#   controls how many simulation runs are used
#   more simulations generally give more stable estimates but take longer
wisca_outcome <- clean |>
  wisca(
    antimicrobials = c(
      "ampicillin",
      "piperacillin_tazobactam",
      "meropenem"
    ),
    simulations = 100
  )

# Print the object
wisca_outcome

# Teaching note:
# During a live webinar, 100 simulations is fast for demonstration.
# In real analysis you may want more simulations for greater stability.


# Plot the simple WISCA model -------------------------------------------------

# Base R plot
plot(wisca_outcome)

# ggplot2-based plot
autoplot(wisca_outcome)

# Interpretation note:
# Look at which regimen appears to provide the highest empirical coverage.
# Also consider uncertainty, not only the point estimate.


# EXERCISE 3 ------------------------------------------------------------------
# Change the number of simulations.
#
# Task:
# Re-run the same model with:
#   - simulations = 50
#   - simulations = 500
#
# Compare the printed results and plots.
#
# Questions:
#   - Do the results change much?
#   - Do they appear more stable with more simulations?



# Stratified WISCA ------------------------------------------------------------

# Often, empiric coverage is not the same across all syndromes.
# We can therefore stratify results by a syndromic grouping variable.

wisca_outcome2 <- clean |>
  wisca(
    antimicrobials = c(
      "ampicillin",
      "piperacillin_tazobactam",
      "meropenem"
    ),
    simulations = 100,
    syndromic_group = "infection_type"
  )

wisca_outcome2
autoplot(wisca_outcome2)

# Interpretation note:
# This lets you compare regimen performance across infection groups.
# A regimen that performs well overall may not be the best for every syndrome.


# EXERCISE 4 ------------------------------------------------------------------
# Compare pooled versus stratified results.
#
# Task:
# Compare:
#   - wisca_outcome
#   - wisca_outcome2
#
# Questions:
#   - Does the best regimen stay the same in every infection type?
#   - Which syndrome seems to have the lowest estimated coverage?
#   - Why might pooled analysis hide useful differences?



# Advanced regimen comparison -------------------------------------------------

# This example uses abbreviated regimen names.
# The AMR package can interpret many antimicrobial names and codes.

wisca_advanced <- clean |>
  wisca(
    antimicrobials = c(
      "AMP + genta",
      "TZP + AMK",
      "mero + vanco",
      "meropenem"
    ),
    simulations = 100
  )

wisca_advanced
autoplot(wisca_advanced)

# Interpretation note:
# This is closer to a real empiric therapy comparison:
# several plausible strategies are evaluated side by side.


# Stratified advanced comparison ----------------------------------------------

wisca_advanced_stratied <- clean |>
  wisca(
    antimicrobials = c(
      "AMP + genta",
      "TZP + AMK",
      "mero + vanco",
      "meropenem"
    ),
    syndromic_group = "infection_type",
    simulations = 100
  )

wisca_advanced_stratied
autoplot(wisca_advanced_stratied)

# Note:
# The object name above is kept as in the original script, although the word
# "stratied" is presumably meant to be "stratified".


# EXERCISE 5 ------------------------------------------------------------------
# Test whether combinations outperform monotherapy.
#
# Task:
# Focus on these regimens:
#   - "meropenem"
#   - "mero + vanco"
#
# Then compare them in:
#   a) an overall model
#   b) a stratified model
#
# Questions:
#   - Does adding vancomycin improve empirical coverage?
#   - Is the effect similar across all syndromic groups?
#   - In which settings might a broader regimen be justified?



# Changing language and grouping ----------------------------------------------

# AMR supports multiple languages for output labelling.
# Here we also define a custom age-based grouping.

wisca_advanced_inggris <- clean |>
  wisca(
    antimicrobials = c(
      "mero + vanco",
      "meropenem"
    ),
    language = "Indonesian",
    syndromic_group = case_when(
      clean$age_at_onset == 0 ~ "Bayi",       # babies
      clean$age_at_onset < 18 ~ "Anak-anak",  # children
      .default                = "Dewasa"      # adults
    ),
    conf_interval = 0.99,
    simulations = 100
  )

wisca_advanced_inggris
autoplot(wisca_advanced_inggris)

# Interpretation note:
# Here we changed three things:
#   - language of the output
#   - grouping variable
#   - confidence interval width
#
# This shows that WISCA can be adapted for teaching, reporting,
# and subgroup analysis.


# EXERCISE 6 ------------------------------------------------------------------
# Create your own custom syndromic grouping.
#
# Task:
# Replace the current age groups with your own grouping rule.
#
# Ideas:
#   - two groups: children vs adults
#   - neonatal vs paediatric vs adult
#   - another categorical variable present in the dataset
#
# Questions:
#   - Do different age groups appear to have different empiric coverage?
#   - Which regimen looks best in each subgroup?



# Access additional information -----------------------------------------------

# A WISCA object contains more than just what is printed to screen.
# You can inspect its attributes for additional stored outputs.

wisca_advanced_stratied
attributes(wisca_advanced_stratied)

# Retrieve model parameters used to create the object
retrieve_wisca_parameters(wisca_advanced_stratied)

# This is useful for:
#   - checking exactly how the model was run
#   - reproducibility
#   - documenting analysis choices


# Export results --------------------------------------------------------------

# Some downstream work is easier outside R, for example:
#   - sharing tables with colleagues
#   - importing into dashboards
#   - building reports in Excel

attributes(wisca_advanced_stratied)$long_numeric

attributes(wisca_advanced_stratied)$long_numeric |>
  writexl::write_xlsx("combinations_stratified.xlsx")

# The exported Excel file contains long-format numeric results that can be
# reused elsewhere.


# Final reflection exercise ---------------------------------------------------

# EXERCISE 7
#
# Write down your conclusion from the training dataset:
#
#   1. Which empiric regimen appears to give the best overall coverage?
#   2. Does that change when stratifying by infection type?
#   3. Does combination therapy always provide a meaningful benefit?
#   4. What limitations should you remember when interpreting WISCA from
#      a fictitious training dataset?
#
# Suggested discussion points:
#   - pathogen mix drives the estimate
#   - subgroup analysis can change conclusions
#   - broader coverage is not automatically better in practice
#   - training data are for learning, not clinical decision-making


###############################################################################
# End of script
###############################################################################
