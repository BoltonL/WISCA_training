###############################################################################
# WISCA Training Script
# Weighted Incidence Syndromic Combination Antibiogram
#
# This script demonstrates how to calculate empiric antimicrobial coverage
# using the `wisca()` function from the AMR package.
# See https://amr-for-r.org/articles/WISCA.html for more info about WISCA.
#
# The script is intentionally heavily annotated so that participants can
# reproduce the analysis at home after the webinar.
#
# Larisse Bolton, Aislinn Cook, Matthijs Berends
###############################################################################

# Load required packages ------------------------------------------------------

# tidyverse: collection of packages for data manipulation and plotting
# AMR: specialised package for antimicrobial resistance data analysis
library(tidyverse)
library(AMR)

# TIP:
# In RStudio you can open a code outline using:
#   Code > Show Document Outline
# or the shortcut:
#   Ctrl/Cmd + Shift + O
# This makes it easier to navigate through sections of this script.


# HAVEN'T SET UP THE ENVIRONMENT YET? ----------------------------------------

# If you have not prepared the R environment yet (installed packages,
# downloaded the repository data, etc.), please follow the instructions in:
#
# https://github.com/BoltonL/WISCA_training
#
# The README explains how to clone or download the repository and install
# the required packages.


# Import data -----------------------------------------------------------------

# The dataset used in this training is a *cleaned* dummy dataset.
# It contains simulated microbiology data with:
#   - pathogen names
#   - patient age
#   - infection type
#   - antimicrobial susceptibility test results
#
# The dataset is stored inside the repository in:
#   data/wisca_dummy_clean.csv

clean <- read_csv("data/wisca_dummy_clean.csv")


# Always inspect your dataset after importing it.

# glimpse() shows:
#   - column names
#   - column types
#   - example values
glimpse(clean)

# View() opens the dataset in a spreadsheet-like viewer in RStudio
View(clean)


# Prepare data ----------------------------------------------------------------

# The AMR package works with specialised data classes for:
#   - antimicrobial interpretations (S/I/R)
#   - microorganisms
#
# These classes allow AMR functions to understand the data correctly.

clean <- clean |>
  mutate(

    # Convert all antimicrobial test result columns to the "sir" class.
    #
    # The SIR class represents:
    #   S = Susceptible
    #   I = Susceptible, increased exposure
    #   R = Resistant
    #
    # This ensures AMR functions recognise these columns as susceptibility data.
    across(
      c(ampicillin:vancomycin),
      as.sir
    ),

    # Convert the pathogen name column into the AMR "mo" class.
    #
    # This links the organism to the internal microbial taxonomy database
    # included in the AMR package (~79,000 species).
    #
    # This enables functions to understand microbial properties and grouping.
    mo = as.mo(pathogen_name)
  )


# Review the data again to confirm the transformations
glimpse(clean)
View(clean)


# Familiarise with the dataset -----------------------------------------------

# Before running WISCA, it is good practice to explore the data.

# Distribution of isolates by hospital
clean |>
  count(hospital_name) |>
  mutate(p = n / sum(n))   # proportion of total isolates


# Distribution by infection syndrome
clean |>
  count(infection_type) |>
  mutate(p = n / sum(n))


# Age distribution of patients
hist(clean$age_at_onset)


# Most common pathogens
clean |>
  count(mo, sort = TRUE) |>
  mutate(p = n / sum(n))


# Run WISCA model -------------------------------------------------------------

# WISCA = Weighted Incidence Syndromic Combination Antibiogram
#
# Instead of calculating susceptibility per organism, WISCA estimates the
# *empirical coverage* of treatment regimens given the observed distribution
# of pathogens in a syndrome.
#
# In other words:
# "If we treat empirically with this regimen, what proportion of infections
# would be covered?"

# Different ways to specify antimicrobial regimens ----------------------------

# 1. Single antimicrobials
clean |>
  wisca(
    antimicrobials = c(
      "ampicillin",
      "piperacillin_tazobactam",
      "meropenem"
    )
  )


# 2. Explicit combination therapies
clean |>
  wisca(
    antimicrobials = c(
      "ampicillin",
      "ampicillin + meropenem",
      "vancomycin + meropenem"
    )
  )


# 3. Using antimicrobial *selectors*
# These automatically select columns belonging to a drug class
clean |>
  wisca(
    antimicrobials = aminopenicillins()
  )


# 4. Adding a drug to a selector group
clean |>
  wisca(
    antimicrobials = aminopenicillins() + "gentamicin"
  )


# 5. Generating multiple combination regimens
clean |>
  wisca(
    antimicrobials = aminopenicillins() + c("", "amikacin", "gentamicin")
  )


# Run a simple WISCA model ----------------------------------------------------

# simulations = number of Monte Carlo simulations used
# More simulations increase precision but take longer to compute.

wisca_outcome <- clean |>
  wisca(
    antimicrobials = c(
      "ampicillin",
      "piperacillin_tazobactam",
      "meropenem"
    ),
    simulations = 100
  )

wisca_outcome


# Stratified WISCA model ------------------------------------------------------

# WISCA can also be stratified by a syndrome group.
#
# Here we calculate separate empiric coverage estimates for each
# infection type.

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


# Plot the results ------------------------------------------------------------

# Base R plotting
plot(wisca_outcome)

# ggplot2 plotting (usually nicer and easier to customise)
autoplot(wisca_outcome)

# Stratified plot
autoplot(wisca_outcome2)


# Advanced regimen definitions ------------------------------------------------

# You can define multiple realistic empiric treatment strategies.

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


# Stratified advanced model ---------------------------------------------------

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


# Changing language and model parameters -------------------------------------

# AMR supports multiple languages for antimicrobial names and output.

wisca_advanced_inggris <- clean |>
  wisca(
    antimicrobials = c(
      "mero + vanco",
      "meropenem"
    ),

    # Change output language
    language = "Indonesian",

    # Create custom syndrome groups based on patient age
    syndromic_group = case_when(
      clean$age_at_onset == 0 ~ "Bayi",       # babies
      clean$age_at_onset < 18 ~ "Anak-anak",  # children
      .default                = "Dewasa"      # adults
    ),

    # Confidence interval level
    conf_interval = 0.99,

    # Number of simulations
    simulations = 100,

    # Interpretation of S/I/R depends on the guideline used.
    #
    # EUCAST:
    #   S = Susceptible (standard dosing)
    #   I = Susceptible, increased exposure
    #   R = Resistant
    #
    # Therefore EUCAST considers both S and I as clinically treatable.
    #
    # Setting combine_SI = TRUE means:
    #   both S and I will be counted as susceptible when estimating
    #   empirical treatment coverage (EUCAST-style interpretation).
    #
    # CLSI:
    #   S = Susceptible
    #   I = Intermediate
    #   R = Resistant
    #
    # CLSI traditionally treats I more cautiously, so analyses sometimes
    # only count S as susceptible. Set combine = FALSE if using CLSI.
    combine_SI = TRUE
  )

wisca_advanced_inggris
autoplot(wisca_advanced_inggris)


# Access additional model information ----------------------------------------

# The printed output only shows summary results,
# but the WISCA object contains much more information.

wisca_advanced_stratied

# Inspect all attributes stored in the WISCA object
attributes(wisca_advanced_stratied)


# Retrieve parameters used in the model
retrieve_wisca_parameters(wisca_advanced_stratied)


# Export results --------------------------------------------------------------

# The object also contains a long-format numeric table
# that can easily be exported.

attributes(wisca_advanced_stratied)$long_numeric


# Export to Excel
attributes(wisca_advanced_stratied)$long_numeric |>
  writexl::write_xlsx("combinations_stratified.xlsx")

# This file can then be used for:
#   - reports
#   - dashboards
#   - additional statistical analysis
