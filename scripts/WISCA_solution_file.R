
###############################################################################
# WISCA Training - Short Solution File
#
# Larisse Bolton, Aislinn Cook, Matthijs Berends
###############################################################################

library(tidyverse)
library(AMR)

# Import ----------------------------------------------------------------------

clean <- read_csv("data/wisca_dummy_clean.csv")

# Prepare ---------------------------------------------------------------------

clean <- clean |>
  mutate(
    across(c(ampicillin:vancomycin), as.sir),
    mo = as.mo(pathogen_name)
  )

# Quick exploration -----------------------------------------------------------

clean |> count(hospital_name) |> mutate(p = n / sum(n))
clean |> count(infection_type) |> mutate(p = n / sum(n))
clean |> count(mo, sort = TRUE) |> mutate(p = n / sum(n))
hist(clean$age_at_onset)

# Basic WISCA -----------------------------------------------------------------

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
plot(wisca_outcome)
autoplot(wisca_outcome)

# Stratified by infection type ------------------------------------------------

wisca_outcome2 <- clean |>
  wisca(
    antimicrobials = c(
      "ampicillin",
      "piperacillin_tazobactam",
      "meropenem"
    ),
    syndromic_group = "infection_type",
    simulations = 100
  )

wisca_outcome2
autoplot(wisca_outcome2)

# Advanced regimen comparison -------------------------------------------------

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

# Custom age-based grouping ---------------------------------------------------

wisca_advanced_inggris <- clean |>
  wisca(
    antimicrobials = c("mero + vanco", "meropenem"),
    language = "Indonesian",
    syndromic_group = case_when(
      clean$age_at_onset == 0 ~ "Bayi",
      clean$age_at_onset < 18 ~ "Anak-anak",
      .default                = "Dewasa"
    ),
    conf_interval = 0.99,
    simulations = 100
  )

wisca_advanced_inggris
autoplot(wisca_advanced_inggris)

# Extra output ----------------------------------------------------------------

attributes(wisca_advanced_stratied)
retrieve_wisca_parameters(wisca_advanced_stratied)

# Export ----------------------------------------------------------------------

attributes(wisca_advanced_stratied)$long_numeric |>
  writexl::write_xlsx("combinations_stratified.xlsx")


###############################################################################
# Suggested answers to exercises
#
# 1. Counts/proportions:
#    Use count(...) |> mutate(p = n / sum(n)) for any grouping variable.
#
# 2. Custom regimens:
#    Example:
#    clean |> wisca(antimicrobials = c("ampicillin", "ampicillin + gentamicin"))
#
# 3. Simulations:
#    Re-run with simulations = 50 and simulations = 500 to compare stability.
#
# 4. Pooled vs stratified:
#    Compare wisca_outcome with wisca_outcome2.
#
# 5. Combination vs monotherapy:
#    Compare "meropenem" with "mero + vanco".
#
# 6. Custom groups:
#    Modify syndromic_group with case_when(...) as needed.
###############################################################################
