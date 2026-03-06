library(tidyverse)
library(AMR)

# Click menu Code > Show Document Outline (Ctrl/Cmd + Shift + O) to show an outline to the right


# HAVEN'T SET UP THE ENVIROMENT YET?
# There are various options, assuming you're now in RStudio:
#
# A. Click File > New Project > Version Control > Git > use https://github.com/BoltonL/WISCA_training as the Repository URL
#
# B. Install `usethis` and run its `use_course()` function to set everything up:
#      - install.packages("usethis")
#      - usethis::use_course("https://github.com/BoltonL/WISCA_training")


# Import data -------------------------------------------------------------

# we use read_csv2 to read a semicolon-separated CSV file that uses a comma as decimal mark
messy <- read_csv2("data/wisca_dummy_clean_vs_messy_200.csv")
View(messy) # review the data


# Clean data --------------------------------------------------------------


# Run WISCA model ---------------------------------------------------------

## Prepare --------------------------------------

antibiotic_in_demo <- c("meropenem","piptaz+amikacin","meropenem+vancomycin")
include_pathogens <- 10
exclude_path_demo <- c("Coagulase-Negative Staphylococcus", "(Unknown Name)","Candida Spp.")
deduped <- "no"
susceptible_I_in <- "S"

## Model 1: manual -------------------------------


## Model 2: AMR package -------------------------

example_isolates_unclean

cleaned <- example_isolates_unclean |>
  mutate(
    # let's make sure the SIR columns have the 'sir' data type
    across(c(AMX, AMC, CIP, GEN),
           as.sir)
  )

# Review
cleaned

cleaned <- cleaned |>
  mutate(
    # let's make the bacteria get the official and right names
    bacteria = as.mo(bacteria)
  )
# 'Microorganism translation was uncertain for four microorganisms'
# do check them!
mo_uncertainties()

# Review again
cleaned

# let's inspect 'hospital'
cleaned |> count(hospital)

# that's enough info for WISCA!

# now we have various ways:
cleaned |> wisca(antimicrobials = c("AMX", "AMC", "CIP"))
cleaned |> wisca(antimicrobials = c("AMX", "AMX + CIP", "AMX + GEN"))
cleaned |> wisca(antimicrobials = aminopenicillins())
cleaned |> wisca(antimicrobials = aminopenicillins() + "GEN")
cleaned |> wisca(antimicrobials = aminopenicillins() + c("", "CIP", "GEN"))

# let's go with this first:
wisca_outcome <- cleaned |>
  wisca(antimicrobials = c("AMX", "CIP", "GEN"))

### Plot the results ----

plot(wisca_outcome)     # base R
autoplot(wisca_outcome) # ggplot2

### Advanced options ----

wisca_outcome <- cleaned |>
  wisca(antimicrobials = c("AMX", "CIP", "GEN"))
