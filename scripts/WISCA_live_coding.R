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

clean <- read_csv("data/wisca_dummy_clean.csv")

# review the data
glimpse(clean)
View(clean)


# Prepare data ------------------------------------------------------------

clean <- clean |>
  mutate(
    # let's make sure the SIR columns have the 'sir' data type
    across(c(ampicillin:vancomycin),
           as.sir),
    # let's make 'mo' get the right class
    mo = as.mo(pathogen_name)
  )

# review
clean


# Familiarise with data set -----------------------------------------------

clean |> count(hospital_name) |> mutate(p = n / sum(n))

clean |> count(infection_type) |> mutate(p = n / sum(n))

hist(clean$age_at_onset)

clean |> count(mo, sort = TRUE) |> mutate(p = n / sum(n))


# Run WISCA model ---------------------------------------------------------

# now we have various ways:
clean |> wisca(antimicrobials = c("ampicillin", "piperacillin_tazobactam", "meropenem"))
clean |> wisca(antimicrobials = c("ampicillin", "ampicillin + meropenem", "vancomycin + meropenem"))
clean |> wisca(antimicrobials = aminopenicillins())
clean |> wisca(antimicrobials = aminopenicillins() + "gentamicin")
clean |> wisca(antimicrobials = aminopenicillins() + c("", "amikacin", "gentamicin"))

# let's go with this first:
wisca_outcome <- clean |>
  wisca(antimicrobials = c("ampicillin", "piperacillin_tazobactam", "meropenem"),
        simulations = 100)

wisca_outcome2 <- clean |>
  wisca(antimicrobials = c("ampicillin", "piperacillin_tazobactam", "meropenem"),
        simulations = 100,
        syndromic_group = "infection_type")

### Plot the results ----

plot(wisca_outcome)     # base R
autoplot(wisca_outcome) # ggplot2

plot(wisca_outcome2)
autoplot(wisca_outcome2)

### Advanced options ----

wisca_advanced <- clean |>
  wisca(antimicrobials = c("AMP", "TZP", "AMP + GEN", "VAN", "VAN + mero"),
        simulations = 100)

wisca_advanced <- clean |>
  wisca(antimicrobials = c("AMX", "pip/tazo", "AMX + GEN", "VAN", "VAN + mero"),
        language = "Indonesian", # instead of English
        conf_interval = 0.99,    # instead of 0.95
        simulations = 250        # instead of 1000
  )

