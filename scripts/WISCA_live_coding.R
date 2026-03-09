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

clean <- read_csv("data/data_wisca_in_demo.csv")

# review the data
glimpse(clean)
View(clean)


# Prepare data ------------------------------------------------------------

clean <- clean |>
  mutate(
    # let's make sure the SIR columns have the 'sir' data type
    across(c(MEM:`MEM+VAN`),
           as.sir),
    # let's make 'mo' get the right class
    mo = as.mo(mo)
  )

# review
clean


# Familiarise with data set -----------------------------------------------

clean |> count(SEX) |> mutate(p = n / sum(n))

clean |> count(INFECTION_TYPE) |> mutate(p = n / sum(n))

clean |> count(hospital) |> mutate(p = n / sum(n))

hist(clean$AGE_AT_ONSET)

clean |> count(mo, sort = TRUE) |> mutate(p = n / sum(n))


# Run WISCA model ---------------------------------------------------------

# now we have various ways:
clean |> wisca(antimicrobials = c("AMX", "AMC", "CIP"))
clean |> wisca(antimicrobials = c("AMX", "AMX + CIP", "AMX + GEN"))
clean |> wisca(antimicrobials = aminopenicillins())
clean |> wisca(antimicrobials = aminopenicillins() + "GEN")
clean |> wisca(antimicrobials = aminopenicillins() + c("", "CIP", "GEN"))

# let's go with this first:
wisca_outcome <- cleaned |>
  wisca(antimicrobials = c("AMX", "CIP", "GEN"))

### Plot the results ----

plot(wisca_outcome)     # base R
autoplot(wisca_outcome) # ggplot2

### Advanced options ----

wisca_advanced <- cleaned |>
  wisca(antimicrobials = c("AMX", "AMX + CIP", "AMX + GEN"),
        language = "Indonesian", # instead of English
        conf_interval = 0.99,    # instead of 0.95
        simulations = 250        # instead of 1000
  )

