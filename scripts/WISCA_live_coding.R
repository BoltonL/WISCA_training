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
#


# Import data -------------------------------------------------------------

# we use read_csv2 to read a semicolon-separated CSV file that uses a comma as decimal mark
messy <- read_csv2("data/wisca_dummy_clean_vs_messy_200.csv")
View(messy) # review the data


# Clean data --------------------------------------------------------------



# Run WISCA model ---------------------------------------------------------


## Model 1 ----------------------------------------------------------------


## Model 2 ----------------------------------------------------------------
