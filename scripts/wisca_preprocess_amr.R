# Copyright (c) 2026 Larisse Bolton
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

wisca_preprocess <- function(x, antibiotic_in,pathogen_in, exclude, susceptible_I, isolate_first){

  # Define antibiotic regimens under investigation from "antibiotics" dataset
  antibiotic_df <- function(abo_in){
    mono_ab <- str_to_title(abo_in[!str_detect(abo_in,fixed("+"))]) #identify single antibiotics

    mono_antibiotic_x_set <- data.frame()
    for (j in 1:length(mono_ab)){ #generate dataframe for single antibiotics with abbreviations and fullnames
      mono_x_name <- ab_name(mono_ab[j], only_first = TRUE)
      mono_x_ab <- as.ab(mono_x_name)
      mono_antibiotic <- subset(antimicrobials, antimicrobials$ab %in% mono_x_ab)[,c("ab","name")]
      mono_antibiotic_x_set <- rbind.data.frame(mono_antibiotic_x_set,mono_antibiotic, row.names = NULL)
    }
    names(mono_antibiotic_x_set) <- c("ab1","namerx")
    mono_antibiotic_x_set$ab1 <- as.character(mono_antibiotic_x_set$ab1)


    comb_ab <- str_to_title(abo_in[str_detect(abo_in,fixed("+"))]) #identify combination regimens
    if(length(comb_ab) > 0){

      comb_antibiotic_x_total <- data.frame()
      for (jj in 1:length(comb_ab)){ #generate dataframe for combination antibiotics with abbreviations and fullnames
        comb_ab_sep <- unlist(str_split(comb_ab[jj], fixed("+")))
        comb_antibiotic_x_pre <- data.frame(x = rep(0, times = 2), y = rep(0, times = 2))
        for (ii in 1:length(comb_ab_sep)){
          comb_x_name <- ab_name(comb_ab_sep[ii], only_first = TRUE)
          comb_x_ab <- as.ab(comb_x_name)
          comb_antibiotic <- subset(antimicrobials, antimicrobials$ab %in% comb_x_ab)[,c("ab","name")]
          comb_antibiotic_x_pre[ii,] <- comb_antibiotic
          names(comb_antibiotic_x_pre) <- names(comb_antibiotic)
        }
        comb_antibiotic_x <- pivot_wider(comb_antibiotic_x_pre, names_from = name, values_from = ab)
        comb_antibiotic_x_set <- comb_antibiotic_x %>%
          mutate(namerx = str_flatten(str_c(names(comb_antibiotic_x), collapse = "+")))
        new_ab <- sapply(X = as.character(1:(ncol(comb_antibiotic_x_set)-1)), FUN = function(x){str_flatten(str_c("ab",x))})
        names(comb_antibiotic_x_set)[1:(ncol(comb_antibiotic_x_set)-1)] <- new_ab
        names(comb_antibiotic_x_set)[ncol(comb_antibiotic_x_set)]<- "namerx"
        comb_antibiotic_x_total <-rbind.data.frame(comb_antibiotic_x_total, comb_antibiotic_x_set, row.names = NULL)
      }



      antibiotic_x_set <- full_join(x= mono_antibiotic_x_set,y = comb_antibiotic_x_total) #generate full antibiotic dataset for analysis
      antibiotic_x_set <- select(antibiotic_x_set, starts_with("ab"),starts_with("namerx"))

    } else {
      antibiotic_x_set <- mono_antibiotic_x_set
    }
    return(antibiotic_x_set)
  }

  antibiotic_rx <- antibiotic_df(abo_in = antibiotic_in) #full dataframe with antibiotic regimens


  #if the dataset names are in different format to built-in set
  names(x) <- str_to_upper(names(x))
  if ((length(names(x)[sapply(names(x), function(x1){any(str_detect(x1,"[A-Z]$") & !(x1 %in% antimicrobials$ab) & !(x1 %in% antimicrobials$name))})])>0)){
    names(x)[(str_detect(names(x),"^ORG") & str_detect(names(x),"NAME$"))| (str_detect(names(x),"PATH") & str_detect(names(x),"NAME$"))] <- "fullname"
    names(x)[(str_detect(names(x),"^ORG") & str_detect(names(x),"CODE$")) | (str_detect(names(x),"^PATH") & str_detect(names(x),"CODE$"))] <- "organism"
    names(x)[str_detect(names(x),"DATE") & str_detect(names(x),"CULTURE") | str_detect(names(x),"DATE") & str_detect(names(x),"SPECIMEN")] <- "date"
    names(x)[str_detect(names(x),"^EPISODE") | str_detect(names(x),"^PAT") & str_detect(names(x),"ID$")] <- "patient"
    names(x)[which(str_detect(names(x),"^WARD") & str_detect(names(x),"TYPE$"))] <- "ward"
    names(x)[which(str_detect(names(x),"^HOSP"))] <- "hospital"
    names(x)[which(str_detect(names(x),"DATE") & (str_detect(names(x),"REGISTRATION")))] <- "date"

  }

  #remove any duplicated variables
  if (any(duplicated(names(x)))){
    x <- x[,-max(which(names(x) == names(x)[which(duplicated(names(x)))]))]
  } else {
    x <- x
  }

  x <- arrange(x,"date") #arrange in increasing date

  #if your dataset does not have an episode identifier
  if (length(names(x)[sapply(names(x), function(x1){any(str_detect(x1,"patient"))})]) == 0){
    x$patient <- seq(1,nrow(x),1)
  }



  #find demographic data within dataset
  micro_df_ind <- which(sapply(x,function(y) any(str_detect(y,"^S$"))|any(str_detect(y,"^R$"))|any(str_detect(y,"SENSITIVE"))| any(str_detect(y,"RESISTANT"))) ==  TRUE)
  micro_df <-  names(x)[micro_df_ind]#column names for susceptibility data
  demograph_df <- names(x)[which(!names(x) %in% micro_df)] #column names for demographic or organism data

  #Dataset required to contain both mo and fullname
  if ("mo" %in% names(x)){
    # If mo included but not fullname:
    x <- x %>%
      mutate(mo = as.mo(mo))
    test_isolates_pre <- x %>%
      mutate(across(.cols = all_of(demograph_df),~ str_to_title(.x))) %>%
      mutate(mo = toupper(mo))
    if ("fullname" %in% names(x)){
      test_isolates <- test_isolates_pre[,c(demograph_df)]
    } else {
      test_isolates <- left_join_microorganisms(x = test_isolates_pre, by = "mo")[,c(demograph_df,"fullname")] #add in fullname
    }
    test_isolates <- test_isolates %>% #remove bugs to be excluded and final style conversions
      filter(!(fullname %in% exclude)) %>%
      mutate(fullname = str_to_title(fullname))

    rx_abo <- x %>% #from orginal dataset
      filter(mo %in% test_isolates$mo)  %>% #select only those bugs that should be included
      distinct(date, patient, mo, .keep_all = TRUE) #ensure that there are no duplications
    rx_regs <- names(rx_abo)[!(names(rx_abo) %in% demograph_df)] #extract susceptibility data

    if (any(str_detect(rx_regs,fixed("+")))){#if combinations have already been accounted for in susceptibility profiling
      rx_regs_ind <- which(!(names(rx_abo) %in% demograph_df)) # extract susceptibility data
      rx_regs_comb <- rx_regs[str_detect(rx_regs,fixed("+"))]# find combination regimen variables
      rx_regs_comb_x_total <- data.frame()

      for (jj in 1:length(rx_regs_comb)){ # for every combination regimen
        rx_regs_comb_sep <- str_remove_all(unlist(str_split(rx_regs_comb[jj], fixed("+"))),"[[:punct:]]") #extract separate antibiotics included in regimen
        rx_regs_comb_pre <- data.frame(x = rep(0, times = 2), y = rep(0, times = 2))
        for (ii in 1:length(rx_regs_comb_sep)){ #for every antibiotic in the combination regimen
          rx_regs_comb_name <- ab_name(rx_regs_comb_sep[ii], only_first = TRUE) #extract their antibiotic name
          rx_regs_comb_ab <- as.ab(rx_regs_comb_name) #convert antibiotic name to abbreviation
          rx_regs_comb_antibiotic <- subset(antimicrobials, antimicrobials$ab %in% rx_regs_comb_ab)[,c("ab","name")]
          rx_regs_comb_pre[ii,] <- rx_regs_comb_antibiotic
          names(rx_regs_comb_pre) <- names(rx_regs_comb_antibiotic)
        }
        rx_regs_comb_x <- pivot_wider(rx_regs_comb_pre, names_from = name, values_from = ab)
        rx_regs_comb_x_set <- rx_regs_comb_x %>%
          mutate(namerx = str_flatten(str_c(names(rx_regs_comb_x), collapse = "+"))) #generate new antibiotics table
        rx_regs_new_ab <- sapply(X = as.character(1:(ncol(rx_regs_comb_x_set)-1)), FUN = function(x){str_flatten(str_c("ab",x))})
        names(rx_regs_comb_x_set)[1:(ncol(rx_regs_comb_x_set)-1)] <- rx_regs_new_ab
        names(rx_regs_comb_x_set)[ncol(rx_regs_comb_x_set)]<- "namerx"
        rx_regs_comb_x_total <-rbind.data.frame(rx_regs_comb_x_total, rx_regs_comb_x_set, row.names = NULL)
      }

      rx_regs_comb_rename <- select(rx_regs_comb_x_total,starts_with("ab")) #extract all antibiotic abbreviations
      #regenerate variable names for dataset
      suppressWarnings({
        for (k in 1:nrow(rx_regs_comb_x_total)){
          rx_regs_comb_rename_2 <- unlist(rx_regs_comb_rename[k,])
          rx_regs_comb_rename_vec <- str_flatten(str_c(rx_regs_comb_rename_2, collapse = "+"))
          names(rx_abo)[names(rx_abo) == rx_regs_comb[k]] <- rx_regs_comb_rename_vec

        }
      })
      names(rx_abo)[!(str_detect(names(rx_abo),fixed("+")))  & !(names(rx_abo) %in% demograph_df)] <-
        as.ab(names(rx_abo)[!(str_detect(names(rx_abo),fixed("+")))  & !(names(rx_abo) %in% demograph_df)])
      names(rx_abo)[which(duplicated(names(rx_abo)))] <- str_c(names(rx_abo)[which(duplicated(names(rx_abo)))],".x")
      rx_abo_upd <- left_join_microorganisms(x = rx_abo, by = "mo")[,c(names(rx_abo),"fullname")]

    } else {
      names(rx_abo)[!(names(rx_abo) %in% demograph_df)] <- as.ab(names(rx_abo)[!(names(rx_abo) %in% demograph_df)])
      names(rx_abo)[which(duplicated(names(rx_abo)))] <- str_c(names(rx_abo)[which(duplicated(names(rx_abo)))],".x")
      rx_abo_upd <- left_join_microorganisms(x = rx_abo, by = "mo")[,c(names(rx_abo),"fullname")]

    }
  } else {
    #If fullname included in dataset:
    test_isolates_pre <- x %>%
      mutate(across(.cols = all_of(demograph_df),~ str_to_title(.x))) %>%
      filter(!(fullname %in% unique(exclude)))
    print(nrow(test_isolates_pre))
    test_isolates <- left_join_microorganisms(x = test_isolates_pre, by = "fullname")[,c(demograph_df,"mo")]
    test_isolates$fullname <- str_to_title(test_isolates$fullname)
    test_isolates$mo <- toupper(test_isolates$mo)

    rx_abo <- x %>%
      mutate(fullname = str_to_title(fullname)) %>%
      filter(!(fullname %in% unique(exclude)))  %>%
      distinct(date, patient, fullname, .keep_all = TRUE)
    rx_regs <- names(rx_abo)[!(names(rx_abo) %in% demograph_df)] #susceptibility data

    if (any(str_detect(rx_regs,fixed("+")))){#if combinations have already been accounted for
      rx_regs_ind <- which(!(names(rx_abo) %in% demograph_df))
      rx_regs_comb <- rx_regs[str_detect(rx_regs,fixed("+"))]
      rx_regs_comb_x_total <- data.frame()

      for (jj in 1:length(rx_regs_comb)){
        rx_regs_comb_sep <- str_remove_all(unlist(str_split(rx_regs_comb[jj], fixed("+"))),"[[:punct:]]")
        rx_regs_comb_pre <- data.frame(x = rep(0, times = 2), y = rep(0, times = 2))
        for (ii in 1:length(rx_regs_comb_sep)){
          rx_regs_comb_name <- ab_name(rx_regs_comb_sep[ii], only_first = TRUE)
          rx_regs_comb_ab <- as.ab(rx_regs_comb_name)
          rx_regs_comb_antibiotic <- subset(antimicrobials, antimicrobials$ab %in% rx_regs_comb_ab)[,c("ab","name")]
          rx_regs_comb_pre[ii,] <- rx_regs_comb_antibiotic
          names(rx_regs_comb_pre) <- names(rx_regs_comb_antibiotic)
        }
        rx_regs_comb_x <- pivot_wider(rx_regs_comb_pre, names_from = name, values_from = ab)
        rx_regs_comb_x_set <- rx_regs_comb_x %>%
          mutate(namerx = str_flatten(str_c(names(rx_regs_comb_x), collapse = "+")))
        rx_regs_new_ab <- sapply(X = as.character(1:(ncol(rx_regs_comb_x_set)-1)), FUN = function(x){str_flatten(str_c("ab",x))})
        names(rx_regs_comb_x_set)[1:(ncol(rx_regs_comb_x_set)-1)] <- rx_regs_new_ab
        names(rx_regs_comb_x_set)[ncol(rx_regs_comb_x_set)]<- "namerx"
        rx_regs_comb_x_total <-rbind.data.frame(rx_regs_comb_x_total, rx_regs_comb_x_set, row.names = NULL)
      }

      rx_regs_comb_rename <- select(rx_regs_comb_x_total,starts_with("ab"))

      suppressWarnings({
        for (k in 1:nrow(rx_regs_comb_x_total)){
          rx_regs_comb_rename_2 <- unlist(rx_regs_comb_rename[k,])
          rx_regs_comb_rename_vec <- str_flatten(str_c(rx_regs_comb_rename_2, collapse = "+"))
          names(rx_abo)[names(rx_abo) == rx_regs_comb[k]] <- rx_regs_comb_rename_vec
        }
      })
      names(rx_abo)[!(str_detect(names(rx_abo),fixed("+"))) & !(names(rx_abo) %in% demograph_df)] <-
        as.ab(names(rx_abo)[!(str_detect(names(rx_abo),fixed("+"))) & !(names(rx_abo) %in% demograph_df)])
      names(rx_abo)[which(duplicated(names(rx_abo)))] <- str_c(names(rx_abo)[which(duplicated(names(rx_abo)))],".x")
      rx_abo_upd <- left_join_microorganisms(x = rx_abo, by = "fullname")[,c(names(rx_abo),"mo")]

    } else {
      names(rx_abo)[!(names(rx_abo) %in% demograph_df)] <- as.ab(names(rx_abo)[!(names(rx_abo) %in% demograph_df)])
      names(rx_abo)[which(duplicated(names(rx_abo)))] <- str_c(names(rx_abo)[which(duplicated(names(rx_abo)))],".x")
      rx_abo_upd <- left_join_microorganisms(x = rx_abo, by = "fullname")[,c(names(rx_abo),"mo")]

    }
  }


  rx_abo_upd$fullname <- str_to_title(rx_abo_upd$fullname)


  #generate new dataset with demographic information and only antibiotics under investigation
  rx_1 <- select(antibiotic_rx,starts_with("ab")) #all antibiotics
  rx_df <- test_isolates
  var_date <- names(rx_df)[str_detect(names(rx_df),"date")]
  rx_df[[var_date]] <- ymd(rx_df[[var_date]])
  var_patient <- names(rx_df)[str_detect(names(rx_df),"patient")]
  rx_df[[var_patient]] <- toupper(rx_df[[var_patient]])

  names_vec <- vector()
  suppressWarnings({
    for (k in 1:nrow(rx_1)){
      rx_2 <- unlist(rx_1[k,][!is.na(str_extract(rx_1[k,],"[A-Z]+"))])
      if ((length(rx_2) > 1) & !(str_flatten(str_c(rx_2, collapse = "+")) %in% names(rx_abo_upd))){ # if more than one antibiotic in regimen and susceptibility of the combination not considered
        rx_df <- merge(x = rx_df, y = rx_abo_upd[,c("date","patient","mo",rx_2)], by = c("date","patient","mo")) #then bind separate antibiotics in regimen
        #if prior individual antibiotic was added, suffix will be generated
        for (jj in 1:length(rx_2)){
          if (any(str_detect(names(rx_df),paste0(rx_2[jj],".y")))){
            rx_df[[paste0(rx_2[jj],".y")]] <- str_trim(rx_df[[paste0(rx_2[jj],".y")]])
          } else {
            rx_df[[rx_2[jj]]] <- str_trim(rx_df[[rx_2[jj]]])
          }
        }
        names(rx_df)[names(rx_df) %in% rx_2] <- str_flatten(str_c(rx_2, collapse = "+")) #rename both variables as the combination
        names(rx_df)[names(rx_df) %in% paste0(rx_2,".y")] <- str_flatten(str_c(rx_2, collapse = "+"))
        names(rx_df)[names(rx_df) %in% paste0(rx_2,".x")] <- rx_2
        names_vec[k] <- str_flatten(str_c(rx_2, collapse = "+"))
      } else if ((length(rx_2) > 1) & (str_flatten(str_c(rx_2, collapse = "+")) %in% names(rx_abo_upd))){#if more than one antibiotic in regimen and susceptibility of the combination considered
        rx_df <- merge(x = rx_df, y = rx_abo_upd[,c("date","patient","mo",str_flatten(str_c(rx_2, collapse = "+")))], by = c("date","patient","mo"))
        rx_df[[str_flatten(str_c(rx_2, collapse = "+"))]] <- str_trim(rx_df[[str_flatten(str_c(rx_2, collapse = "+"))]])
        names_vec[k] <- str_flatten(str_c(rx_2, collapse = "+"))
      } else {
        rx_df <- merge(x = rx_df, y = rx_abo_upd[,c("date","patient","mo",rx_2)], by = c("date","patient","mo")) #then bind separate antibiotics in regimen
        rx_df[[rx_2]] <- str_trim(rx_df[[rx_2]])
        names(rx_df)[names(rx_df) %in% rx_2] <- str_flatten(str_c(rx_2, collapse = "+")) #rename both variables as the combination
        names(rx_df)[names(rx_df) %in% paste0(rx_2,".y")] <- str_flatten(str_c(rx_2, collapse = "+"))
        names(rx_df)[names(rx_df) %in% paste0(rx_2,".x")] <- rx_2
        names_vec[k] <- rx_2
      }
      if (any(duplicated(names(rx_df)))){
        rx_df <- rx_df[,-max(which(names(rx_df) == names(rx_df)[which(duplicated(names(rx_df)))]))]
      } else {
        rx_df <- rx_df
      }
    }
  })


  if (any(is.na(rx_df$mo))){
    missing1 <- which(is.na(rx_df$mo))
    if (rx_df[missing1,"fullname"] != "NA"){
      rx_df[missing1,"mo"] <- as.mo(rx_df[missing1,]$fullname)
    } else {
      rx_df <- rx_df[-missing1,]
    }
  }

  rx_df <- rx_df %>%
    mutate(across(.cols = everything(),.fns = toupper))

  #Generate adapted first isolate dataset on antimicrobial level
  if (isolate_first == "yes"){ #if dataset has already been deduplicated for first isolates
    test_isolates_first <- rx_df
  } else {
    test_isolates_first <- rx_df[first_isolate(x = rx_df,
                                               method = "e",
                                               episode_days = 14,
                                               col_date = var_date,
                                               col_mo = "mo") == TRUE,] #Only include isolates for microorganisms < 14 days apart
  }


  test_isolates_first_upd <- test_isolates_first %>%
    pivot_longer(cols = names(test_isolates_first)[which(sapply(test_isolates_first,
                                                                function(y) any(str_detect(y,"^S$"))|any(str_detect(y,"^R$"))|any(str_detect(y,"^SENSITIVE$"))| any(str_detect(y,"^RESISTANT$"))))],
                 names_to = "keyantimicrobials", values_to = "antibiogram") %>% #combine all regimens and antibiograms into single variables
    mutate(antibiogram = ifelse(antibiogram == "I"|antibiogram == "INTERMEDIATE",susceptible_I,antibiogram)) %>% #simplify susceptibility outputs - 90-60 rule
    mutate(antibiogram = ifelse(antibiogram == "SENSITIVE","S", antibiogram)) %>%
    mutate(antibiogram = ifelse(antibiogram == "RESISTANT","R",antibiogram)) %>%
    mutate(antibiogram = factor(antibiogram),
           keyantimicrobials = as.factor(keyantimicrobials))


  #recode susceptibility of regimes if not already in the dataset: if any S, then all S; if any NA and not S, then combination NA; if all R, then combination R.
  if (any(duplicated(names(test_isolates_first_upd)))){
    test_isolates_temp3 <- test_isolates_upd %>%
      mutate(antibiogram = as.character(antibiogram)) %>%
      mutate(antibiogram = replace_na(antibiogram,"U")) %>%
      mutate(combiantibiogram = "A")
    test_isolates_temp4 <- test_isolates_temp3 %>%
      group_by(date,patient,mo,keyantimicrobials) %>%
      mutate(combiantibiogram = ifelse(any(str_detect(antibiogram,"S")), "S",combiantibiogram)) %>%
      mutate(combiantibiogram = ifelse(any(str_detect(antibiogram,"U")) & !any(str_detect(antibiogram,"S")), "U",combiantibiogram)) %>%
      mutate(combiantibiogram = ifelse(!any(str_detect(antibiogram,"U")) & !any(str_detect(antibiogram,"S")), "R", combiantibiogram)) %>%
      distinct(date,patient,mo,.keep_all = TRUE) %>%
      ungroup() %>%
      select(-antibiogram) %>%
      rename("antibiogram" = "combiantibiogram") %>%
      mutate(antibiogram = ifelse(antibiogram == "U",NA_character_, antibiogram)) %>%
      mutate(antibiogram = as.factor(antibiogram))

    test_isolates_temp5 <- test_isolates_upd %>% #isolate the single regimens
      group_by(date,patient,mo) %>%
      filter(!str_detect(keyantimicrobials,fixed("+"))) %>%
      ungroup()

    test_isolates_final <- bind_rows(test_isolates_temp4,test_isolates_temp5)
    test_isolates_final <- test_isolates_final %>%
      arrange(patient)

  } else {
    test_isolates_final <- test_isolates_first_upd
  }

  test_isolates_final$fullname <- as.factor(test_isolates_final$fullname)
  test_isolates_final$date <- ymd(test_isolates_final$date)

  #Calculate pathogen incidence and susceptible by regimen
  mo_name_total <- test_isolates_final %>% #deduplicate pathogen list
    filter(!duplicated(mo)) %>%
    select(mo,fullname)

  organism_count_total <- test_isolates_final %>% #count the number of occurrences of pathogens within first isolates set
    group_by(mo) %>%
    distinct(date,patient,.keep_all = TRUE) %>%
    count(mo) %>%
    ungroup() %>%
    arrange(desc(n))

  organism_count_total_upd <- organism_count_total %>%
    mutate(mo = as.character(mo)) %>%
    mutate(mo = ifelse(str_detect(mo, "^B_KLBSL"),"B_KLBSL",mo)) %>%
    mutate(mo = as.factor(mo))
  klebsiella_combine <- organism_count_total_upd %>%
    filter(mo == "B_KLBSL") %>%
    summarise(n_kleb = sum(n))
  organism_count_total_upd <- organism_count_total_upd %>%
    mutate(n = ifelse(mo == "B_KLBSL", klebsiella_combine$n_kleb, n)) %>%
    distinct(mo,.keep_all = TRUE)
  organism_count_total_upd <- rename(organism_count_total_upd,
                                     "n_inc" = "n")

  regimen_organism <- slice_max(organism_count_total_upd, order_by = n_inc, n = pathogen_in) #choose only the top n number of occurring pathogens to consider in WISCA
  total_pathogens <- sum(regimen_organism$n_inc) #total number of occurrences of pathogens
  print(paste0("Total pathogens = ",sum(organism_count_total_upd$n_inc)))
  print(paste0("Number of pathogens contributing to WISCA = ", total_pathogens))
  print(paste0("Percentage of pathogens contributing to WISCA = ",round((total_pathogens/sum(organism_count_total_upd$n_inc))*100,0),"%"))

    test_isolates_final <- test_isolates_final %>%
    mutate(mo = as.character(mo)) %>%
    mutate(mo = ifelse(str_detect(mo, "^B_KLBSL"), "B_KLBSL",mo)) %>% #combine all Klebsiella species together
    mutate(mo = as.factor(mo)) %>%
    mutate(fullname = as.character(fullname)) %>%
    mutate(fullname = ifelse(str_detect(fullname, "Klebsiella"), "Klebsiella sp.", fullname)) %>% #rename all Klebsiella species
    mutate(fullname = as.factor(fullname)) %>%
    filter(mo %in% regimen_organism$mo)  #only keep top 10 pathogens

  test_isolates_final_upd <- test_isolates_final %>%
    mutate(antibiogram = as.sir(antibiogram)) %>%
    pivot_wider(names_from = "keyantimicrobials", values_from = "antibiogram")

  preprocess_out <- list(names_vec,test_isolates_final_upd)
  return(preprocess_out)
}
