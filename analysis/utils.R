
# General utils -----------------------------------------------------------

addSchoolTerm <- function(df,
                          date_col = "date_rdv") {
  df  %>% 
    mutate(term = case_when(
      !!rlang::sym(date_col) < "2021-09-01" ~ "20/21",
      T ~ "21/22")
    )
}

convertDateNum <- function(x) {
  suppressWarnings(as.Date(as.numeric(x)-2, origin = as.Date("1900-01-01")))
}


#' Determine variant
#' @description Determine the variant based on data
#'
#' @param date_pos suspected data of positive result
#'
#' @return a string
#' 
determineVariant <- function(date_pos) {
  alpha_date <- "2021-01-25"
  delta_date <- "2021-06-14"
  omicron_date <- "2021-12-25"
  
  case_when(date_pos < alpha_date ~ "wildtype",
            date_pos < delta_date ~ "alpha",
            date_pos < omicron_date ~ "delta",
            date_pos >= omicron_date ~ "omicron",
            is.na(date_pos) ~ NA_character_,
            T ~ NA_character_)
}

#' Flatten character
#' Removes accents and spetial charcters from string
#' @param chr the character string to flatten
#'
#' @return flattened character string
#' 
flattenChr <- function(chr) {
  stringi::stri_trans_general(str = chr, 
                              id = "Latin-ASCII")
}

#' List to Vector
#'
#' @param s 
#'
#' @return a vector
#'
listToVec <- function(s) {
  s %>% 
    str_split("\n") %>% 
    .[[1]]
}

getIDDict <- function(path = "data/06_sugar_data/20220316-archive_codbars.csv") {
  id_dict <- read_csv(path,
                      col_types =  cols(.default = "c")) %>%
    janitor::clean_names() %>% 
    pivot_longer(cols = contains("code"),
                 names_to = "source",
                 values_to = "codbar") %>% 
    filter(!is.na(codbar), !str_detect(source, "foyer")) %>% 
    select(sugar_id, codbar) 
  
  id_dict
}



#' list2ThisEnv
#' similar to list2env but to current environment
#' 
#' @param l list to pass to environment
#'
#' @return
#' @export
#'
#' @examples
list2ThisEnv <- function(l) {
  for (i in names(l)) {
    eval(parse(text = str_glue("{i} <- l[['{i}']]")))
  }
}

multiPivotLonger <- function(df, 
                             id_cols,
                             col_sets) {
  map(1:length(col_sets), 
      function(x) {
        colset <- col_sets[[x]]
        df %>% 
          select(one_of(id_cols), one_of(colset)) %>% 
          pivot_longer(cols = colset,
                       names_to = "what",
                       values_to = names(col_sets)[x]) %>% 
          mutate(what = str_extract(what, "min|max"))
      }) %>% 
    reduce(.f = inner_join)
}

#' printTicTocTime
#' Based on example in tictoc package help
#' 
#' @param tic 
#' @param toc 
#' @param msg 
#' @param info 
#'
#' @return
#' @export
#'
#' @examples
printTicTocTime <- function(tic, toc, msg, info = "") {
  
  # The time elapsed
  dt <- toc - tic
  
  # Define time units to print
  units <- dplyr::case_when(
    dt > 3600 ~ "hours",    # hours
    dt > 60 ~ "minutes",    # minutes
    TRUE ~ "seconds"        # seconds
  )
  
  # to time units
  dt_units <- dplyr::case_when(
    units == "hours" ~ dt/3600,
    units == "minutes" ~ dt/60,        
    TRUE ~ dt               
  )
  
  if (is.null(msg) || is.na(msg) || length(msg) == 0) {
    outmsg <- paste(round(dt_units, 2), units, "elapsed")
  } else {
    outmsg <- paste0(msg, ": ", round(dt_units, 2), " ", units, " elapsed")
  }
  outmsg
}


#' getSeed
#' Get a random number generator seed value based on time in ms
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
getSeed <- function(x = 0) {
  # Current time in ms, from https://stackoverflow.com/questions/40059573/get-current-time-in-milliseconds
  ms <- round(as.numeric(Sys.time())*1000)
  
  # Remove millions
  seed <- abs(round(ms/1e6)*1e6 - ms) + x
  
  cat(" -- Seed is", seed, "\n")
  
  seed
}

# Questionnaire utils -----------------------------------------------------

cleanSymptomeNames <- function() {
  
}

#' Get columns of interest
#' Get coi based on the sugar id
#' 
#' @param sid the sugar id
#'
#' @return vector of columns of interest
#'
getCOI <- function(sid) {
  # Get data dictionary name
  fname <- getDicNameFromSID(sid)
  
  if (length(fname) == 0) {
    return("all")
  } else {
    res <- readRDS("generated_data/columns_of_interest.rds") %>% 
      .[[fname]]
    if (is.null(res)) {
      return("all")
    }
  }
  res
}

#' get sugar id from data dictionary name
#'
#' @param fname
#'
#' @return the sugar id
#'
getDicNameFromSID <- function(sid) {
  readRDS("generated_data/data_dic_quest_dict.rds") %>% 
    filter(sugar_id == sid) %>% 
    pull(data_dic)
}

#' Get header columns
#' Get the head column names for sugar questionnaires. These column names should alwyas be exported
#'
#' @return
#'
getHeaderCols <- function() {
  c("participant_id",
    "submission_id",
    "codbar",
    "codbar_foyer",
    "testeur",
    "round_code", 
    "date_soumission")
}

getGroup <- function(str) {
  # Masked for data protection
}

#' get sugar id from data dictionnary name
#'
#' @param fname
#'
#' @return the sugar id
#'
getSIDFromDicName <- function(fname) {
  readRDS("generated_data/data_dic_quest_dict.rds") %>% 
    filter(data_dic == fname) %>% 
    pull(sugar_id)
}

getRenamDict <- function() {
  readRDS("generated_data/rename_dictionary.rds")
}

#' Ingest schools questionnaire
#' The function usese the covquest package and select columns of interest
#' 
#' @param path the path to the file to ingest
#'
#' @return data frame with data
#' 
ingestSchoolsQuest <- function(path,
                               keep_multi = T,
                               drop_by_name = NULL,
                               drop_by_regexp = NULL) {
  
  sugar_id <- covquest::getSugarID(path)
  
  covquest::ingestQuestionnaire(path,
                                keep_multi = keep_multi,
                                drop_by_name = drop_by_name,
                                drop_by_regexp = drop_by_regexp) %>% 
    selectCOI(., sugar_id) %>% 
    mutate(source = getQuestID(path)) %>% 
    {
      if("sugar_id" %in% colnames(.)) {
        select(., -sugar_id)
      } else {
        .
      }
    } %>% 
    rename(sugar_id = participant_id)#%>% 
  # renameCOI(., sugar_id)
  
}

mapIngest <- function(files, 
                      keep_multi = T,
                      drop_by_name = NULL,
                      drop_by_regexp = NULL) {
  
  map_df(files, function(x) {
    res <- try(ingestSchoolsQuest(x, 
                                  keep_multi = keep_multi,
                                  drop_by_name = drop_by_name,
                                  drop_by_regexp = drop_by_regexp), silent = F)
    if (!inherits(res, "try-error")) {
      return(res %>% 
               mutate(questionnaire = getQuestID(x)))
    }
  })
}

#' Parse multiple column entries
#' This function parses multiple column entries of the form XXX_1 in questionnaires
#'
#' @param df dataframe with data
#' @param cols columns to parse, will be used with regexp
#' @param header_cols columns to keep from original dataframe
#'
#' @return a long format dataframe with the parsed columns
#'
parseMulti <- function(df, 
                       cols,
                       header_cols = c("sugar_id"),
                       suffix = "") {
  
  # Set row id to keep track of information
  df <- df %>% 
    mutate(row_id = row_number()) 
  
  # Rename columns to match
  df_colnames <- df %>% 
    # select(contains(cols)) %>% 
    colnames()
  
  pattern <- str_glue("_{suffix}[0-9]+")
  indexed_cols <- str_detect(df_colnames, pattern)
  
  # Make sure index appreas at the end
  df_colnames[indexed_cols] <- str_c(
    df_colnames[indexed_cols] %>% str_remove_all(pattern),
    map_chr(df_colnames[indexed_cols], ~str_extract(., pattern)) %>% 
      str_remove_all(str_glue("(?<=_){suffix}(?=[0-9])"))
  )
  
  colnames(df) <- df_colnames
  
  # Get unique types
  all_columns <-  df %>% 
    select(contains(cols)) %>% 
    colnames() %>% 
    str_remove_all("_[0-9]+") %>% 
    unique()
  
  map(all_columns, function(x) {
    
    df %>% 
      select(row_id, one_of(str_c(x, "_", 1:5))) %>%
      pivot_longer(cols = contains(cols),
                   values_to = "value",
                   names_to = "var") %>% 
      filter(!is.na(value)) %>% 
      mutate(index = map_chr(var, ~str_extract(., "(?<=_)[0-9]+")) %>% as.numeric()) %>% 
      rename(!!x := value) %>% 
      select(-var)
  }) %>% 
    reduce(full_join, by = c("row_id", "index")) %>% 
    # add back header columns
    inner_join(
      df %>% 
        select(row_id, one_of(header_cols)),
      ., 
      by = "row_id"
    ) %>% 
    select(-row_id)
}

#' Rename COI
#'
#' @param df 
#' @param sugar_id 
#'
#' @return
#'
renameCOI <- function(df, sugar_id) {
  this_sugar_id <- sugar_id
  
  var_dict <- getRenamDict() %>% 
    filter(sugar_id == this_sugar_id) %>% 
    select(new, old) %>% 
    deframe()
  
  if (length(var_dict) == 0) {
    return(df)
  } else {
    df %>% 
      rename(!!var_dict)
    
  }
}

#' Select columns of interest
#' Select columns of interest in a dataframe based on the sugar id
#'
#' @param df the dataframe to select from
#' @param sugar_id the sugar id
#'
#' @return
#'
selectCOI <- function(df, sugar_id) {
  coi <- getCOI(sugar_id)
  
  if (coi[1] == "all") {
    warning("Could not find data dictionary name for COI for sugar id: ", sugar_id, "\n|--> Keeping all columns")
    df
  } else {
    df %>% 
      select(one_of(getHeaderCols()), one_of(coi))
  }
}

checkDate <- function(x) {
  if (is.na(x)) {
    return(NA_character_)
  } else if (as.character(x) == "77") {
    return(NA_character_)
  } else {
    return(as.character(x))
  }
}

makeDateFromNum <- function(y, m, d) {
  c(y, m, d) %>% 
    map_chr(checkDate) %>% 
    str_c(collapse = "-") %>% 
    as.Date()
}

mapMakeDates <- function(df, col_ids) {
  map(col_ids, function(x) {
    df %>% 
      select(sugar_id, contains(x)) %>% 
      rowwise() %>% 
      mutate(
        "{x}_date" := makeDateFromNum(
          y = .data[[str_c(x, "_y")]],
          m = .data[[str_c(x, "_m")]],
          d = .data[[str_c(x, "_d")]]
        )
      ) 
  }) %>% 
    reduce(.f = inner_join)
}

checkVaccDates <- function(x) {
  n_nonna <- sum(!is.na(x))
  if (n_nonna == 0) {
    date <- as.Date(NA_character_)
    doubt <- FALSE
    bool <- FALSE
  } else {
    freq <- table(x) 
    
    if (length(freq) > 1) {
      if (max(freq) != min(freq)) {
        date <- freq %>% .[which.max(.)]
      } else {
        date <- max(x)
      }
      doubt <- TRUE 
      bool <- TRUE
    } else {
      date <- x[1]
      doubt <- FALSE
      bool <- TRUE
    }
  }
  
  list(date = date, doubt = doubt, bool = bool)
}

checkVaccDatesAcross <- function(x) {
  res <- checkVaccDates(x)
  res <- as_tibble(res) %>% 
    set_colnames(str_c(cur_column(), "_", colnames(.)))
  res
}


# Adults for remove from analysis as these are not in contact with children
getNonEduStaff <- function() {
  c("72045820", "72030220", "72045620", "72066620", "72065510", "72072520") %>% 
    map_chr(~getSIDFromHUGCodbar(.))
}

# Data fields -------------------------------------------------------------
getAllSchools <- function() {
  # This assumes school_groups is defined
  unique(getDataGroups()$school) %>% sort()
}

getAllSchoolsIdentified <- function() {
  c("Mosaic", "La Decouverte", "Lina Stern", "Allobroges", "Magnolias")
}

getAllGroups <- function() {
  # This assumes school_groups is defined
  unique(getDataGroups()$group)  
}

getSchoolTermGroups <- function(sch, term) {
  getDataGroups() %>%
    filter(school %in% sch, term %in% term) %>% 
    distinct(group) %>% 
    pull(group)
}

getAllTerms <- function() {
  # This assumes school_groups is defined
  unique(getDataGroups()$term)  
}

# This assumes that all participant sugar_ids are eithzer in the groups or the 
# households datasets
getAllParticipantIDs <- function() {
  if (exists("all_school_sids")) {
    return(all_school_sids)
  } else {
    all_school_sids <<- unique(c(
      unique(getDataGroups()$sugar_id),
      unique(getDataHouseholds()$sugar_id)
    ))
  }
  all_school_sids
}

getAllOutbreaks <- function() {
  sort(unique(getOutbreakMetadata()$outbreak_id))
}

studyStartDate <- function() {
  as.Date("2021-01-01")
}

studyEndDate <- function() {
  as.Date("2022-08-31")
}

# Data pulling ------------------------------------------------------------

getDataSerologies <- function(data_dir = "data/06_sugar_data/SCS_all_PJ-202308221210/") {
  
  dat_files <- dir(data_dir, full.names = T)
  
  # Start working on serology questionnaires
  sero_files <- dat_files %>% str_subset("sero")
  
  all_sero <- map_df(sero_files, function(x) {
    res <- try(ingestSchoolsQuest(here::here(x), keep_multi = T), silent = T)
    if (!inherits(res, "try-error")) {
      return(res)
    }
  })
  
  all_sero
}


#' buildAdjacency
#'
#' @param full_sugar_ids 
#'
#' @return
#' @export
#'
#' @examples
#' 
buildAdjacencyDF <- function(full_sugar_ids,
                             u_sugar_ids,
                             do_checks = TRUE) {
  
  U <- length(u_sugar_ids)    # number of distinct participants
  u_sugar_ids_num <- makeUIdsNum(u_sugar_ids)
  
  res <- map_df(
    1:nrow(full_sugar_ids), 
    function(i) {
      sid <- full_sugar_ids$sugar_id[i]
      s <- full_sugar_ids$school[i]
      g <- full_sugar_ids$group[i]
      t <- full_sugar_ids$term[i]
      
      # Household adjacency (assumed to be eq. to group)
      in_hh <- full_sugar_ids %>% 
        filter(hh_id %in% full_sugar_ids$hh_id[i])
      
      # Group adjacency
      in_group <- subsetData(full_sugar_ids, 
                             groups = g,
                             schools = s,
                             terms = t,
                             do_checks = do_checks) %>% 
        filter(!(sugar_id %in% in_hh$sugar_id))
      
      
      hh_adj <- tibble(from = full_sugar_ids$sugar_id[i],
                       to = in_hh$sugar_id,
                       type = -1,
                       term = t)
      
      group_adj <- tibble(from = full_sugar_ids$sugar_id[i],
                          to = in_group$sugar_id,
                          type = 0,
                          term = t)
      
      # School adjacency
      in_school <- subsetData(full_sugar_ids, 
                              schools = s,
                              terms = t,
                              do_checks = do_checks) %>% 
        filter(!(sugar_id %in% in_group$sugar_id))
      
      school_adj <- tibble(from = full_sugar_ids$sugar_id[i],
                           to = in_school$sugar_id,
                           type = 1,
                           term = t)
      
      bind_rows(hh_adj,
                group_adj,
                school_adj) %>% 
        filter(from != to)
      
    }) 
  
  res %>% 
    mutate(from = map_chr(from, ~ u_sugar_ids_num[which(u_sugar_ids == .)]),
           to = map_chr(to, ~ u_sugar_ids_num[which(u_sugar_ids == .)])) %>% 
    arrange(from, to, term)
}


# Participant characteristics
cleanGroup <- function(str) {
  # Masked for data protection
}

# Participant characteristics
cleanSchool <- function(str) {
  # Masked for data protection
}

getDataVisitPCRs <- function(data_file = "data/06_sugar_data/scs_visit_results_20220303.csv") {
  read_csv(data_file) %>% 
    select(-1) %>% 
    rename(sugar_id = Contact_Id)
}


getSchoolUMapping <- function() {
  mapping <- str_c("S", LETTERS[1:length(getAllSchoolsIdentified())]) %>% 
    magrittr::set_names(getAllSchoolsIdentified())
  
  mapping
}

getGroupUMapping <- function() {
  mapping <- str_c("G", 1:length(getGroupLevelsIdentified())) %>% 
    magrittr::set_names(getGroupLevelsIdentified())
  
  mapping
}

getDataGroups <- function(data_path = here::here("data/04_pre_processed_data/scs_groups_20231031.csv")) {
  if (exists("school_groups") & Sys.getenv("REDO_DATA") == "FALSE") {
    return(school_groups)
  } else {
    school_groups <<- read_delim(data_path) %>% 
      janitor::clean_names() %>% 
      rename(sugar_id = id, 
             hh_id = code_barre_du_foyer) %>%
      mutate(nom = flattenChr(nom),
             term = str_extract(nom, "[0-9]+/[0-9]+"),
             group = str_remove(nom, term) %>% str_trim(),
             school = cleanSchool(nom),
             school = getSchoolUMapping()[school] %>% as.vector(),
             group = str_remove(group, school) %>% str_trim(),
             group = cleanGroup(group),
             group = getGroupUMapping()[group] %>% as.vector(),
             group = factor(group, levels = getGroupLevels())
      ) %>% 
      filter(!str_detect(nom, "Test")) %>% 
      ungroup()
  }
  return(school_groups)
}

getDataHouseholds <- function(data_path = "data/04_pre_processed_data/scs_w_families-20231109.csv") {
  school_groups <- getDataGroups()
  
  if (exists("households") & Sys.getenv("REDO_DATA") == "FALSE") {
    return(households)
  } else {
    households <<- read_delim(data_path)  %>% 
      janitor::clean_names() %>% 
      rename(sugar_id = id, 
             hh_id = hug_codbar_house,
             sex = hug_sex_contact) %>% 
      mutate(
        # Use Jan 1st 2021 as reference for computing age
        age = difftime("2021-01-01", birthdate, units = "days") %>% {floor(as.numeric(.)/365)},
        age_cat = ifelse(age < 18, "child", "adult"),
        age = ifelse(age_cat == "adult" & birthdate == "1970-01-01", NA, age),
        hug_codbar = as.character(hug_codbar),
        hh_id = as.character(hh_id)
      ) %>% 
      {
        subset <- filter(., sugar_id %in% school_groups$sugar_id) %>% distinct(hh_id)
        inner_join(., subset)
      }
  }
  return(households)
}

#' getTimelineMetadata
#'
#' @param data_path 
#'
#' @return
#' @export
#'
#' @examples
getTimelineMetadata <- function(data_path = here::here("data/03_study_info/SCS_outbreak-overview_javier.xlsx")) {
  # Masked for data protection
}

#' getOutbreakMetadata
#' Get the information on outbreak start and end dates by group
#'
#' @param data_path 
#'
#' @return
#' @export
#'
#' @examples
getOutbreakMetadata <- function(data_path = here::here("data/03_study_info/SCS_outbreak-overview_javier.xlsx")) {
  # Masked for data protection
}

getPartConsentDates <- function(data_path = "data/04_pre_processed_data/scs_consents_date_20231024.csv") {
  dformat <- "%d.%m.%Y"
  
  households <- getDataHouseholds()
  
  if (exists("consent_dates") & Sys.getenv("REDO_DATA") == "FALSE") {
    return(consent_dates)
  } else {
    
    consent_dates <- read_delim(data_path) %>% 
      janitor::clean_names() %>% 
      rename(sugar_id = id) %>%
      # Fix wrong years
      mutate(consentement_school_date = str_replace(consentement_school_date, "2020", "2021")) %>% 
      mutate(consent_date = as.Date(consentement_school_date, format = dformat),
             consent = consentement_school == "Oui") %>% 
      select(sugar_id, consent_date, consent) %>%
      # !! Drop participants without consent
      filter(consent) %>%
      inner_join(
        households %>% 
          select(sugar_id, hug_codbar, hh_id, sex, age, age_cat))
  }
  
  consent_dates
}

addOutbreakID <- function(df, 
                          date_col = "date") {
  
  outbreak_metadata <- getOutbreakMetadata()
  
  res <- df %>% 
    left_join(outbreak_metadata %>% 
                select(outbreak_id, school, contains("visite")) %>% 
                mutate(visite_2 = ifelse(is.na(visite_2), visite_1 + 5, visite_2))) %>% 
    filter(
      (!!rlang::sym(date_col) >= (visite_1 - 5) &
         !!rlang::sym(date_col) <= (visite_3 + 5)) |
        is.na(school) | is.na(sugar_id) | is.na(group)
    ) %>% 
    select(-contains("visite"))
  
  # Second pass
  res2 <- df %>% 
    filter(!(hug_codbar %in% res$hug_codbar)) %>% 
    left_join(outbreak_metadata %>% 
                select(outbreak_id, school, visite_1, visite_2, visite_3) %>% 
                mutate(visite_2 = ifelse(is.na(visite_2), visite_1 + 5, visite_2))) %>% 
    filter(
      (!!rlang::sym(date_col) >= (visite_1 - 5) &
         !!rlang::sym(date_col) <= (visite_2 + 15)) |
        is.na(school) | is.na(sugar_id) | is.na(group)
    ) %>% 
    select(-contains("visite")) %>% 
    add_count(hug_codbar)
  
  # Add back no outbreaks
  no_outbreaks <- df %>% 
    filter(!(seq_id %in% res$seq_id)) 
  
  bind_rows(
    res,
    no_outbreaks
  )
}

#' getDataVaccination
#'
#' @param data_path 
#'
#' @return
#' @export
#'
#' @examples
getDataVaccination <- function(data_path = here::here("generated_data/unified_vaccination_dates.rds")) {
  # Masked for data protection
}


getDataCombined <- function() {
  # Masked for data protection
}

# Get data for a given outbreak
#' Title
#'
#' @param ob_ids 
#'
#' @return
#' @export
#'
#' @examples
getDataOutbreak <- function(ob_ids = getAllOutbreaks()) {
  
  comb_data <- getDataCombined()
  
  # Get outbreak metadata for filtering
  outbreak_metadata <- getOutbreakMetadata() %>% 
    # Here complete visite_2 dates when is NA by replacing with visite_1 + 5 days
    mutate(visite_2 = coalesce(visite_2, visite_1 + 5)) %>% 
    select(outbreak_id, school, group, contains("visite"))
  
  comb_data %>% 
    # Add metadata information
    inner_join(outbreak_metadata %>% 
                 filter(outbreak_id %in% ob_ids)) %>% 
    # Filter for dates compatible with this outbreaks
    filter(date >= (visite_1 - 7),
           date <= (visite_3 + 7))
}


#' getDataStartEnd
#' Get the data for the study start/end timepoints (basline, end of study, and start
#' of school year)
#' 
#' @return
#' @export
#'
#' @examples
getDataStartEnd <- function() {
  
  comb_data <- getDataCombined()
  
  # Get start/end timeline of the study
  timeline_metadata <- getTimelineMetadata()
  
  start_end_data <- map_df(
    unique(timeline_metadata$event), function(x) {
      timeline_metadata %>% 
        filter(event == x) %>% 
        inner_join(comb_data) %>% 
        # Filter for dates compatible with this outbreaks
        filter(date >= (TL - 7),
               date <= (TR + 7))
    })
  
  start_end_data
}

# Get data for a given outbreak
#' Title
#'
#' @param ob_ids 
#'
#' @return
#' @export
#'
#' @examples
getDataBaseline <- function(school, 
                            term) {
  
  comb_data <- getDataCombined()
  
  # Get outbreak metadata for filtering
  outbreak_metadata <- getOutbreakMetadata() %>% 
    # Here complete visite_2 dates when is NA by replacing with visite_1 + 5 days
    mutate(visite_2 = coalesce(visite_2, visite_1 + 5)) %>% 
    select(outbreak_id, school, group, contains("visite"))
  
  comb_data %>% 
    # Add metadata information
    inner_join(outbreak_metadata %>% 
                 filter(outbreak_id %in% ob_ids)) %>% 
    # Filter for dates compatible with this outbreaks
    filter(date >= (visite_1 - 7),
           date <= (visite_3 + 7))
}


#' getBaselineSeroPopData
#' Data from Stringhini et al., 2021. Seroprevalence of anti-SARS-CoV-2 antibodies 
#' after the second pandemic peak. The Lancet Infectious Diseases, 21(5), pp.600-601.
#' https://doi.org/10.1016/S1473-3099(21)00054-2
#' @param data_path 
#'
#' @return
#' @export
#'
#' @examples
getBaselineSeroPopData <- function(data_path = here::here("data/09_other_data/73527f99_HH_silviaAges_results_table.csv")) {
  data <- read_csv(data_path) %>% 
    janitor::clean_names() %>% 
    select(age_cat = category,
           n_tot = obs,
           n_pos = test_positive) %>% 
    mutate(age_L = age_cat %>% str_extract("(?<=\\[)[0-9]+(?=,)") %>% as.numeric(),
           age_R = age_cat %>% str_extract("(?<=,)[0-9]+(?=[\\)|\\]])") %>% as.numeric(),
           age_R = age_R - 1,
           n_pos = str_extract(n_pos, "[0-9]+(?= )") %>% as.numeric()) %>% 
    filter(!is.na(age_L))
  
  data  
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
getPartWithConsent <- function() {
  
  school_groups <- getDataGroups()
  
  # Add information to consent dates
  with_consent <- getPartConsentDates() %>% 
    # !! Keep household data only for participants with consent
    # !! Keep only children and educators in school groups
    inner_join(school_groups %>% 
                 distinct(sugar_id))
  
  with_consent
}


getGroupLevelsIdentified <- function() {
  # Masked for data protection
}

getGroupLevels <- function() {
  getGroupUMapping()[getGroupLevelsIdentified()]
}


getVisitData <- function(outbreak_data, 
                         data,
                         time_lim_lo,
                         time_lim_hi,
                         date_col = "date_rdv") {
  
  date_col_sym <- rlang::sym(date_col)
  
  outbreak_data %>% 
    group_by(school, group, outbreak_id) %>% 
    group_modify(function(x, y) {
      # A. define baseline for each participant
      labeled_dates <- data %>% 
        inner_join(y) %>% 
        # Keep data "close" to visit dates
        filter(
          !!date_col_sym > (min(c(x$visite_1[1], x$visite_2[1]), na.rm = T) - time_lim_lo),
          !!date_col_sym < (x$visite_3[1] + time_lim_hi)
        ) %>%
        # Label dates
        mutate(
          date_label = case_when(
            abs(difftime(x$visite_1[1], !!date_col_sym, "days")) <= 2 ~ "visite_1",
            abs(difftime(x$visite_3[1], !!date_col_sym, "days")) <= 2 ~ "visite_3",
            !!date_col_sym < x$visite_1[1] ~ "prior_visite_1",
            !!date_col_sym > x$visite_3[1] ~ "post_visite_3",
            difftime(x$visite_1[1], !!date_col_sym, "days") < 10 ~ "post_visite_1",
            T ~ "unknown"
          )) %>% 
        rowwise() %>% 
        mutate(
          date_ref_diff = case_when(
            str_detect(date_label, "_1") ~ difftime(!!date_col_sym, x$visite_1[1], "days") %>% as.numeric(),
            str_detect(date_label, "_3") ~ difftime(!!date_col_sym, x$visite_3[1], "days") %>% as.numeric(),
            T ~ NA_real_
          )
        ) %>% 
        ungroup() %>% 
        select(one_of(colnames(y)), sugar_id, !!date_col_sym, date_label, date_ref_diff)
      
      # Select baseline and followup dates
      minmax_dates <- labeled_dates %>% 
        group_by(sugar_id) %>% 
        group_modify(function(.x, .y) {
          
          date_df_min <- getVisitDate(
            .x, 
            visit = "visite_1",
            date_col = date_col) %>% 
            magrittr::set_colnames(str_replace(colnames(.), "visit", "min"))
          
          date_df_max <- getVisitDate(
            .x, 
            "visite_3",
            date_col = date_col)  %>% 
            magrittr::set_colnames(str_replace(colnames(.), "visit", "max"))
          
          bind_cols(
            date_df_min,
            date_df_max
          )
        })
      
      minmax_dates
      
    }) %>% 
    distinct()
}


getVisitDate <- function(.x, 
                         visit,
                         date_col = "date_rdv") {
  
  
  pre_visit <- str_c("prior_", visit)
  post_visit <- str_c("post_", visit)
  
  if (any(.x$date_label == visit)) {
    date_visit <- .x[[date_col]][.x$date_label == visit]
    date_visit_label <- visit
    date_visit_ref_diff <- 0
  } else if (any(.x$date_label == pre_visit | .x$date_label == post_visit)) {
    subx <- filter(.x, date_label == pre_visit | date_label == post_visit)
    ind <- which.min(abs(subx$date_ref_diff))
    date_visit <- subx[[date_col]][ind]
    date_visit_label <- subx$date_label[ind]
    date_visit_ref_diff <- subx$date_ref_diff[ind]
  } else {
    date_visit <- as.Date(NA_character_)
    date_visit_label <- NA_character_
    date_visit_ref_diff <- NA_real_
  }
  
  return(
    tibble(
      date_visit = date_visit,
      date_visit_label = date_visit_label,
      date_visit_ref_diff = date_visit_ref_diff
    )
  )
}

readSympData <- function(x, 
                         dict) {
  dat <- ingestSchoolsQuest(x) %>% 
    janitor::clean_names() %>% 
    select(sugar_id, date_soumission, source, contains("symp"),
           -contains("contact"), -contains("sibling"))
  
  colnames(dat) <- str_replace_all(colnames(dat), "__", "_")
  
  # Rename columns if present
  for (i in 1:length(dict)) {
    old <- dict[i]
    new <- names(dict)[i]
    
    if (old %in% colnames(dat)) {
      dat <- rename(dat, !!dict[i])
    }
  }
  
  dat  
}


getEndlineDates <- function() {
  # Masked for data protection
}

getBaselineDates <- function() {
  # Masked for data protection
}



#' Get number of children/adults in each group
#'
#' @return
#' @export
#'
#' @examples
getSchoolCounts <- function() {
  school_stats_fixed <- readRDS(file = "generated_data/school_event_stats2.rds")
  
  school_stats_fixed %>%
    filter(event_type == "any_test") %>% 
    group_by(school, group, term, who) %>% 
    summarise(n_tot = sum(number)) %>% 
    ungroup() %>% 
    rename(age_cat = who)
}

#' Expand sugar ids
#' Expands sugar ids to create fake ids for group members that were not observed
#'
#' @param obs dataframe as produced by getAllDataGroups, with columns sugar_id, school, group, term
#' @param all_counts dataframe with counts per school/group/term/age_cat as returned by getSchoolCounts
#'
#' @return
#' @export
#'
#' @examples
expandIds <- function(obs,
                      all_counts) {
  
  obs_counts <- obs %>% 
    group_by(school, group, term, age_cat) %>% 
    summarise(n_obs = n()) %>% 
    ungroup()
  
  full_ids <- map_df(1:nrow(obs_counts), function(i) {
    
    these_obs <- obs_counts %>% slice(i)
    these_counts <- all_counts %>% 
      inner_join(these_obs %>% select(-n_obs), by = join_by(school, group, term, age_cat))
    
    res <- obs %>% 
      inner_join(these_obs %>% select(-n_obs), by = join_by(school, group, term, age_cat))
    
    if (these_counts$n_tot[1] <= these_obs$n_obs[1]) {
      cat("-- Matched counts for\t", these_obs %>% unlist(), "target:", these_counts$n_tot[1], "\n")
    } else {
      # Number of missing participants
      n_diff <- these_counts$n_tot[1] - these_obs$n_obs[1]
      cat("-- Missing", n_diff, "for\t","i:", i, "data:", these_obs %>% unlist(), "target:", these_counts$n_tot[1], "\n")
      # Create sequence of unique missing IDS
      hash <- str_c("M",
                    str_sub(these_counts$age_cat[1], 0, 1),
                    digest::digest(these_obs %>% select(-n_obs), algo = "md5") %>% 
                      str_sub(1, 6),
                    1:n_diff, 
                    sep = "-")
      
      res <- bind_rows(
        res,
        tibble(sugar_id = hash) %>% 
          bind_cols(these_obs %>% select(-n_obs))
      )
    }
    
    res
  })
  
  cat("-- Added", nrow(full_ids) - nrow(obs), "participants.\n")
  full_ids
}

checkSchools <- function(schools) {
  if (is.null(schools)) {
    return(schools)
  }
  
  for (i in schools) {
    if (!(i %in% getAllSchools())) {
      stop("School '", i, "' not in available schools: ", str_c(getAllSchools(), collapse = " "))
    }
  }
  schools
}

checkTerms <- function(terms) {
  if (is.null(terms)) {
    return(terms)
  }
  
  for (i in terms) {
    if (!(i %in% getAllTerms())) {
      stop("Term '", i, "' not in available terms: ", str_c(getAllTerms(), collapse = " "))
    }
  }
  terms
}

checkGroups <- function(groups) {
  if (is.null(groups)) {
    return(groups)
  }
  
  for (i in groups) {
    if (!(i %in% getAllGroups())) {
      stop("Group '", i, "' not in available groups: ", str_c(getAllGroups(), collapse = " "))
    }
  }
  groups
}

subsetData <- function(df, 
                       schools = NULL, 
                       groups = NULL, 
                       terms = NULL,
                       verbose = F,
                       do_checks = TRUE) {
  
  if (is.null(schools) & is.null(groups) & is.null(terms) & verbose) {
    cat("-- No data subsetting required \n")
    return(df)
  } 
  
  if (do_checks) {
    # Check validity of choices
    schools <- checkSchools(schools)
    groups <- checkGroups(groups)
    terms <- checkTerms(terms)
  }
  
  if (!is.null(schools)){
    df <- filter(df, school %in% schools)
  } else {
    schools <- "all"
  }
  
  if (!is.null(groups)){
    df <- filter(df, group %in% groups)
  } else {
    groups <- "all"
  }
  
  if (!is.null(terms)){
    df <- filter(df, term %in% terms)
  } else {
    terms <- "all"
  }
  
  filter_str <- str_c("schools: [", str_c(schools, collapse = ","),
                      "], groups: [", str_c(groups, collapse = ","),
                      "], terms: [", str_c(terms, collapse = ","), "]")
  
  if (nrow(df) == 0) {
    stop("No data matched required filters ", filter_str)
  } else {
    if (verbose) {
      cat("-- Subsetting data for", filter_str, "\n")
    }
  }
  
  return(df)
}

makeUIdsNum <- function(x) {
  U <- length(x)    # number of distinct participants
  map_chr(seq_len(U), ~ idToChr(.))
}

computeTerm <- function(date) {
  ifelse(date < getTermSplit(), "20/21", "21/22")
}

getTermSplit <- function() {
  as.Date("2021-09-01")
}

addTerms <- function(df, date_col = "date") {
  df %>% 
    mutate(term = map_chr(!!rlang::sym(date_col), ~ computeTerm(.)))
}

addEpiWeek <- function(df, 
                       date_col = "date") {
  
  df %>% 
    mutate(week = lubridate::epiweek(!!rlang::sym(date_col)),
           year = lubridate::epiyear(!!rlang::sym(date_col)),
           epiweek = str_c(year, week, sep = "-"),
           epiweek_date = as.Date(str_c(epiweek, "-1"), format = "%Y-%U-%u"))
}

getAllPreschools <- function() {
  # Masked for data protection
}

#' makeMaskMat
#' Make matrix of times when individual participants were not in school.
#' This accounts both for detected infection and quaranteen as well as 
#' vacation closure.
#'
#' @param full_sugar_ids 
#' @param test_data
#' @param u_sugar_ids 
#' @param group_closures 
#' @param vacations 
#' @param quarantine_duration
#'
#' @return
#' @export
#'
makeMaskDF <- function(full_sugar_ids,
                       test_data,
                       u_sugar_ids,
                       group_closures = getGroupColsures(),
                       vacations = getVacations(),
                       quarantine_duration = getQuarantinePolicy()) {
  
  
  U <- length(u_sugar_ids)    # number of distinct participants
  u_sugar_ids_num <- makeUIdsNum(u_sugar_ids)
  
  all_masks <- map_df(
    1:nrow(full_sugar_ids),
    function(x) {
      
      g <- full_sugar_ids$group[x]
      s <- full_sugar_ids$sugar_id[x]
      school <- full_sugar_ids$school[x]
      is_preschool <- school %in% getAllPreschools()
      
      # Get group closures
      grp_close <- filter(group_closures, group == g) %>% 
        select(TL, TR)
      
      # Get infection quaranteens
      pos_tests <- filter(test_data, sugar_id == s, n_pos == 1) %>% 
        ungroup() 
      
      if (nrow(pos_tests) > 0) {
        pos_tests <- pos_tests %>% 
          mutate(TL = date,
                 TR = map_chr(date, function(x) {
                   # Get the duration corresponding to the pcr+ date
                   dt <- filter(quarantine_duration, TL <= x, TR >= x) %>% 
                     pull(duration) %>% 
                     first()
                   
                   as.character(x + dt)
                 }) %>% as.Date()) %>% 
          select(TL, TR)
      } else {
        pos_tests <-  tibble(TL = NA, 
                             TR = NA)
      }
      
      # Combine periods
      res <- bind_rows(
        grp_close %>% 
          mutate(type = "group_closure"),
        vacations %>% 
          filter((institution == "preschool") == is_preschool) %>% 
          select(TL, TR) %>% 
          mutate(type = "vacations"), 
        pos_tests %>% 
          mutate(type = "quarantine")) %>% 
        arrange(TL) %>% 
        mutate(row = x,
               sugar_id = s) %>% 
        filter(!is.na(TL))
      
      res
    }) %>% 
    distinct() %>% 
    mutate(id = map_chr(sugar_id, ~u_sugar_ids_num[which(u_sugar_ids == .)])) %>% 
    arrange(id, TL) %>% 
    group_by(id) %>% 
    mutate(num = row_number()) %>% 
    ungroup() 
  
  
  all_masks
}


makeMaskInfo <- function(masks_df) {
  
  masks_df2 <- masks_df %>% 
    arrange(as.numeric(id), num) %>% 
    # Dates to time since ref date
    mutate(TL = map_dbl(TL, ~dateToTime(.)),
           TR = map_dbl(TR, ~dateToTime(.))) %>% 
    distinct(id, TL, TR) %>%
    group_by(id) %>% 
    mutate(num = row_number()) %>% 
    ungroup()
  
  n_masks <- masks_df2 %>% 
    distinct(id, TL, TR) %>% 
    arrange(as.numeric(id), TL) %>% 
    count(id) %>% 
    pull(n)
  
  start_mask <- masks_df2 %>% 
    select(id, TL, num) %>% 
    pivot_wider(values_from = "TL",
                names_from = "num") %>% 
    select(-id) %>% 
    as.matrix()
  
  start_mask[is.na(start_mask)] <- 0
  
  
  end_mask <- masks_df2 %>% 
    select(id, TR, num) %>% 
    pivot_wider(values_from = "TR",
                names_from = "num") %>% 
    select(-id) %>% 
    as.matrix()
  
  end_mask[is.na(end_mask)] <- 0
  
  list(
    n_masks = n_masks,
    start_mask = start_mask,
    end_mask = end_mask
  )
}


#' getNextTermGroup
#' Get the school group in the next term accounting for aging
#' 
#' @param group 
#'
#' @return
#' @export
#'
#' @examples
getNextTermGroup <- function(group) {
  # Masked for data protection
}

#' getPerviousGroup
#' Get school group in the previous term
#' @param group 
#'
#' @return
#' @export
#'
#' @examples
getPerviousGroup <- function(group) {
  res <- getAllGroups()[which(map_chr(getAllGroups(), ~ getNextTermGroup(.)) == group)]
  if (length(res) == 0) {
    "first"
  } else {
    first(as.character(res))
  }
}

correctOutbreakID <- function(df) {
  df %>% 
    mutate(outbreak_id = as.numeric(outbreak_id),
           outbreak_id = case_when(outbreak_id == 6 ~ 8,
                                   outbreak_id == 7 ~ 6,
                                   outbreak_id == 8 ~ 7,
                                   TRUE ~ outbreak_id),
           outbreak_id = factor(outbreak_id))
}


getSchoolCorrectedCounts <- function() {
  # Masked for data protection
}

# Network utils -----------------------------------------------------------
filterVertices <- function(vertices, 
                           schools = getAllSchools(),
                           groups = getAllGroups(),
                           terms = getAllTerms(),
                           keep_hh = TRUE) {
  # Keep all vertices
  subvert <- vertices %>% 
    filter(school %in% schools | is.na(school),
           group %in% groups | is.na(group),
           term %in% terms | is.na(term)) %>% 
    select(id, school, term, group, hh_id, what)
  
  if (keep_hh) {
    
    # Keep only households for given subset
    u_hh <- getDataGroups() %>% 
      filter(school %in% schools,
             group %in% groups,
             term %in% terms) %>% 
      distinct(hh_id) %>% 
      pull(hh_id)
    
    subvert %>% 
      filter(hh_id %in% u_hh[!is.na(u_hh)] |
               what == "group")
  } else {
    subvert %>% filter(what != "hh")
  }
}

getHH <- function(sid) {
  households <- getDataHouseholds()
  hh <- with(households, hh_id[sugar_id == sid])
  if (length(hh) == 0) {
    return(NA_character_)
  } else {
    return(hh)
  }
}

#' Title
#'
#' @param sid 
#' @param term 
#' @param ignore_groups 
#'
#' @return
#' @export
#'
#' @examples
getSchoolSID <- function(sid, 
                         term, 
                         ignore_groups = FALSE) {
  
  households <- getDataHouseholds()
  groups <- getDataGroups()
  
  if (sid %in% groups$sugar_id & !ignore_groups) {
    first(groups$school[groups$sugar_id == sid & groups$term == term])
  } else {
    hh <- with(households, hh_id[sugar_id == sid])
    ref_id <- households$sugar_id[
      households$hh_id == hh &
        map_lgl(households$sugar_id, ~ . %in% groups$sugar_id)
    ] %>% 
      first()
    
    if (length(ref_id) == 0) {
      return(NA_character_)
    } else {
      return(first(groups$school[groups$sugar_id == ref_id]))
    }
  }
}

getSIDFromHUGCodbar <- function(hug_codbar, 
                                verbose = F) {
  households <- getDataHouseholds()
  
  res <- households$sugar_id[households$hug_codbar == hug_codbar]
  if (length(res) == 0) {
    if (verbose) {
      cat("-- Could not find sugar_id for hug_codbar: ", hug_codbar, "\n")
    }
    res <- NA_character_
  }
  
  res
}

getHHId <- function(sid) {
  households <- getDataHouseholds()
  
  if (is.na(sid)) {
    return(NA_character_)
  }
  
  res <- households$hh_id[households$sugar_id == sid]
  
  if (length(res) == 0) {
    return(NA_character_)
  }
  
  res
}

getGroupSID <- function(sid,
                        term,
                        match_household = T) {
  households <- getDataHouseholds()
  groups <- getDataGroups()
  
  if (is.na(sid)) {
    return(NA_character_)
  }
  
  if (sid %in% groups$sugar_id) {
    first(as.character(groups$group[groups$sugar_id == sid & groups$term == term]))
  } else if (match_household) {
    hh <- with(households, hh_id[sugar_id == sid])
    ref_id <- households$sugar_id[
      households$hh_id == hh &
        map_lgl(households$sugar_id, ~ . %in% groups$sugar_id)
    ] %>% 
      first()
    
    if (length(ref_id) == 0) {
      return(NA_character_)
    } else {
      return(as.character(first(groups$group[groups$sugar_id == ref_id])))
    }
  } else {
    return(NA_character_)
  }
}

getRefParent <- function(sid) {
  households <- getDataHouseholds()
  # Get houshold id 
  hh <- with(households, hh_id[sugar_id == sid])
  # Get all parent ids
  parents <- filter(households, hh_id == hh) 
  # Filter for ref parent if present
  ref_id <- with(parents, sugar_id[referent_famille_c == "Oui"])
  ref_id <- ref_id[!is.na(ref_id)]
  
  if (length(ref_id) == 0) {
    return(NA_character_)
  } else {
    return(ref_id)
  }
}

getAllParent <- function(sid) {
  households <- getDataHouseholds()
  # Get houshold id 
  hh <- with(households, hh_id[sugar_id == sid])
  # Get all parent ids
  parents <- filter(households, hh_id == hh) 
  
  setdiff(parents$sugar_id, sid)
}

getAllHousehold <- function(sid) {
  households <- getDataHouseholds()
  # Get houshold id 
  hh <- with(households, hh_id[sugar_id == sid])
  # Get all parent ids
  parents <- filter(households, hh_id == hh) 
  
  parents$sugar_id
}

getGroupForParent <- function(sid, school_net_people, households) {
  # Get household
  hh <- getHH(sid)
  # Get children
  chldrn <- filter(households, hh_id == hh) %>% 
    distinct(sugar_id) %>% 
    inner_join(school_net_people) %>% 
    mutate(sugar_id = sid,
           name = str_c(sugar_id, term, sep = "-")) 
  
  chldrn
}

makeNetEdges <- function(school_groups,
                         households,
                         schools = getAllSchools(),
                         groups = getAllGroups(),
                         terms = getAllTerms(),
                         keep_hh = TRUE) {
  
  subdat <- school_groups %>% 
    filter(school %in% schools,
           group %in% groups,
           term %in% terms) %>% 
    mutate(uid =  str_c(sugar_id, term, sep = "-"))
  
  school_edges <- subdat %>% 
    select(uid, sugar_id, group) %>% 
    rename(unit = group)
  
  if (keep_hh) {
    hh_edges <- filter(households, hh_id %in% subdat$hh_id) %>%
      # Expand to all terms
      mutate(pin = "x") %>% 
      inner_join(tibble(pin = "x", term = terms)) %>% 
      mutate(uid =  str_c(sugar_id, term, sep = "-")) %>% 
      select(uid, sugar_id, hh_id) %>% 
      rename(unit = hh_id)
    
    bind_rows(school_edges,
              hh_edges)
  } else {
    school_edges
  }
}

makeNetVertices <- function(school_groups,
                            households) {
  
  # First make network of kids and school staff
  school_net_people <- school_groups %>% 
    distinct(sugar_id, enseignant, term, school) %>% 
    mutate(what = case_when(enseignant == "Oui" ~ "teacher",
                            T ~ "child")) %>% 
    select(-enseignant) %>% 
    # make name if full info
    mutate(id = str_c(sugar_id, term, sep = "-")) %>% 
    distinct() %>% 
    mutate(hh_id = map_chr(sugar_id, ~ getHH(.)))
  
  school_net_groups <- school_groups %>% 
    distinct(group, term, school) %>% 
    mutate(what = "group",
           id = group) %>% 
    distinct()
  
  # Make household net
  hh_net_people <- households %>% 
    distinct(sugar_id, hh_id, age_cat) %>% 
    mutate(what = case_when(age_cat == "child" ~ "child", 
                            sugar_id %in% school_groups$sugar_id ~ "teacher",
                            T ~ "parent")) %>% 
    select(-age_cat) %>% 
    distinct()
  
  # Add only parents
  hh_net_parents <- hh_net_people %>% 
    filter(what == "parent" | !(sugar_id %in% school_net_people$sugar_id)) %>% 
    mutate(pin = "x") %>% 
    inner_join(
      tibble(pin = "x",
             term = getAllTerms())) %>%
    mutate(id = str_c(sugar_id, term, sep = "-")) %>% 
    select(-pin)
  
  # Households
  hh_net_hh <- households %>% 
    distinct(hh_id) %>% 
    mutate(what = "hh",
           id = hh_id) %>% 
    distinct()
  
  # Combine the two
  bind_rows(
    school_net_people, 
    school_net_groups,
    hh_net_parents,
    hh_net_hh
  )
}


# Genetic data utils ------------------------------------------------------


getPhyloMetadata <- function(data_path = "data/05_genomic_data/SCS_sequences_coverage50.xlsx") {
  households <- getDataHouseholds()
  school_groups <- getDataGroups()
  
  readxl::read_xlsx(data_path) %>% 
    janitor::clean_names() %>% 
    mutate(seq_id = str_remove(gisaid_name, "hCoV-19/"),
           collection_date = as.Date(collection_date),
           participant = as.character(participant)) %>% 
    rename(hug_codbar = participant) %>% 
    left_join(households %>% 
                select(hug_codbar, sugar_id)) %>% 
    rowwise() %>% 
    mutate(term = computeTerm(collection_date),
           group = getGroupSID(sugar_id, term),
           school = getSchoolSID(sugar_id, term),
           who = ifelse(sugar_id %in% school_groups$sugar_id, "school", "hh_member")) %>% 
    ungroup() %>% 
    mutate(collection_week = str_c(lubridate::year(collection_date), lubridate::epiweek(collection_date), sep = "-")) %>% 
    arrange(collection_week, group, collection_date) %>%  
    mutate(seq_num = row_number())
}


getTreeData <- function(term) {
  
  if (term == "20/21") {
    variants <- "alpha"
  } else {
    variants <- "delta_omicron"
  }
  
  phylo <- ape::read.nexus(str_glue("data/05_genomic_data/nextstrain__tree_{variants}.nexus"))
  phylo_context <- ape::read.nexus(str_glue("data/05_genomic_data/nextstrain__tree_{variants}_context.nexus"))
  
  time_phylo <- ape::read.nexus(str_glue("data/05_genomic_data/nextstrain__timetree_{variants}.nexus"))
  time_phylo_context <- ape::read.nexus(str_glue("data/05_genomic_data/nextstrain__timetree_{variants}_context.nexus"))
  
  
  return(
    list(
      phylo = phylo,
      phylo_context = phylo_context,
      time_phylo = time_phylo,
      time_phylo_context = time_phylo_context
    )
  )
}

getSeqBaseLength <- function() {
  28e3
}

toLongDistMat <- function(dist_mat, 
                          phylo_metadata, 
                          target_set = NULL) {
  dist_mat %>% 
    as.data.frame() %>% 
    mutate(from = rownames(.)) %>% 
    {
      if (!is.null(target_set)) {
        filter(., from %in% target_set) 
      } else {
        .
      }
    } %>% 
    pivot_longer(cols = contains("Swi"),
                 names_to = "to",
                 values_to = "dist")  %>% 
    # Extract school/group info
    mutate(
      across(from:to, 
             ~ map_chr(.x, ~ phylo_metadata$sugar_id[phylo_metadata$seq_id == .]),
             .names = "{.col}_sugar_id"),
      across(from:to, 
             ~ map_chr(.x, ~ computeTerm(phylo_metadata$collection_date[phylo_metadata$seq_id == .])),
             .names = "{.col}_term"),
      from_group = map2_chr(to_sugar_id, to_term, ~getGroupSID(.x, term = .y, match_household = F)),
      to_group = map2_chr(to_sugar_id, to_term, ~getGroupSID(.x, term = .y, match_household = F)),
      across(contains("sugar_id"), 
             list(
               hh_id = ~ map_chr(.x, ~getHHId(.))
             ),
             .names = "{.col %>% str_remove('_sugar_id')}_{.fn}"
      )
    ) %>% 
    filter(from != to)
}


getTreeLabels <- function(tree) {
  tree_labels <- c(tree$tip.label, 
                   tree$node.label)  
  tree_labels
}


getAncestors <- function(tree) {
  
  tree_labels <- getTreeLabels(tree)
  
  phylo_ancestors <- phangorn::Ancestors(tree) %>% 
    magrittr::set_names(tree_labels)
  
  phylo_ancestors_labelled <- lapply(phylo_ancestors, function(ind.vec) tree_labels[ind.vec])
  
  phylo_ancestors_labelled
}


computeDistToRoot <- function(tree, dist_mat) {
  
  tree_labels <- getTreeLabels(tree)
  phylo_ancestors <- getAncestors(tree)
  
  root_node <- which(map_dbl(phylo_ancestors, length) == 0)
  
  dist_to_root <- dist_mat[root_node, ] %>% 
    set_names(tree_labels) %>% 
    .[names(.) %in% tree$tip.label]
  
  dist_to_root
}

getPossibleAncestors <- function(tree,
                                 tree_time,
                                 phylo_metadata,
                                 target_set) {
  
  dist_mat <- ape::cophenetic.phylo(tree)
  dist_mat_time <- ape::cophenetic.phylo(tree_time)
  dist_mat_time_full <- ape::dist.nodes(tree_time)
  
  dist_long <- toLongDistMat(dist_mat, 
                             phylo_metadata, 
                             target_set = target_set)
  
  dist_time_long <- toLongDistMat(dist_mat_time, 
                                  phylo_metadata, 
                                  target_set = target_set)
  
  dist_to_root <- computeDistToRoot(tree = tree, 
                                    dist_mat = dist_mat_time_full)
  
  # Use dist to root to determine possible ancestors
  possible_ancestors <- purrr::map(
    tree$tip.label, 
    ~ names(which(dist_to_root <= dist_to_root[.]))
  ) %>% 
    set_names(tree$tip.label)
  
  possible_ancestors_df <- map_df(
    1:length(possible_ancestors), 
    function(x) {
      tibble(from = tree$tip.label[x],
             to = possible_ancestors[[tree$tip.label[x]]])
    }) %>% 
    mutate(from_dist_to_root = dist_to_root[from],
           to_dist_to_root = dist_to_root[to]) %>% 
    mutate(
      across(from:to, 
             ~ map_chr(.x, ~ phylo_metadata$sugar_id[phylo_metadata$seq_id == .]),
             .names = "{.col}_sugar_id"),
      across(from:to, 
             ~ map_chr(.x, ~ computeTerm(phylo_metadata$collection_date[phylo_metadata$seq_id == .])),
             .names = "{.col}_term")
    ) %>% 
    inner_join(dist_long) %>% 
    inner_join(dist_time_long %>% rename(time_dist = dist)) %>% 
    # Add sampling date from
    left_join(phylo_metadata %>% select(seq_id, collection_date), 
              by = c("from" = "seq_id"))  %>% 
    mutate(term = computeTerm(collection_date),
           from_school = map2_chr(from_sugar_id, term, ~ getSchoolSID(.x, .y)),
           to_school = map2_chr(to_sugar_id, term, ~ getSchoolSID(.x, .y)))
  
  possible_ancestors_df
}

#' unpackPhyloAdj
#'
#' @param mat 
#'
#' @return
#' @export
#'
#' @examples
unpackPhyloAdj <- function(mat) {
  purrr::map(1:2, function(i) {
    ind <- which(mat[,i] != 1000)
    if (length(ind) > 0) {
      mat[ind, i]
    }
  }) %>% 
    unlist()
}


#' getPhyloData
#'
#' @param phylo_metadata 
#' @param sugar_ids_mapping 
#' @param terms 
#' @param adj_list_df 
#'
#' @return
#' @export
#'
#' @examples
getPhyloData <- function(phylo_metadata,
                         sugar_ids_mapping,
                         terms,
                         adj_list_df) {
  
  seq_base_len <- getSeqBaseLength()
  
  u_seq_all <- c()
  n_seq_all <- 0
  possible_ancestors_df_all <- tibble()
  possible_ancestors_context_df_all <- tibble()
  dist_seq_all <- matrix()
  
  for (term in terms) {
    
    tree_data <- getTreeData(term)
    
    seq_in_run <- phylo_metadata %>% 
      filter(
        sugar_id %in% sugar_ids_mapping$sugar_id, 
        seq_id %in% tree_data$phylo$tip.label,
        !is.na(group)
      )
    
    if (nrow(seq_in_run) == 0) {
      next()
    }
    
    u_seq <- seq_in_run$seq_id %>% sort()
    u_seq_all <- c(u_seq_all, u_seq)
    n_seq_all <- n_seq_all + length(u_seq)
    
    # Get dist to root
    possible_ancestors_df <- getPossibleAncestors(tree = tree_data$phylo, 
                                                  tree_time = tree_data$time_phylo,
                                                  phylo_metadata = phylo_metadata,
                                                  target_set = u_seq) %>% 
      mutate(term = term)
    
    possible_ancestors_df_all <- bind_rows(
      possible_ancestors_df_all,
      possible_ancestors_df
    )
    
    # Distance to community sequences
    possible_ancestors_context_df <- getPossibleAncestorsContext(tree = tree_data$phylo_context, 
                                                                 tree_time = tree_data$time_phylo_context,
                                                                 phylo_metadata = phylo_metadata,
                                                                 target_set = u_seq) %>% 
      mutate(term = term)
    
    possible_ancestors_context_df_all <- bind_rows(
      possible_ancestors_context_df_all,
      possible_ancestors_context_df
    )
    
    # Matrix of sequence distances
    dist_mat <- ape::cophenetic.phylo(tree_data$phylo)
    dist_seq <- dist_mat[u_seq, u_seq]
    dist_seq_all <- magic::adiag(dist_seq_all, dist_seq, pad = 999)
    
  }
  
  # Case where no genetic data is available for data subset
  if (nrow(possible_ancestors_df_all) == 0) {
    return(
      list(
        do_gen = T,
        adj_list_ancestors = c(), 
        adj_list_infectors = c(), 
        adj_type_infectors = c(), 
        dist_closest_comm_seq = c(), 
        closest_comm_seq_len = c(), 
        seq_len_mat = c(), 
        dist_seq = c(), 
        start_list_seq = c(), 
        end_list_seq = c(), 
        n_max_ancestors = c(), 
        temp_diff_seq = c(),
        closest_comm_seq_temp_diff = c(),
        u_seq_all = c()
      )
    )
  }
  
  
  dist_seq_all <- dist_seq_all[-1, -1]
  
  # Ancestors within school groups
  possible_ancestors_in_school <- possible_ancestors_df_all %>% 
    filter(!is.na(from_sugar_id),  !is.na(to_sugar_id), 
           !is.na(from_group), !is.na(to_group)) %>% 
    filter(from_sugar_id %in% sugar_ids_mapping$sugar_id,
           to_sugar_id %in% sugar_ids_mapping$sugar_id) %>% 
    # !! Keep only within-school ancestry
    filter(from_school == to_school)  %>%
    # Here from and to are in terms of the sequence number
    mutate(across(from:to, 
                  function(x) map_dbl(x, ~ which(u_seq_all == .)),
                  .names = "{.col}_num")) %>% 
    arrange(from_num, to_num) %>% 
    mutate(row_id = row_number())
  
  
  # Build adjacency
  adj_list_ancestors_df <- possible_ancestors_in_school %>%
    select(from, to, year = term, row_id) %>% 
    mutate(across(from:to, 
                  function(x) map_dbl(x, ~ which(u_seq_all == .)))) %>% 
    arrange(from, to) 
  
  adj_list_infectors_df <- possible_ancestors_in_school %>% 
    select(from_sugar_id, to_sugar_id, term, row_id) %>% 
    mutate(across(from_sugar_id:to_sugar_id, 
                  function(x) map_dbl(x, ~ which(sugar_ids_mapping$sugar_id == .)),
                  .names = "{.col %>% str_remove('_sugar_id')}")) %>% 
    mutate(from_num = from,
           from = sugar_ids_mapping$id[from_num],
           to_num = to,
           to = sugar_ids_mapping$id[to_num]) %>% 
    select(from, to, from_num, to_num, term, row_id) %>% 
    inner_join(adj_list_df, by = c("from" = "from_chr", "to" = "to_chr", "term" = "year")) %>% 
    arrange(from, to, term) 
  
  adj_list_ancestors_df <- adj_list_ancestors_df %>% filter(row_id %in% adj_list_infectors_df$row_id)
  possible_ancestors_in_school <- possible_ancestors_in_school %>% filter(row_id %in% adj_list_infectors_df$row_id)
  
  seq_len_mat <- dist_seq_all * 0 + seq_base_len
  temp_diff_seq <- makeTimeDiffMatrix(possible_ancestors_in_school,
                                      u_seq = u_seq_all)
  
  start_list_seq <- getStartEndAdj(adj_list_ancestors_df, what = "start") %>% 
    unpackPhyloAdj() %>% 
    {. - 1}
  start_list_seq[is.na(start_list_seq)] <- 999
  end_list_seq <- getStartEndAdj(adj_list_ancestors_df, what = "end") %>% 
    unpackPhyloAdj() %>% 
    {. - 1}
  end_list_seq[is.na(end_list_seq)] <- 999
  
  # Max number of ancestors
  n_max_ancestors <- adj_list_ancestors_df %>% count(from) %>% pull(n)
  
  adj_list_ancestors <- adj_list_ancestors_df$to - 1
  # Here we have one sequence per sample and assume both are the same
  adj_list_infectors <- adj_list_infectors_df$to_num - 1
  adj_type_infectors <- map_dbl(1:nrow(adj_list_infectors_df), function(x) {
    adj_list_df %>% 
      filter(from_chr == adj_list_infectors_df$from[x],
             to_chr == adj_list_infectors_df$to[x],
             year == adj_list_infectors_df$term[x]) %>%
      pull(type)
  })
  
  # The earliest sample cannot have any ancestors within the sampled sequences of interest
  # so we loose one
  
  possible_ancestors_community <- possible_ancestors_context_df_all %>% 
    filter(from %in% u_seq_all, !(to %in% u_seq_all), 
           # !! Filtering for spatial proximity
           str_detect(to, "Switzerland|France")) %>% 
    group_by(from) %>% 
    slice_min(dist, with_ties = F) %>% 
    ungroup() %>% 
    mutate(from_num = map_dbl(from, ~ which(u_seq_all == .))) %>% 
    arrange(from_num, to)
  
  dist_closest_comm_seq <- possible_ancestors_community$dist
  closest_comm_seq_len <- rep(seq_base_len, n_seq_all)
  closest_comm_seq_temp_diff <- map_dbl(
    1:nrow(possible_ancestors_community),
    function(x) {
      ceiling(365*(possible_ancestors_community$from_dist_to_root[x]  - 
                     possible_ancestors_community$to_dist_to_root[x]))
    })
  
  
  gen <- list(
    do_gen = T,
    adj_list_ancestors = adj_list_ancestors, 
    adj_list_infectors = adj_list_infectors, 
    adj_type_infectors = adj_type_infectors, 
    dist_closest_comm_seq = dist_closest_comm_seq, 
    closest_comm_seq_len = closest_comm_seq_len, 
    seq_len_mat = seq_len_mat, 
    dist_seq = dist_seq, 
    start_list_seq = start_list_seq, 
    end_list_seq = end_list_seq, 
    n_max_ancestors = n_max_ancestors, 
    temp_diff_seq = temp_diff_seq,
    closest_comm_seq_temp_diff = closest_comm_seq_temp_diff,
    u_seq_all = u_seq_all
  )
  
  gen
}

makeWMat <- function(max_kappa,
                     max_delay = 50,
                     w) {
  
  w_mat <- matrix(0, max_kappa, max_delay)
  w <- w/sum(w)
  input <- c(1, rep(0, max_delay-1))
  
  for (i in 1:max_kappa) {
    tmp <- input
    for (k in 1:i) {
      tmp <- convolve(tmp, rev(w), type = "open")
    }
    w_mat[i, ] <- pmax(1e-12, tmp[1:max_delay])
  }
  
  w_mat <- apply(w_mat, 2, function(x) x/sum(x)) %>% 
    matrix(nrow = max_kappa)
  
  w_mat
}


getPossibleAncestorsContext <- function(tree,
                                        tree_time,
                                        phylo_metadata,
                                        target_set) {
  
  dist_mat <- ape::cophenetic.phylo(tree)
  dist_mat_time <- ape::cophenetic.phylo(tree_time)
  dist_mat_time_full <- ape::dist.nodes(tree_time)
  
  dist_to_root <- computeDistToRoot(tree = tree, 
                                    dist_mat = dist_mat_time_full)
  
  # Use dist to root to determine possible ancestors
  possible_ancestors <- purrr::map(
    target_set, 
    ~ names(which(dist_to_root <= dist_to_root[.]))
  ) %>% 
    set_names(target_set)
  
  possible_ancestors_df <- map_df(
    1:length(possible_ancestors), 
    function(x) {
      tibble(from = target_set[x],
             to = possible_ancestors[[target_set[x]]]) %>% 
        mutate(dist = dist_mat[from[1], to],
               from_dist_to_root = dist_to_root[from],
               to_dist_to_root = dist_to_root[to])
    }) 
  
  possible_ancestors_df
}




makeTimeDiffMatrix <- function(possible_ancestors_in_school,
                               u_seq) {
  
  temp_diff_seq <- map_df(
    1:nrow(possible_ancestors_in_school),
    function(x) {
      f <- possible_ancestors_in_school$from[x]
      t <- possible_ancestors_in_school$to[x]
      tibble(
        from  = f,
        to = t,
        dt = round((possible_ancestors_in_school$from_dist_to_root[x] - 
                      possible_ancestors_in_school$to_dist_to_root[x]) * 365)
      )
    }) %>% 
    complete(from = u_seq,
             to = u_seq) %>% 
    pivot_wider(names_from = "to",
                values_from = "dt") %>% 
    {
      df <- .
      mat <- as.matrix(df %>% select(-from))
      rownames(mat) <- df$from
      colnames(mat) <- colnames(df)[-1]
      mat[u_seq, u_seq]
    }
  
  temp_diff_seq[is.na(temp_diff_seq)] <- 999
  temp_diff_seq
}

dateStartAlpha <- function() {
  as.Date("2021-01-25")
}

dateStartDelta <- function() {
  as.Date("2021-06-14")
}

dateStartOmicron <- function() {
  as.Date("2021-12-17")
}


getVariant <- function(date) {
  case_when(date < dateStartAlpha() ~ "wildtype",
            date < dateStartDelta() ~ "alpha",
            date < dateStartOmicron() ~ "delta",
            date >= dateStartOmicron() ~ "omicron",
            is.na(date) ~ NA_character_,
            T ~ NA_character_)
}


getVariantDates <- function(var) {
  if (is.na(var)) {
    stop("Variant is NA")
  } else if (var == "wildtype") {
    TL <- as.Date("2021-01-01")
    TR <- studyStartDate()
  } else if (var == "alpha") {
    TL <- dateStartAlpha()
    TR <- dateStartDelta()
  } else if (var == "delta") {
    TL <- dateStartDelta()
    TR <- dateStartOmicron()
  } else if (var == "omicron") {
    TL <- dateStartOmicron()
    TR <- studyEndDate()
  }
  return(
    tibble(variant = var,
           TL = TL,
           TR = TR)
  )
}

#' parseDate
#' Parse the date of a sequence from json info
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
parseDate <- function(x) {
  
  tibble(node_name = x$name,
         time = x$node_attrs$num_date$value)
}

#' Step children
#' Step through tree to extract metadata from json
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
stepChildren <- function(x) {
  
  if ("children" %in% names(x)) {
    bind_rows(parseDate(x), 
              map_df(x$children, ~stepChildren(.)))
  } else {
    parseDate(x)
  }
}

#' getTreeNodeDates
#' Get the tree node dates from json metadata
#' 
#' @param variants 
#' @param redo 
#'
#' @return
#' @export
#'
#' @examples
getTreeNodeDates <- function(variants,
                             redo = FALSE) {
  
  node_times_file <- str_glue("generated_data/node_times_{variants}.rds")
  
  if (!file.exists(node_times_file) | redo) {
    
    if (variants == "alpha") {
      tree_metadata <- rjson::fromJSON(file = "data/05_genomic_data/ncov_geneva-alpha-2021-06-10.json")
    } else if (variants == "delta_omicron") {
      tree_metadata <- rjson::fromJSON(file = "data/05_genomic_data/ncov_geneva_focus-2022-02-26.json")
    } else {
      stop("Unknown variants")
    }
    
    node_times <- stepChildren(tree_metadata$tree) %>% 
      mutate(date = timeToDate(time - 2020, ref_date = "2020-01-01"))
    
    saveRDS(node_times, file = node_times_file)
    
  } else {
    node_times <- readRDS(node_times_file)
  }
  
  node_times
}


#' viewClade2
#' based on code in ggtree
#' @param tree_view 
#' @param node 
#' @param xmax_adjust 
#'
#' @return
#' @export
#'
#' @examples
viewClade2 <- function (tree_view = NULL, node, xmax_adjust = 0) 
{
  # tree_view %<>% ggtree:::get_tree_view
  cpos <- get_clade_position(tree_view, node = node)
  xmax <- aplot::xrange(tree_view)[2]
  attr(tree_view, "viewClade") <- TRUE
  attr(tree_view, "viewClade_node") <- node
  tree_view + coord_cartesian(xlim = c(cpos$xmin, cpos$xmax + xmax_adjust), 
                              ylim = c(cpos$ymin, cpos$ymax), expand = FALSE)
}



getOutbreakTimesToMRCA <- function(ob_id,
                                   phylo_meta) {
  
  submeta <- phylo_meta %>% 
    filter(outbreak_id == ob_id) %>% 
    select(-school, -group) %>% 
    # !! get sequences in this school
    inner_join(school_groups)
  
  if (nrow(submeta) < 2) {
    cat("Not enough sequences found for outbreak", ob_id, ", returning NULL\n")
    return(NULL)
  }
  
  this_variant <- submeta$variant[1]
  
  if (this_variant == "alpha") {
    variants <- this_variant
  } else {
    variants <- "delta_omicron"
  }
  
  node_times <- getTreeNodeDates(variants = variants)
  ob_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__timetree_{variants}.nexus"))
  context_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__timetree_{variants}_context.nexus"))
  
  # Split delta/omicron
  if (this_variant == "alpha") {
    root_node <- "NODE_0002967"
  } else if (this_variant == "delta") {
    root_node <- "NODE_0001867"
  } else {
    root_node <- "NODE_0001726"
  }
  
  context_tree <- ape::extract.clade(context_tree, node = root_node)
  
  node_times <- node_times %>% filter(node_name %in% c(context_tree$node.label, root_node))
  this_tip_labels <-  submeta$seq_id[submeta$seq_id %in% context_tree$tip.label]
  
  # Get the MRCA of outbreak sequence
  mrca_node <- MRCA(context_tree, this_tip_labels) 
  mrca_node_label <- c(context_tree$tip.label, context_tree$node.label)[mrca_node]
  mrca_node_date <- node_times$date[node_times$node_name == mrca_node_label]
  
  # Get the distances
  dist_mat_time_full <- ape::dist.nodes(context_tree)
  dist_to_mrca <- dist_mat_time_full[mrca_node, ]
  names(dist_to_mrca) <- getTreeLabels(context_tree)
  this_dist_to_mrca <- dist_to_mrca[this_tip_labels]
  
  return(this_dist_to_mrca)
}

# Plotting utils ----------------------------------------------------------

plotParticipantData <- function(sid) {
  long_dat <- readRDS("generated_data/combined_time_data.rds") %>% 
    filter(sugar_id %in% sid)
  
  outbreak_dat <- readRDS("generated_data/oubtreak_data.rds") %>% 
    filter(school == getSchool(sid)) %>% 
    pivot_longer(cols = c("visite_1", "visite_2", "visite_3"),
                 values_to = "date",
                 names_to = "visite") %>% 
    mutate(visite_full = str_c("outbreak_", outbreak_id, ":", visite))
  
  if (nrow(long_data) == 0) 
    stop("No data found to plot")
  
  ggplot(long_data, aes(x = event_date, y = event_type)) +
    geom_vline(data = outbreak_dat, 
               aes(xintercept = date,
                   color = outbreak_id,
                   lty = visite), alpha = .5) +
    geom_line(alpha = .5, size = .7) +
    geom_point(aes(pch = event), size = 2.5) +
    scale_shape_manual(values = c(15, 7, 8)) +
    facet_wrap(~ sugar_id) +
    theme_bw()  +
    scale_x_date(limits = as.Date(c("2021-03-01", "2022-07-01")),
                 date_breaks = "2 months",
                 date_labels = "%Y-%b")
  
}

plotEvenData <- function(df, 
                         event_types = getAllEventTypes(),
                         schools = getAllSchools()) {
  sdf <- df %>% 
    filter(event_type %in% event_types, 
           school %in% schools)
  
  if ("serology" %in% event_types) {
    cols <- c("blue", "purple", "red")
  } else {
    cols <- c("blue", "red")
  }
  
  p <- sdf %>% 
    ggplot(aes(x = event_date, y = sugar_id)) +
    geom_line(size = .2, alpha = .5)  +
    geom_point(data = sdf %>% filter(event != "positif"), 
               aes(pch = event_type, color = event), alpha = 1) +
    geom_point(data = sdf %>% filter(event == "positif"), 
               aes(pch = event_type, color = event), alpha = 1) +
    facet_grid(school ~., scales =  "free_y", space = "free_y") +
    # scale_color_manual(values = cols) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())
  p
}

plotHHcomp <- function(x, y, what_num = c("home_n", "home_n_adult", "home_n_child")) {
  if (x == y) {
    n_tot <- hh_comp %>% 
      ungroup() %>% 
      filter(what == what_num, !is.na(!!sym(x))) %>% 
      distinct(hh_id) %>% 
      nrow()
    
    p_hh_comp <-  hh_comp %>% 
      ungroup() %>% 
      filter(what == what_num) %>% 
      ggplot(aes(x = !!sym(x))) +
      geom_histogram() +
      coord_cartesian(xlim = c(0, 10), ylim = c(0, 120)) +
      scale_x_continuous(breaks = seq(0, 10, by = 2)) +
      theme_bw() +
      ggtitle(str_glue("{x}, n = {n_tot}")) +
      theme(plot.title = element_text(size = 10))
  } else {
    
    p_hh_comp <- hh_comp %>% 
      ungroup() %>% 
      count(what, !!sym(x), !!sym(y)) %>%
      filter(what %in% what_num, !is.na(!!sym(x)) | !is.na(!!sym(y))) %>% 
      ggplot(aes(x = !!sym(x), y = !!sym(y), size = n)) +
      geom_abline(aes(intercept = 0, slope = 1), lty = 2, col = "darkgray") +
      geom_point() +
      scale_size_continuous(breaks = seq(0, 200, by = 30), range = c(.5, 5), limits = c(1, 130)) +
      theme_bw() +
      coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
      scale_x_continuous(breaks = seq(0, 10, by = 2)) +
      scale_y_continuous(breaks = seq(0, 10, by = 2)) +
      guides(size = "none")
  }
  
  p_hh_comp
}


plotOutbreakData <- function(df) {
  school <- unique(df$school)
  term <- unique(df$term)
  id <- unique(df$outbreak_id)
  
  title <- str_glue("Outbreak {id}: {school} {term}")
  
  pch_levels <- c("pcr/antigen", "sero", "vacc dose 1", "vacc dose 2", "vacc dose 3")
  
  df <- df %>% 
    group_by(full_id, date) %>% 
    arrange(full_id, date, what) %>% 
    mutate(tid = row_number()) %>% 
    ungroup() %>% 
    mutate(full_id_2 = str_c(full_id, tid, sep = "_"))  %>% 
    mutate(what = factor(what, levels = pch_levels))
  
  
  vacc_df <- df %>% 
    select(full_id, school, group, d1_date, d2_date, d3_date) %>% 
    pivot_longer(cols = contains("date"),
                 names_to = "dose",
                 values_to = "date") %>% 
    filter(!is.na(date)) %>% 
    mutate(dose = str_replace(dose, "d", "vacc dose ") %>% str_remove("_date")) %>% 
    filter(date >= min(df$date), date <= max(df$date)) %>% 
    inner_join(df %>% 
                 filter(!is.na(date), !is.na(result)) %>% 
                 distinct(full_id, school, group)) %>% 
    mutate(dose = factor(dose, levels = pch_levels))
  
  df %>% 
    ggplot(aes(x = date, y = full_id)) +
    geom_rect(data = distinct(df, school, group, visite_1, visite_2, visite_3) %>% 
                pivot_longer(cols = contains("visite")), 
              inherit.aes = F,
              aes(xmin = value - 1, xmax = value + 1, ymin = -Inf, ymax = Inf, fill = name), 
              alpha = .4) +
    geom_line(linewidth = .2, color = "black", alpha = .7) +
    geom_point(aes(pch = what, color = factor(result)), alpha = .8, 
               position = position_dodge(.5)) +
    geom_point(inherit.aes = F,
               data = vacc_df, 
               aes(x = date, y = full_id, pch = dose)) +
    facet_grid(group ~ ., space = "free", scales = "free_y") +
    theme_bw() +
    theme(
      # axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) +
    scale_color_manual(values = c("blue", "red")) +
    scale_fill_manual(values = c("#E3CB7B", "#9EA82B", "#8F7A24")) +
    guides(color = guide_legend("Test result")) +
    ggtitle(title)
}


# Plot variant trees
plot_variant_tree <- function(this_variant, 
                              phylo_meta,
                              what = "time") {
  
  if (this_variant == "alpha") {
    variants <- this_variant
  } else {
    variants <- "delta_omicron"
  }
  
  node_times <- getTreeNodeDates(variants = variants)
  
  if (what == "phylo") {
    ob_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__tree_{variants}.nexus"))
    context_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__tree_{variants}_context.nexus"))
  } else {
    ob_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__timetree_{variants}.nexus"))
    context_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__timetree_{variants}_context.nexus"))
    
    # Split delta/omicron
    if (this_variant == "alpha") {
      root_node <- "NODE_0002967"
    } else if (this_variant == "delta") {
      root_node <- "NODE_0001867"
    } else {
      root_node <- "NODE_0001726"
    }
    context_tree <- ape::extract.clade(context_tree, node = root_node)
  }
  
  node_times <- node_times %>% filter(node_name %in% c(context_tree$node.label))
  
  # Make treeio
  tree_dat <- as.treedata(context_tree) %>% 
    full_join(phylo_meta %>% 
                filter(variant == this_variant), by = c("label" = "seq_id"))
  
  if (what == "time") {
    
    p_tree <- ggtree(tree_dat,
                     as.Date = what == "time", 
                     mrsd = max(node_times$date),
                     color = "darkgray",
                     lwd = .2) +
      theme_tree2() +
      scale_x_date(date_labels = "%Y-%m-%d",
                   date_breaks = "2 months",
                   limits = c(min(node_times$date) - 20, max(node_times$date))) 
  } else {
    # Tree
    p_tree <- ggtree(tree_dat,
                     color = "darkgray",
                     lwd = .2)+
      theme_tree2()
  }
  
  p_tree  +
    geom_point2(aes(subset = (isTip & !is.na(school)), 
                    fill = school),
                color = "white",
                pch = 21,
                size = 3.5) +
    scale_color_manual(values = getSchoolColors(), drop = FALSE) +
    scale_fill_manual(values = getSchoolColors(), drop = FALSE) +
    ggtitle(str_glue("{str_to_title(this_variant)}"))
}


plot_outbreak_tree <- function(ob_id, 
                               phylo_meta,
                               what = "phylo") {
  submeta <- phylo_meta %>% 
    filter(outbreak_id == ob_id) %>% 
    select(-school, -group) %>% 
    inner_join(school_groups)
  
  this_variant <- submeta$variant[1]
  
  if (this_variant == "alpha") {
    variants <- this_variant
  } else {
    variants <- "delta_omicron"
  }
  
  node_times <- getTreeNodeDates(variants = variants)
  
  if (what == "phylo") {
    ob_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__tree_{variants}.nexus"))
    context_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__tree_{variants}_context.nexus"))
  } else {
    ob_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__timetree_{variants}.nexus"))
    context_tree <- read.nexus(str_glue("data/05_genomic_data/nextstrain__timetree_{variants}_context.nexus"))
  }
  
  # Split delta/omicron
  if (this_variant == "alpha") {
    root_node <- "NODE_0002967"
  } else if (this_variant == "delta") {
    root_node <- "NODE_0001867"
  } else {
    root_node <- "NODE_0001726"
  }
  
  # Extract the tree for this root 
  context_tree <- ape::extract.clade(context_tree, node = root_node)
  node_times <- node_times %>% filter(node_name %in% c(context_tree$node.label, root_node))
  
  # Get outbreak tips that are in trees
  ob_tips <- submeta$seq_id[(submeta$seq_id %in% ob_tree$tip.label) & (submeta$seq_id %in% context_tree$tip.label)]
  
  # Get the MRCA of outbreak sequence
  mrca_node <- MRCA(context_tree, ob_tips) 
  mrca_node_label <- c(context_tree$tip.label, context_tree$node.label)[mrca_node]
  mrca_node_date <- node_times$date[node_times$node_name == mrca_node_label]
  
  context_tree <- ape::extract.clade(context_tree, node = mrca_node_label)
  
  # Make treeio
  #  data
  tree_dat <- as.treedata(context_tree) %>% 
    full_join(phylo_meta %>% 
                mutate(this_outbreak = outbreak_id == ob_id,
                       this_school = school == submeta$school[1]), 
              by = c("label" = "seq_id"))
  
  if (what == "time") {
    # Tree
    p_tree <- ggtree(tree_dat,
                     as.Date = what == "time", 
                     mrsd = max(node_times$date),
                     color = "darkgray") +
      theme_tree2() +
      scale_x_date(date_labels = "%Y-%m-%d",
                   date_breaks = "1 months",
                   limits = c(min(node_times$date) - 20, max(node_times$date)))
    
  } else {
    # Tree
    p_tree <- ggtree(tree_dat,
                     color = "darkgray") +
      theme_tree2() 
  }
  
  # Get the dates of all nodes in clade
  # p_tree_zoom <- viewClade2(p_tree, node = mrca_node)
  p_tree_zoom <- p_tree
  
  # these_tips <- which(p_tree$data$label %in% submeta$seq_id[submeta$seq_id %in% ob_tree$tip.label])
  these_tips <- which(p_tree$data$label %in% phylo_meta$seq_id[phylo_meta$collection_date >= (min(submeta$collection_date) - 7)])
  all_tips <- which(p_tree$data$label %in% context_tree$tip.label)
  
  # Get all ancestor nodes
  these_ancestor_nodes <- purrr::map(these_tips, ~phangorn::Ancestors(context_tree, .))
  all_ancestor_nodes <- purrr::map(all_tips, ~phangorn::Ancestors(context_tree, .))
  
  # Get all tips that are not linked to the data
  diff_nodes <- setdiff(unlist(all_ancestor_nodes), unlist(these_ancestor_nodes))
  
  # Get all nodes that are not ancestors of others
  keep_nodes <- map_lgl(diff_nodes, function(x) {
    # Nodes that are ancestors of this node and that are among the diff nodes
    a <- intersect(phangorn::Ancestors(context_tree, x), diff_nodes)
    if (length(a) == 0) {
      TRUE
    } else {
      FALSE
    }
  })
  
  for (i in diff_nodes[keep_nodes]) {
    p_tree_zoom <- p_tree_zoom %>%
      collapse(node = i, mode = "max", fill = "#E0E0E0", color = "darkgray", lwd = .2, scale = .1)
  }
  
  p_tree_zoom +
    geom_point2(aes(subset = (isTip & !is.na(school) & this_school), 
                    fill = school),
                color = "white",
                pch = 21,
                size = 4.5) +
    geom_point2(aes(subset = (isTip & !is.na(school) & !this_school), 
                    color = school, fill = school),
                # fill = "gray",
                pch = 21,
                size = 2) +
    scale_color_manual(values = getSchoolColors(), drop = FALSE) +
    scale_fill_manual(values = getSchoolColors(), drop = FALSE) +
    guides(color = "none", fill = "none") +
    ggtitle(str_glue("Outbreak {ob_id}")) +
    coord_cartesian(xlim = c(0, 20)) +
    labs(x = "# mutations from MRCA") +
    theme(title = element_text(size = 8))
}


# Seroconversion ----------------------------------------------------------

joinSeroData <- function(df,
                         sero_data,
                         time_lim_lo,
                         time_lim_hi) {
  df %>% 
    filter(!is.na(date_min), !is.na(date_max),
           abs(date_min_ref_diff) < time_lim_lo,
           abs(date_max_ref_diff) < time_lim_hi) %>% 
    left_join(sero_data %>% 
                select(sugar_id, 
                       date_min = date_rdv, 
                       interp_date_min = nia_epfl_interp_simple)) %>% 
    left_join(sero_data %>% 
                select(sugar_id, 
                       date_max = date_rdv, 
                       interp_date_max = nia_epfl_interp_simple))  
}

countSerochanges <- function(seroconv,
                             sero_data,
                             time_lim_hi = 180,
                             time_lim_lo = 180,
                             by_pcr = FALSE) {
  
  x <- seroconv %>% 
    group_by(sugar_id, outbreak_id, group, school, age_cat) 
  
  if (by_pcr) {
    x <- seroconv %>% 
      group_by(sugar_id, outbreak_id, group, school, age_cat, any_pos_test) 
  }
  
  x %>% 
    filter(!prior_pcrpos, !prior_vacc) %>% 
    joinSeroData(sero_data = sero_data,
                 time_lim_hi = time_lim_hi,
                 time_lim_lo = time_lim_lo) %>% 
    count(interp_date_min, interp_date_max) %>% 
    inner_join(x %>% 
                 summarise(date_min = median(date_min, na.rm = T),
                           date_max = median(date_max, na.rm = T))) %>% 
    # !! Drop duplicate data for educators in multiple groups
    group_by(sugar_id, outbreak_id) %>% 
    arrange(outbreak_id, sugar_id, desc(date_min)) %>% 
    slice(1) %>% 
    ungroup()
  
}

summarizeSerochanges <- function(seroconv_counts) {
  seroconv_counts %>% 
    summarise(n_seroconv = sum(n[interp_date_min == "non-positive" &
                                   interp_date_max == "positive"]),
              n_neg_baseline = sum(n[interp_date_min == "non-positive"]),
              n_serorev = sum(n[interp_date_min == "positive" &
                                  interp_date_max == "non-positive"]),
              n_stable = sum(n[interp_date_min == interp_date_max]),
              n_tot = sum(n),
              date_min = date_min[1],
              date_max = date_max[1]) %>% 
    rowwise() %>% 
    mutate(frac_conv = n_seroconv/n_neg_baseline,
           lo = Hmisc::binconf(n_seroconv, n_neg_baseline)[2],
           hi = Hmisc::binconf(n_seroconv, n_neg_baseline)[3]) %>% 
    ungroup()
}


# Natural history parameters ----------------------------------------------

# These functions are all based on Hart 2022a/b:
# Hart, William S., et al. "Generation time of the alpha and delta SARS-CoV-2 variants: an epidemiological analysis." The Lancet Infectious Diseases 22.5 (2022): 603-610.
# Hart, William S., et al. "Inference of the SARS-CoV-2 generation time using UK household data." ELife 11 (2022): e70767.

computeC <- function(k_inc, k_P, gamma, mu, alpha_P_I) {
  k_inc * gamma * mu/(alpha_P_I * k_P * mu + k_inc * gamma)
}

computeEystar <- function(k_P, k_inc, gamma, alpha_P_I, k_I, mu) {
  E_P <- k_P * 1/(k_inc * gamma)
  E_P2 <- mean(rgamma(1e4, shape = k_P, scale = 1/(gamma*k_inc))^2)
  E_I <- mu
  E_I2 <- mean(rgamma(1e4, shape = k_I, scale = mu/k_I)^2)
  C <- computeC(k_inc, k_P, gamma, mu, alpha_P_I)
  E_ystar <- C/2 * (alpha_P_I*E_P2 + 2 * E_P * E_I + E_I2)
  E_ystar
}

computeEgen <- function(k_P, k_inc, gamma, alpha_P_I, k_I, mu) {
  E_ystar <- computeEystar(k_P, k_inc, gamma, alpha_P_I, k_I, mu)
  E_E <- (k_inc - k_P)/(k_inc*gamma)
  E_ystar + E_E
}

computeVystar <- function(k_P, k_inc, gamma, alpha_P_I, k_I, mu) {
  E_P <- k_P * 1/(k_inc * gamma)
  E_P2 <- mean(rgamma(1e4, shape = k_P, scale = 1/(gamma*k_inc))^2)
  E_P3 <- mean(rgamma(1e4, shape = k_P, scale = 1/(gamma*k_inc))^3)
  E_I <- mu
  E_I2 <- mean(rgamma(1e4, shape = k_I, scale = mu/k_I)^2)
  E_I3 <- mean(rgamma(1e4, shape = k_I, scale = mu/k_I)^3)
  C <- computeC(k_inc, k_P, gamma, mu, alpha_P_I)
  e_g <- computeEystar(k_P = k_P, 
                       k_inc = k_inc, 
                       gamma = gamma, 
                       alpha_P_I = alpha_P_I, 
                       k_I = k_I,
                       mu = mu) 
  V_ystar <- C/3 * (alpha_P_I*E_P3 + 3 * E_P2 * E_I +
                      + 3 * E_P * E_I2 +  E_I3) - e_g^2
  V_ystar
}

computeVgen <- function(k_P, k_inc, gamma, alpha_P_I, k_I, mu)  {
  V_ystar <- computeVystar(k_P, k_inc, gamma, alpha_P_I, k_I, mu)
  V_E <- (k_inc - k_P)/(k_inc*gamma)^2
  V_ystar + V_E
}


# We here assume that all of the infectious period is a single period
objfunErlangs <- function(param,
                          k_P,
                          k_I,
                          k_inc,
                          alpha_P_I,
                          gamma,
                          e_gen,
                          v_gen) {
  
  
  # Assume same
  mu_I <- param[1]
  
  # Compute the generation time expectation
  e_g <- computeEgen(k_P = k_P, 
                     k_inc = k_inc, 
                     gamma = gamma, 
                     alpha_P_I = alpha_P_I, 
                     k_I = k_I,
                     mu = mu_I) 
  
  v_g <- computeVgen(k_P = k_P, 
                     k_inc = k_inc, 
                     gamma = gamma, 
                     alpha_P_I = alpha_P_I, 
                     k_I = k_I,
                     mu = mu_I) 
  
  if (is.nan(v_g) | v_g < 0) {
    error <- 1e2
  } else {
    error <- ((e_gen - e_g)^2 + (sqrt(v_gen) - sqrt(v_g))^2)
  }
  
  error
}


optimizeErlangs <- function(k_inc,
                            gamma,
                            alpha_P_I,
                            e_gen,
                            v_gen) {
  map_df(
    1:(k_inc-1), 
    function(x) 
    {
      map_df(
        1:3, 
        function(y) 
        {
          res <- optimize(
            f = objfunErlangs,
            interval = c(0, 6),
            k_I = y,
            k_P = x,
            alpha_P_I = alpha_P_I,
            k_inc = k_inc,
            gamma = gamma,
            e_gen = e_gen,
            v_gen = v_gen
          )
          
          e_g <- computeEgen(k_P = x, 
                             k_inc = k_inc, 
                             gamma = gamma, 
                             alpha_P_I = alpha_P_I, 
                             k_I = y,
                             mu = res$minimum) 
          
          v_g <- computeVgen(k_P = x, 
                             k_inc = k_inc, 
                             gamma = gamma, 
                             alpha_P_I = alpha_P_I, 
                             k_I = y,
                             mu = res$minimum) 
          
          tibble(
            err = res$objective,
            k_P = x,
            k_I = y,
            mu = res$minimum,
            e_g = e_g,
            v_g = v_g
          )
        })
    })
}

getNatHistParams <- function() {
  
  # Load estimates
  nathist_params <- readRDS("generated_data/variant_nathis_erlang_param.rds")
  
  # Add dummy estimate for ancestral
  nathist_params <- nathist_params %>% 
    bind_rows(nathist_params %>% slice(1) %>% mutate(variant = "ancestral")) %>% 
    mutate(variant = factor(variant, levels = c("ancestral", "alpha", "delta", "omicron"))) %>% 
    arrange(variant)
  
  # Susceptibility once infected
  # Matrix based on:
  # Stein et al. (2023), “Past SARS-CoV-2 infection protection against re-infection: 
  # a systematic review and meta-analysis”, The Lancet, Vol. 401 No. 10379, 
  # pp. 833–842, doi: 10.1016/S0140-6736(22)02465-5.
  # 
  # Note that the protection from Omicron infection to prior-to-omicron 
  # variants is set arbitrarily as it is not used
  
  nathist_params$delta <- 1 - matrix(
    c(
      rep(c(.84, .9, .82, .45), 3),
      c(.6, .6, .6, .8)
    ),
    nrow = 4,
    byrow = TRUE)
  
  colnames(nathist_params$delta) <- c("ancestral", "alpha", "delta", "omicron")
  rownames(nathist_params$delta) <- colnames(nathist_params$delta)
  
  nat_hist <- list(
    scov2_gamma = nathist_params$scale_inc * 365,    # rate of E -> P
    scov2_eta = nathist_params$scale_inc * 365,      # rate of Pn -> Pn+1
    scov2_phi = nathist_params$scale_I * 365,        # rate of I -> Sprime 
    scov2_delta = nathist_params$delta,          # fraction of susceptibility once infected
    alpha = nathist_params$alpha_P_I       # fraction of susceptibility once infected
  )
  
  nat_hist
}


# POMP --------------------------------------------------------------------

#' aggregateSim
#'
#' @param sim 
#' @param seis_unit_statenames 
#' @param by_school 
#' @param by_group 
#' @param simplify_states 
#'
#' @return
#' @export
#'
#' @examples
aggregateSim <- function(sim, 
                         seis_unit_statenames,
                         by_school = FALSE,
                         by_group = FALSE,
                         simplify_states = FALSE) {
  
  if (simplify_states) {
    seis_unit_statenames <- c("S", "Sprime", "E", "Itot", "C")
    
    sim <- sim %>% 
      # Make the totals of currently infectious (Itot)
      mutate(Itot = Pone + Ptwo + Pthree + I) %>% 
      select(time, date, .id, unitname, 
             all_of(c(seis_unit_statenames, "pcr_pos", "sero")), 
             any_of(c("group", "school")))
  }
  
  # Pivot longer
  long_sim <- sim %>%
    pivot_longer(cols = any_of(c(seis_unit_statenames, "pcr_pos", "sero")),
                 values_to = "value",
                 names_to = "state")
  
  # And aggregate
  long_sim %>% 
    {
      x <- .
      if (by_group) {
        x %>% group_by(date, .id, state, school, group)
      } else if (by_school) {
        x %>% group_by(date, .id, state, school)
      } else {
        x %>% group_by(date, .id, state)
      }
    } %>% 
    summarise(
      value = sum(value, na.rm = T),
      n_tot = n()  # total number of individuals in each group
    ) %>% 
    ungroup() %>% 
    mutate(
      n_tot = case_when(state %in% c("C", "S", "Sprime", "sero") ~ n_tot,
                        T ~ NA_real_)
    )
}

#' computeParamStats
#' This functions computes..
#' @param params 
#' @param U 
#'
#' @return
#'
#' 
computeParamStats <- function(params, U) {
  # Get unique parameters
  param_names <- names(params)
  param_cat <- str_remove(param_names, str_c(rev(1:U), collapse = "|"))
  u_param_cat <- unique(param_cat)
  # Compute mean and sd for each parameter
  
  map_df(u_param_cat, function(x) {
    val <- params[param_cat == x]
    tibble(param = x,
           mean = mean(val),
           sd = sd(val),
           q5 = quantile(val, .05),
           q95 = quantile(val, .95))
  })
}

expandParams <- function(df, 
                         params, 
                         U) {
  resmat <- df %>% select(-one_of(params))
  
  for (p in params) {
    tmpmat <- df %>% select(-one_of(p))
    tmpcol <- df %>% select(one_of(p))
    
    resmat <- bind_cols(resmat,
                        map_dfc(1:U, function(x) {
                          tmpcol %>% 
                            magrittr::set_colnames(str_c(p, x))
                        }))
  }
  
  resmat
}

#' computeBlockLikStats
#' Compute block-wise log-likelihood mean and variance from ibpf
#' 
#' @param bpf_list 
#'
#' @return
#' @export
#'
#' @examples
computeBlockLikStats <- function(bpf_list) {
  bpf_bll <- sapply(bpf_list,function(x)x@block.cond.loglik)
  dim(bpf_bll) <- c(dim(bpf_list[[1]]@block.cond.loglik), length(bpf_list))
  bcll <- apply(bpf_bll, c(1,2), mean)
  bcllv <- apply(bpf_bll, c(1,2), var)
  
  list(
    bcll_mean = bcll,
    bcll_var = bcllv
  )
}



coolingFraction <- function(a, N, M) {
  a^((N-1 + ((1:M)-1)*N)/(50*N))
}


dateToDay <- function(date, 
                      ref_date = "2021-01-01") {
  difftime(date, ref_date, units = "days") %>% as.numeric()
}

dateToTime <- function(date, 
                       ref_date = "2021-01-01") {
  dateToDay(date, ref_date)/365
}


dayToDate <- function(day,
                      ref_date = "2021-01-01") {
  as.Date(ref_date) + round(day)
}

timeToDate <- function(x,
                       ref_date = "2021-01-01") {
  
  dayToDate(x * 365, ref_date = ref_date) 
}

pompTrajToDf <- function(x) {
  times <- c(0, unique(fake_data$date)) %>% sort()
  vars <- dimnames(x)$variable
  
  a <- map_df(1:dim(x)[3], function(i) {
    res <- x[,1,i]
    tibble(var = names(res),
           val = res,
           time = times[i]) %>% 
      mutate(
        id = map_chr(var, ~str_extract(., "[0-9]+")),
        state = str_remove(var, id))
  })
  a
}



pompSatesToDf <- function(x) {
  N <- length(x)
  times <- names(x) %>% as.numeric()
  map_df(1:N, function(i) {
    states <- x[[i]]
    states %>% 
      as.data.frame() %>% 
      as_tibble() %>% 
      mutate(state = rownames(states)) %>% 
      pivot_longer(cols = contains("V"),
                   names_to = "particle",
                   values_to = "val") %>% 
      mutate(particle = str_remove(particle, "V") %>% as.numeric(),
             time = times[i],
             id = str_extract(state, "[0-9]+") %>% as.numeric(),
             state = str_remove(state, as.character(id)))
  })
}

getRunDateBounds <- function(terms) {
  
  if ("20/21" %in% terms) {
    t_min <- as.Date("2021-03-11")
  } else if ("21/22" %in% terms) {
    t_min <- as.Date("2021-08-31")
  } else {
    stop("Pass valid terms")
  }
  
  if ("21/22" %in% terms) {
    # This excludes the last visit at Decouverte
    t_max <- as.Date("2022-02-28")
    # t_max <- as.Date("2022-06-30")
  } else {
    t_max <- as.Date("2021-06-30")
  }
  
  c("t_min" = t_min, "t_max" = t_max)
}

filterRunData <- function(data,
                          terms) {
  
  run_bounds <- getRunDateBounds(terms = terms)
  
  # Data that matches date filter
  sub_data <- data %>% 
    filter(date >= run_bounds["t_min"], date <= run_bounds["t_max"])
  
  sub_data
}

filterRunDataForPomp <- function(data,
                                 terms) {
  
  run_bounds <- getRunDateBounds(terms = terms)
  
  # Data that matches date filter
  sub_data <- data %>% 
    filter(date >= run_bounds["t_min"], date <= run_bounds["t_max"])
  
  # Previous data to build initial values with fake serological observation
  prev_data <- data %>% 
    filter(date < run_bounds["t_min"]) %>% 
    replace_na(list(sero_result = 0,
                    n_pos = 0))
  
  if (nrow(prev_data) > 0) {
    sub_data <- sub_data %>% 
      bind_rows(
        prev_data %>% 
          group_by(participant) %>%
          summarise(numpos = sum(n_pos, sero_result, na.rm = T) > 0,
                    n = n()) %>% 
          filter(numpos > 0) %>% 
          select(participant, 
                 sero_result = numpos) %>% 
          # Set time to t0 in pomp model
          mutate(date = run_bounds["t_min"],
                 time = dateToTime(run_bounds["t_min"]),
                 n_pcr = 1, 
                 pcr_pos = 0) %>% 
          filter(!is.na(sero_result))
      )
  }
  
  sub_data
}


getFixedParams <- function(base_pars, 
                           terms) {
  
  if ("20/21" %in% terms & "21/22" %in% terms) {
    fixed_params <- base_pars[str_detect(base_pars, "wt")]
  } else if ("20/21" %in% terms) {
    fixed_params <- base_pars[str_detect(base_pars, "_d|_o")]
  } else if ("21/22" %in% terms) {
    fixed_params <- base_pars[str_detect(base_pars, "wt|_a")]
  } else {
    stop("Provide valide terms")
  }
  
  # Assume no current infection at beginning
  c(fixed_params, "S_0")
}


#' getSearchResults
#'
#' @param search_res_dir 
#' @param redo 
#'
#' @return
#' @export
#'
#' @examples
getSearchResults <- function(search_res_dir,
                             redo = FALSE) {
  
  out_file <- str_c(search_res_dir, "/", 
                    str_remove(search_res_dir, "generated_data/") %>% 
                      str_remove_all("/"),
                    "_compiled.rds")
  
  if (!file.exists(out_file) | redo) {
    
    search_res_files <- dir(search_res_dir, full.names = T) %>% 
      str_subset("res_ibpf", negate = T) %>% 
      str_subset("compiled", negate = T)
    
    res <- map_df(
      1:length(search_res_files), 
      function(x) {
        print(search_res_files[x])
        out <- readRDS(search_res_files[x])
        out$res_df %>% 
          mutate(run = str_extract(search_res_files[x], "[0-9]+(?=\\.rds)") %>% as.numeric(),
                 file = x,
                 path = search_res_files[x])
      })
    
    saveRDS(res, file = out_file)
  } else {
    res <- readRDS(out_file)
  }
  
  res
}


parseFilter <- function(str) {
  str %>% str_split(",") %>% .[[1]] %>% str_trim()
}

getSimpleParname <- function(name){
  if (name == "loglik") {
    name
  } else if (str_detect(name, "lambda")) {
    str_extract(name, "lambdaA_|lambdaB_|lambda_")
  } else if (str_detect(name, "S_|C_")) {
    str_extract(name, "S|C") %>% str_c("0")
  }
}

getVarDict <- function() {
  var_dict <- c("w" = "wildtype",
                "a" = "alpha", 
                "d" = "delta",
                "o" = "omicron")
  
  var_dict
}

#' getMaxLikStarts
#' Get starting parameters based on maximum likelihood fits
#'
#' @param run_file 
#' @param lik_thresh 
#'
#' @return
#'
getMaxLikStarts <- function(hash, 
                            lik_thresh = 15,
                            redo = FALSE) {
  
  # Find previous run results
  search_res_dir <- getResDir(hash = hash,
                              what = "fit")
  
  if (any(str_detect(search_res_dir, "refined"))) {
    search_res_dir <- search_res_dir[str_detect(search_res_dir, "refined")]
  }
  
  cat("-- Getting max likelihood starts from: ", search_res_dir, "\n")
  
  # Corresponding result files, excluding the results of the ibpf objects only
  # as these were being saved separately at some point
  # Load the output parameter vectors
  search_res_df <- getSearchResults(search_res_dir = search_res_dir,
                                    redo = redo)
  
  # Get results within lik_thresh ll units
  max_lik_starts <- search_res_df %>% 
    {
      x <- .
      if (is.null(lik_thresh)) {
        # If no threshold specified no filtering is applyed
        x
      } else {
        filter(x, (ll + (1.96/2) * ll_se) > (max(ll) - lik_thresh)) 
      }
    } %>% 
    # Compute likelihood weights for random sampling
    mutate(weight = exp(-0.5 * (max(ll) - ll)),
           weight = weight/sum(weight)) %>% 
    arrange(desc(ll)) %>% 
    mutate(fit_id = row_number())
  
  cat("-- Returning max likelihood parameters for", nrow(max_lik_starts), "of", 
      nrow(search_res_df), "fits at threshold of", lik_thresh ,"\n")
  
  max_lik_starts
}


getDeltaT <- function() {
  1/365
}


getNonNA <- function(vec) {
  if (any(!is.na(vec))) {
    first(vec[!is.na(vec)])
  } else {
    NA
  }
}

getAgeComposition <- function(ids, 
                              age_cat) {
  
  ages <- filter(age_cat, id %in% idToChr(ids)) %>% 
    count(age_cat) %>% 
    complete(age_cat = c(0, 1)) %>% 
    replace_na(list(n = 0))
  
  
  ages
}



getResDir <- function(hash, 
                      what = "fit") {
  
  dir("generated_data/", 
      pattern = str_c("pomp_", what), 
      full.names = T) %>% 
    str_subset(hash) %>% 
    str_subset("rds", negate = T) %>% 
    str_subset("bckup", negate = T) %>% 
    cleanFileName()
}

cleanFileName <- function(x) {
  x %>% 
    str_replace("//", "/")
}

getRunDataFile <- function(hash) {
  file <- dir("generated_data/", pattern = "run_data", full.names = T) %>% 
    str_subset(hash) %>%
    str_subset("test_data", negate = T)
  
  if (length(file) == 0) {
    stop("No run data file found for hash", hash)
  }
  
  file
}

extractHash <- function(s) {
  s %>% 
    cleanFileName() %>% 
    str_remove("generated_data/") %>% 
    str_remove("pomp_|run_data_") %>%
    str_remove("fit_|traj_|profile_") %>% 
    str_sub(1, 6)
}

makeRunDataFile <- function(filter_schools = "all", 
                            filter_terms = "all",
                            model_version = "gen",
                            use_gen_lik = TRUE) {
  
  # make time hash for tile name
  time_hash <- str_sub(digest::digest(c(now())), 1, 6)
  
  str_glue("run_data_{time_hash}_{filter_schools %>% 
                     str_c(collapse = '-')}_{filter_terms  %>% 
                     str_c(collapse = '-')}_{model_version}_genlik{as.character(use_gen_lik)}") %>% 
    janitor::make_clean_names() %>% 
    str_c("generated_data/", ., ".rds")
}

makePompObjFile <- function(run_file) {
  f <- str_replace(run_file, "run_data", "pomp_obj")
  
  cleanFileName(f)
}

makePompFitFile <- function(run_file,
                            refine = FALSE) {
  f <- str_replace(run_file, "run_data", "pomp_fit")
  
  if (refine) {
    f <- str_replace(f, "\\.rds", "_refined.rds")
  }
  
  cleanFileName(f)
}

makePompTrajFile <- function(run_file,
                             refine = FALSE,
                             outbreak_id = NULL) {
  
  f <- str_replace(run_file, "run_data", "pomp_traj")
  
  if (!is.null(outbreak_id)) {
    f <- str_replace(f, "\\.rds", str_glue({"_states_outbreak_{outbreak_id}.rds"}))
  }
  
  cleanFileName(f)
}

makeTraceplotFile <- function(run_file) {
  str_glue("figures/traceplot_fit_{str_replace(run_file, '.rds', '.png') %>% str_remove('generated_data/')}")
}

makeSimplotFile <- function(run_file) {
  str_glue("figures/simplot_fit_{str_replace(run_file, '.rds', '.png') %>% str_remove('generated_data/')}")
}

makeMarkdownFile <- function(run_file) {
  str_glue("data_repot_{str_replace(run_file, '.rds', '.html') %>% str_remove('generated_data/')}")
}

makePompProfileFile <- function(run_file,
                                param,
                                lower,
                                upper,
                                n_prof,
                                n_rep,
                                maxlik_start) {
  
  str_replace(run_file, "run_data", "pomp_profile") %>% 
    str_remove(".rds") %>% 
    str_c(
      str_glue("_{param}_l{lower}_to_h{upper}_{n_prof}_by_{n_rep}_max_lik_{maxlik_start}"),
      ".rds"
    )
}

makeSharedParams <- function(init_params,
                             fixed_params) {
  names(init_params)[str_detect(names(init_params), "lambda")] %>% 
    setdiff(fixed_params)
}



makeBlockAges <- function(block_list,
                          age_cat) {
  map_df(block_list, 
         ~ getAgeComposition(., age_cat = age_cat) %>% 
           pivot_wider(names_from = "age_cat",
                       values_from = "n")) %>% 
    mutate(block = names(block_list)) %>% 
    inner_join(
      map_df(seq_along(block_list), function(x) {
        tibble(
          participant = idToChr(block_list[[x]]),
          block = names(block_list)[x]
        )
      }) %>% 
        inner_join(age_cat, by = c("participant" = "id"))
    ) %>% 
    mutate(n_block = case_when(age_cat == 0 ~ `0`, 
                               TRUE ~ `1`)) %>% 
    select(participant, age_cat, n_block) %>% 
    mutate(n_block = as.numeric(n_block))
}



makeMapParticipantBlock <- function(block_list) {
  blocks <- as.integer(names(block_list))
  map_to_block <- rep(NA, length(unlist(block_list)))
  
  for (i in seq_along(block_list)) {
    b <- block_list[[i]]
    for (j in seq_along(b)) {
      map_to_block[b[j]] <- blocks[i]
    }
  }
  
  if (any(is.na(map_to_block))) {
    stop("Some participants were not assigned to blocks")
  }
  
  map_to_block
}


#' makeInitParams
#'
#' @param n_start 
#' @param refined 
#' @param fixed_params 
#'
#' @return
#' @export
#'
#' @examples
makeInitParams <- function(n_start,
                           refined,
                           fixed_params,
                           init_params) {
  if (refined) {
    
    # Take initial parameters from reasonable bounds based on previous runs
    init_params_df <- map_df(1:n_start, function(x) {
      tibble(
        lambda_wt = 10^runif(1, -1, 1.5),
        lambdaA_wt = 10^runif(1, -1, 1.5),
        lambdaB_wt = lambdaA_wt * runif(1, .1, .5),
        lambda_a = runif(1, .01, .05),
        lambdaA_a = runif(1, 1, 1.2),
        lambdaB_a = runif(1, 0.05, 0.15),
        lambda_d = runif(1, .2, .5),
        lambdaA_d = runif(1, 1, 1.2),
        lambdaB_d = lambdaB_a * runif(1, 1.5, 1.7),
        lambda_o = runif(1, 5, 8),
        lambdaA_o = runif(1, 1.1, 1.4),
        lambdaB_o = lambdaB_d * runif(1, 1.5, 1.7),
      ) %>%
        mutate(C_0 = runif(1, .3, .5),
               S_0 = .9)
    })
  } else {
    
    # Take initial parameters from reasonable bounds based on previous runs
    init_params_df <- map_df(1:n_start, function(x) {
      tibble(
        lambda_wt = 10^runif(1, -1, 1.5),
        lambdaA_wt = 10^runif(1, -1, 1.5),
        lambdaB_wt = lambdaA_wt * runif(1, .1, .5),
        lambda_a = runif(1, .5, 1.5),
        lambdaA_a = runif(1, .3, .6),
        lambdaB_a = runif(1, 0.01, 0.3),
        lambda_d = runif(1, .3, 1.5),
        lambdaA_d = runif(1, .1, 10),
        lambdaB_d = runif(1, 0.1, 0.3),
        lambda_o = runif(1, 3, 9),
        lambdaA_o = runif(1, 1, 3),
        lambdaB_o = runif(1, 0.1, 1),
      ) %>% 
        mutate(C_0 = runif(1, .3, .4),
               S_0 = .9)
      
    }) 
  }  
  
  # Keep only parameters of interest for this run 
  init_params_df <- init_params_df %>%
    select(-one_of(fixed_params)) %>% 
    bind_cols(
      init_params[c(names(init_params) %in% fixed_params)] %>% 
        enframe() %>% 
        pivot_wider()
    )
  
  init_params_df
}


#' makeBlockList
#'
#' @param full_sugar_ids 
#' @param u_sugar_ids 
#' @param u_sugar_ids_num 
#'
#' @return
#' @export
#'
#' @examples
makeBlockList <- function(full_sugar_ids,
                          u_sugar_ids,
                          u_sugar_ids_num) {
  
  full_sugar_ids %>% 
    inner_join(
      tibble(sugar_id = u_sugar_ids,
             participant = u_sugar_ids_num)
    ) %>% 
    arrange(participant, term) %>% 
    group_by(participant) %>% 
    slice(1) %>% 
    ungroup() %>% 
    group_by(group) %>% 
    group_map(function(x, y) {x$participant}) %>% 
    purrr::map(~ as.numeric(.)) %>% 
    set_names(seq(0, (length(.)-1)))
}


#' mcap2
#' Computation of confidence intervals. Copied from
#' Ionides, Edward L., et al. "Monte Carlo profile confidence intervals for dynamic systems." Journal of The Royal Society Interface 14.132 (2017): 20170126.
#'
#' @param logLik 
#' @param parameter 
#' @param level 
#' @param span 
#' @param Ngrid 
#'
#' @return
#' @export
#'
#' @examples
mcap2 <- function (logLik, 
                   parameter, 
                   level = 0.95, 
                   span = 0.75, 
                   Ngrid = 1000) {
  
  smooth_fit <- loess(logLik ~ parameter, span = span)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_logLik <- predict(smooth_fit, newdata = parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_logLik)]
  dist <- abs(parameter - smooth_arg_max)
  included <- dist <= sort(dist)[trunc(span * length(dist))]
  maxdist <- max(dist[included])
  weights <- numeric(length(parameter))
  weights[included] <- (1 - (dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(logLik ~ a + b, 
                      weights = weights, 
                      data = data.frame(logLik = logLik, 
                                        b = parameter, a = -parameter^2))
  b <- unname(coef(quadratic_fit)["b"])
  a <- unname(coef(quadratic_fit)["a"])
  m <- vcov(quadratic_fit)
  var_b <- m["b", "b"]
  var_a <- m["a", "a"]
  cov_ab <- m["a", "b"]
  se_mc_squared <- (1/(4 * a * a)) * (var_b - (2 * b/a) * 
                                        cov_ab + (b * b/a/a) * var_a)
  se_stat_squared <- 1/2/a
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(level, df = 1) * (a * se_mc_squared + 0.5)
  logLik_diff <- max(smoothed_logLik) - smoothed_logLik
  ci <- range(parameter_grid[logLik_diff < delta])
  list(logLik = logLik, 
       parameter = parameter, 
       level = level, 
       span = span, 
       quadratic_fit = quadratic_fit, 
       quadratic_max = b/(2 *a), 
       smooth_fit = smooth_fit, 
       fit = data.frame(parameter = parameter_grid, 
                        smoothed = smoothed_logLik, 
                        quadratic = predict(quadratic_fit, 
                                            newdata = list(b = parameter_grid, 
                                                           a = -parameter_grid^2))), 
       mle = smooth_arg_max,
       max_ll = max(smoothed_logLik),
       ci = ci, 
       delta = delta, 
       se_stat = sqrt(se_stat_squared), 
       se_mc = sqrt(se_mc_squared), 
       se = sqrt(se_total_squared))
}



idToChr <- function(id) {
  str_pad(id, width = 3, side = "left", pad = "0")
}


#' runSimulations
#'
#' @param pomp_obj 
#' @param nsim 
#' @param delta_t 
#' @param tmin 
#' @param tmax 
#' @param school_groups
#'
#' @return
#' @export
#'
#' @examples
runSimulations <- function(pomp_obj,
                           params = NULL,
                           nsim = 10,
                           delta_t = 1/365,
                           step = 1,
                           tmin,
                           tmax,
                           school_groups = NULL) {
  
  if (!is.null(params)) {
    pomp_obj@params <- params
  } else {
    params <- pomp_obj@params
  }
  
  sim <- spatPomp::simulate(
    pomp_obj,
    params = params,
    nsim = nsim,
    format = "data.frame",
    times = seq(tmin - delta_t, tmax, by = step*delta_t)) %>% 
    mutate(date = timeToDate(time))
  
  
  if (is.null(school_groups)) {
    sim <- sim %>% 
      mutate(school = 1, group = 1)
  } else {
    sim <- sim %>% 
      inner_join(school_groups, by = c("unitname" = "id")) %>% 
      group_by(.id, unitname, time) %>% 
      slice(1) %>% 
      ungroup()
  }
  
  sim
}

addVariantRibbons <- function() {
  list(geom_rect(inherit.aes = FALSE,
                 data = map_df(setdiff(getVarDict(), "wildtype"), ~ getVariantDates(.)),
                 aes(xmin = TL, xmax = TR, ymin = -Inf, ymax = Inf,
                     fill = variant),
                 alpha = .1),
       scale_fill_manual(values = getVariantColors()))
}

#' plotSimulations
#'
#' @param sim 
#'
#' @return
#' @export
#'
#' @examples
plotSimulations <- function(sim,
                            by_school = FALSE,
                            by_group = FALSE,
                            simplify_states = FALSE) {
  
  sim <- sim %>% 
    aggregateSim(seis_unit_statenames = makeStateNames(version = "gen"),
                 by_school = by_school,
                 by_group = by_group,
                 simplify_states = simplify_states)
  
  if (by_group) {
    sim <- sim %>% 
      mutate(school_group = str_c(school, group, sep = "_"))
    
    # Make subplot for each school
    plot_list <- map(
      unique(sim$school), 
      function(x) {
        sim %>% 
          filter(school == x) %>% 
          ggplot(aes(x = date, y = value)) +
          addVariantRibbons() +
          geom_line(alpha = .3, aes(group = .id))  + 
          geom_hline(aes(yintercept = n_tot), lty = 2, lwd = .3) +
          facet_grid(state ~ group, scales = "free") +
          theme_bw() +
          ggtitle(x)
      })
    
    p_sim_gen <- cowplot::plot_grid(plotlist = plot_list,
                                    ncol = 1)
  } else {
    p_sim_gen <- sim %>% 
      ggplot(aes(x = date, y = value)) +
      addVariantRibbons() +
      geom_line(alpha = .3, aes(group = .id)) +
      geom_hline(aes(yintercept = n_tot), lty = 2, lwd = .3) +
      theme_bw()
    
    if (by_school) {
      p_sim_gen <- p_sim_gen + 
        facet_grid(state ~ school, scales = "free")
    } else {
      p_sim_gen <- p_sim_gen + 
        facet_wrap(~ state, scales = "free")
    }
  }
  
  p_sim_gen +
    coord_cartesian(xlim = c(min(sim$date), max(sim$date)))
}

#' plotIPFtraces
#' Plot traces of iterated particle filter
#'
#' @param res_ibpf 
#'
#' @return
#' @export
#'
#' @examples
plotIPFtraces <- function(res_ibpf,
                          params = NULL) {
  
  traces_df <- res_ibpf@traces %>% 
    as.data.frame() %>% 
    mutate(iteration = row_number()) %>% 
    # Mutate lambdas on log scale for plotting
    mutate(across(contains("lambda"), log10)) %>%
    select(iteration, loglik, contains("lambda"), contains("_")) %>% 
    pivot_longer(cols = c(loglik, contains("lambda"), contains("_"))) %>%
    mutate(name_simple = map_chr(name, ~getSimpleParname(.)),
           name_simple = coalesce(name_simple, name),
           name_simple = case_when(
             str_detect(name_simple, "lambda") ~ str_c("log10_", name_simple),
             T ~ name_simple
           ),
           variant = case_when(str_detect(name, "lambda") ~ getVarDict()[str_extract(name, "_[w|a|d|o]") %>% str_remove("_")],
                               T ~ NA_character_)
    )
  
  if (!is.null(params)) {
    traces_df <- filter(traces_df, name_simple %in% params)
  }
  
  ybounds <- traces_df %>% 
    distinct(name_simple) %>% 
    mutate(y_min = case_when(str_detect(name_simple, "log10_lambda") ~ -1.5,
                             name_simple == "loglik" ~ -Inf,
                             T ~ 0),
           y_max = case_when(str_detect(name_simple, "log10_lambda") ~ 1.5,
                             name_simple == "loglik" ~ Inf,
                             T ~ 1)) %>% 
    pivot_longer(cols = contains("y"),
                 names_to = "bound",
                 values_to = "y")
  
  p_fit_gen <- traces_df %>% 
    mutate(alpha = ifelse(str_detect(name, "loglik"), 1, .4)) %>% 
    ggplot(aes(x = iteration, y = value, group = name, color = variant)) +
    geom_point(data = ybounds, inherit.aes = F, aes(x = 0, y = y), alpha = 0) +
    geom_line(aes(alpha = alpha), size = .4) +
    theme_bw() +
    facet_wrap(~name_simple, scales = "free_y")
  
  p_fit_gen
}

perturbInitParam <- function(x, w = .2) {
  n_elem <- length(x)
  x * runif(n_elem, min = 1 - w/2, max = 1 + w/2)
}


#' Title
#'
#' @param run_data 
#'
#' @return
#' @export
#'
#' @examples
makeSchoolGroupMap <- function(run_data) {
  
  run_data$full_sugar_ids %>% 
    select(sugar_id, school, group, term) %>% 
    inner_join(run_data$sugar_ids_mapping, .) %>% 
    arrange(id)
}



#' filterTrajToDF
#'
#' @param filter.traj 
#' @param to_keep 
#' @param simplify 
#'
#' @return
#' @export
#'
#' @examples
filterTrajToDF <- function(filter.traj,
                           to_keep = c("I", "Pone", "Ptwo", "Pthree"),
                           simplify = FALSE) {
  # Dim: nstates x 1 x ntimes
  state_names <- attributes(filter.traj)$dimnames$name
  
  # unpack
  traj <- map_df(
    # Loop over times
    1:dim(filter.traj)[3],
    function(x) {
      tibble(
        variable = state_names,
        value = filter.traj[, 1, x],
        time = x
      ) %>% 
        mutate(compartment = str_remove(variable, "[0-9]+"),
               id = str_extract(variable, "[0-9]+") %>%  
                 idToChr()) %>% 
        filter(compartment %in% to_keep)
    })
  
  traj
}


#' Title
#'
#' @param res_list 
#' @param to_keep 
#' @param school_groups 
#' @param times 
#'
#' @return
#' @export
#'
#' @examples
extractTraj <- function(res_list,
                        to_keep = "C",
                        school_groups,
                        times,
                        time_min = NULL,
                        time_max = NULL,
                        school = NULL,
                        group = NULL) {
  
  trajs <- purrr::map_df(
    seq_along(res_list), 
    function(y) {
      out <- filterTrajToDF(filter.traj = res_list[[y]],
                            to_keep = to_keep) %>% 
        mutate(draw = y, 
               time = times[time])
      
      if (!is.null(time_min)) {
        out <- filter(out, time >= time_min)
      }
      
      if (!is.null(time_max)) {
        out <- filter(out, time <= time_max)
      }
      
      out
    }) %>% 
    inner_join(school_groups, by = "id")
  
  if (!is.null(school)) {
    this_school <- school
    trajs <- filter(trajs, school %in% this_school)
  }
  
  if (!is.null(group)) {
    this_group <- group
    trajs <- filter(trajs, group %in% this_group)
  }
  
  trajs
}

#' Title
#'
#' @param traj_df 
#' @param this_compartment 
#' @param epsilon 
#'
#' @return
#' @export
#'
#' @examples
computeTrajStats <- function(traj_df, 
                             this_compartment,
                             epsilon = 1e-3) {
  
  # determine what function to use to summarise
  if (this_compartment %in% c("hazard", "ratio")) {
    this_fun <- function(x, ...) {
      mean(x)
    }
  } else {
    this_fun <- function(x, n, ...) {
      sum(x)/n
    }
  }
  
  traj_df %>% 
    filter(compartment == this_compartment) %>% 
    group_by(school, age_cat, time, draw) %>%
    summarise(value = this_fun(x = value, n = n())) %>%
    group_by(school, age_cat, time) %>% 
    summarise(mean = mean(value) + epsilon,
              median = median(value) + epsilon,
              q25 = quantile(value, .25) + epsilon,
              q75 = quantile(value, .75) + epsilon,
              q025 = quantile(value, .025) + epsilon,
              q975 = quantile(value, .975) + epsilon) %>% 
    ungroup() %>% 
    mutate(date = timeToDate(time))
  
}


plotBlockStats <- function(block_mean, 
                           block_var,
                           data) {
  
  times <- data %>% pull(time) %>% unique() %>% sort()
  time_cnts <- data %>% count(time) %>% mutate(date = timeToDate(time))
  
  bstats <- map_df(
    1:nrow(block_var), 
    function(x) {
      tibble(block = x,
             var = block_var[x, ],
             mean = block_mean[x, ]) %>% 
        mutate(date = timeToDate(times[row_number()]))
    }) 
  
  cowplot::plot_grid(
    time_cnts %>% 
      ggplot(aes(x = date, y = n)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank()) +
      labs(y = "n observations"),
    bstats %>% 
      ggplot(aes(x = date, y = mean, group = block)) +
      geom_line(alpha = .3) +
      theme_bw() +
      labs(y = "block mean") +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank()),
    bstats %>% 
      ggplot(aes(x = date, y = var, group = block)) +
      geom_line(alpha = .3) +
      theme_bw() +
      labs(y = "block variance") +
      scale_x_date(date_labels = "%Y-%m"),
    align = "v",
    axis = "lr",
    ncol = 1)
  
}




computeProbFromParams <- function(df) {
  df %>% 
    mutate(
      alpha = case_when(variant == "alpha" ~ 2.8,
                        T ~ 3.9),
      mean_prob = case_when(param == "lambda" ~ 1-exp(-mle*7/365),
                            T ~ 1-exp(-mle*alpha*5/365)),
      lo_prob = case_when(param == "lambda" ~ 1-exp(-lo*7/365),
                          T ~ 1-exp(-lo*alpha*5/365)),
      hi_prob = case_when(param == "lambda" ~ 1-exp(-hi*7/365),
                          T ~ 1-exp(-hi*alpha*5/365)),
      what = case_when(param == "lambda" ~ "community",
                       param == "lambdaA" ~ "within-group",
                       T ~ "between-group"),
      where = case_when(param == "lambda" ~ "community",
                        param == "lambdaA" ~ "school",
                        T ~ "school")) 
}


# Adjacency Pomp ----------------------------------------------------------

# Get start and end dates
getStartEndAdj <- function(x, what = "start") {
  res <- x %>% 
    mutate(row = row_number()) %>% 
    group_by(from, year) %>%
    {
      if (what == "start") {
        slice_min(., row) 
      } else {
        slice_max(., row) 
      }
    } %>% 
    ungroup() %>% 
    arrange(year, from) %>% 
    select(from, row, year) %>% 
    pivot_wider(values_from = "row",
                names_from = "year") %>% 
    ungroup() %>% 
    select(-from)  %>% 
    as.matrix()
  
  if(length(unique(x$year)) == 1) {
    if (unique(x$year) == "21/22") {
      res <- cbind(rep(0, nrow(res)), res)
    } else {
      res <- cbind(res, rep(0, nrow(res)))
    }
  }
  
  res[is.na(res)] <- 1000
  res
}


# Genetic data ------------------------------------------------------------

#' getMutationRate
#' Mutation rate of SARS-Cov-2 per day in mutations / site/ year
#' Prior from:
#' Wang, S.,et al., 2022. Molecular evolutionary characteristics of SARS‐CoV‐2 
#' emerging in the United States. Journal of medical virology, 94(1), pp.310-317.
#' https://doi.org/10.1002/jmv.27331
#' 
#' New value from:
#' Kremer, C., et al. (2023), Reconstruction of SARS-CoV-2 outbreaks in a primary 
#' school using epidemiological and genomic data
#' Epidemics, Vol. 44, p. 100701, doi: 10.1016/j.epidem.2023.100701. (Supplementary table S1)
#' 
#' @export
#'
getMutationRate <- function() {
  # Wang et al. 
  # 6.677e-4
  
  # Kremer et al. Mutation rate site/generation
  # (2.5e-6) * 5
  2e-5
}

#' simOutbreak2
#' This is a sligthly modified version of simOutbreak from the outbreakr package
#'
#' @param R0 
#' @param infec.curve 
#' @param n.hosts 
#' @param duration 
#' @param seq.length 
#' @param mu.transi 
#' @param mu.transv 
#' @param rate.import.case 
#' @param diverg.import 
#' @param group.freq 
#' @param spatial 
#' @param disp 
#' @param area.size 
#' @param reach 
#' @param plot 
#' @param stop.once.cleared 
#'
#' @return
#' @export
#'
#' @examples
simOutbreak2 <- function (R0, infec.curve, n.hosts = 200, duration = 50, seq.length = 10000, 
                          mu.transi = 1e-04, mu.transv = mu.transi/2, rate.import.case = 0.01, 
                          diverg.import = 10, group.freq = 1, spatial = FALSE, disp = 0.1, 
                          area.size = 10, reach = 1, plot = spatial, stop.once.cleared = TRUE) 
{
  if (any(group.freq < 0)) 
    stop("negative group frequencies provided")
  group.freq <- group.freq/sum(group.freq)
  K <- length(group.freq)
  R0 <- rep(R0, length = K)
  infec.curve <- infec.curve/sum(infec.curve)
  infec.curve <- c(infec.curve, rep(0, duration))
  t.clear <- which(diff(infec.curve < 1e-10) == 1)
  NUCL <- as.DNAbin(c("a", "t", "c", "g"))
  TRANSISET <- list(a = as.DNAbin("g"), g = as.DNAbin("a"), 
                    c = as.DNAbin("t"), t = as.DNAbin("c"))
  TRANSVSET <- list(a = as.DNAbin(c("c", "t")), g = as.DNAbin(c("c", 
                                                                "t")), c = as.DNAbin(c("a", "g")), t = as.DNAbin(c("a", 
                                                                                                                   "g")))
  seq.gen <- function() {
    res <- sample(NUCL, size = seq.length, replace = TRUE)
    class(res) <- "DNAbin"
    return(res)
  }
  substi <- function(snp) {
    res <- sapply(1:length(snp), function(i) sample(setdiff(NUCL, 
                                                            snp[i]), 1))
    class(res) <- "DNAbin"
    return(res)
  }
  transi <- function(snp) {
    res <- unlist(TRANSISET[as.character(snp)])
    class(res) <- "DNAbin"
    return(res)
  }
  transv <- function(snp) {
    res <- sapply(TRANSVSET[as.character(snp)], sample, 
                  1)
    class(res) <- "DNAbin"
    return(res)
  }
  seq.dupli <- function(seq, T) {
    n.transi <- rbinom(n = 1, size = seq.length * T, prob = mu.transi)
    if (n.transi > 0) {
      idx <- sample(1:seq.length, size = n.transi, replace = FALSE)
      seq[idx] <- transi(seq[idx])
    }
    n.transv <- rbinom(n = 1, size = seq.length * T, prob = mu.transv)
    if (n.transv > 0) {
      idx <- sample(1:seq.length, size = n.transv, replace = FALSE)
      seq[idx] <- transv(seq[idx])
    }
    return(seq)
  }
  choose.group <- function(n) {
    out <- sample(1:K, size = n, prob = group.freq, replace = TRUE)
    return(out)
  }
  if (plot && !spatial) 
    warning("Plot only available with spatial model")
  dynam <- data.frame(nsus = integer(duration + 1), ninf = integer(duration + 
                                                                     1), nrec = integer(duration + 1))
  rownames(dynam) <- 0:duration
  res <- list(n = 1, dna = NULL, onset = NULL, id = NULL, 
              ances = NULL, dynam = dynam)
  res$dynam$nsus[1] <- n.hosts - 1
  res$dynam$ninf[1] <- 1
  res$onset[1] <- 0
  res$id <- 1
  res$ances <- NA
  res$group <- choose.group(1)
  EVE <- seq.gen()
  res$dna <- matrix(seq.dupli(EVE, diverg.import), nrow = 1)
  class(res$dna) <- "DNAbin"
  if (spatial) {
    res$xy <- matrix(runif(n.hosts * 2, min = 0, max = area.size), 
                     ncol = 2)
    res$inf.xy <- res$xy[1, , drop = FALSE]
  }
  
  # Modify function to save statuses
  res$status <- matrix(0, nrow = duration+1, ncol = n.hosts)
  res$status[1,] <- c("I", rep("S", ncol(res$status) - 1))
  
  for (t in 1:duration) {
    # Initialize status
    res$status[t+1,] <- res$status[t,]
    
    indivForce <- infec.curve[t - res$onset + 1]
    if (spatial) {
      res$xy <- disperse(res$xy, disp = disp, area.size = area.size)
      if (plot) {
        myCol <- rep("black", nrow(res$xy))
        myCol[res$status[t,] == "I"] <- "red"
        myCol[res$status[t,] == "R"] <- "royalblue"
        plot(res$xy, pch = 20, cex = 6, col = transp(myCol), 
             main = paste("time:", t), xlab = "", ylab = "")
      }
      k.spa <- .kernel.expo(res$xy, mean = reach)
      spa.force <- apply(k.spa[, res$status[t,] == "S", drop = FALSE], 
                         1, sum)
      spa.force <- spa.force[res$id]
      if (length(indivForce) != length(spa.force)) 
        warning("temporal and spatial forces of infection have different length")
      indivForce <- indivForce * spa.force
      indivForce[indivForce > 1] <- 1
    }
    indivForce <- indivForce * R0[res$group]
    N <- res$dynam$nrec[t] + res$dynam$ninf[t] + res$dynam$nsus[t]
    globForce <- sum(indivForce)/N
    
    if (stop.once.cleared && (globForce < 1e-12)) {
      break
    }
    
    p <- 1 - exp(-globForce)
    nbNewInf <- rbinom(1, size = res$dynam$nsus[t], prob = p)
    if (nbNewInf > 0) {
      res$onset <- c(res$onset, rep(t, nbNewInf))
      newAnces <- sample(res$id, size = nbNewInf, replace = TRUE, 
                         prob = indivForce)
      res$ances <- c(res$ances, newAnces)
      newGroup <- choose.group(nbNewInf)
      res$group <- c(res$group, newGroup)
      if (!spatial) {
        areSus <- which(res$status[t,] == "S")
        newId <- sample(areSus, size = nbNewInf, replace = FALSE)
        res$id <- c(res$id, newId)
        cat("-- Trans: t:", t, "newID:",newId, "\n")
        res$status[t+1,newId] <- "I"
      }
      else {
        for (i in 1:nbNewInf) {
          areSus <- which(res$status[t,] == "S")
          newId <- sample(areSus, 1, prob = k.spa[newAnces[i], 
                                                  areSus])
          res$id <- c(res$id, newId)
          cat("-- Trans: t:", t, "newID:",newId, "\n")
          res$status[t+1,newId] <- "I"
          res$inf.xy <- rbind(res$inf.xy, res$xy[newId])
        }
      }
      newSeq <- t(sapply(match(newAnces, res$id), function(i) seq.dupli(res$dna[i, 
      ], t - res$onset[match(newAnces, res$id)])))
      res$dna <- rbind(res$dna, newSeq)
    }
    nbImpCases <- rpois(1, rate.import.case)
    if (nbImpCases > 0) {
      res$onset <- c(res$onset, rep(t, nbImpCases))
      res$ances <- c(res$ances, rep(NA, nbImpCases))
      newId <- seq(N + 1, by = 1, length = nbImpCases)
      res$id <- c(res$id, newId)
      # Add new cases to matrix
      res$status <- cbind(res$status, matrix("0", nrow = duration + 1, ncol = nbImpCases))
      res$status[1:t, newId] <- "S"
      
      cat("-- Import: t:", t, "newID:",newId, "\n")
      res$status[t+1,newId] <- "I"
      if (spatial) {
        newXy <- matrix(runif(nbImpCases * 2, min = 0, 
                              max = area.size), ncol = 2)
        res$xy <- rbind(res$xy, newXy)
        res$inf.xy <- rbind(res$inf.xy, newXy)
      }
      res$group <- c(res$group, choose.group(nbImpCases))
      newSeq <- t(sapply(1:nbImpCases, function(i) seq.dupli(EVE, 
                                                             diverg.import)))
      res$dna <- rbind(res$dna, newSeq)
    }
    res$status[t+1,res$id[(t - res$onset) >= t.clear]] <- "R"
    res$dynam$nrec[t + 1] <- sum(res$status[t+1,] == "R", na.rm = T)
    res$dynam$ninf[t + 1] <- sum(res$status[t+1,] == "I", na.rm = T)
    res$dynam$nsus[t + 1] <- sum(res$status[t+1,] == "S", na.rm = T)
  }
  
  # keep only sampled
  res$status <- res$status[, res$id]
  
  res$n <- nrow(res$dna)
  res$ances <- match(res$ances, res$id)
  res$id <- 1:res$n
  res$xy <- res$inf.xy
  res$inf.xy <- NULL
  # res$status <- NULL
  findNmut <- function(i) {
    if (!is.na(res$ances[i]) && res$ances[i] > 0) {
      out <- dist.dna(res$dna[c(res$id[i], res$ances[i]), 
      ], model = "raw") * ncol(res$dna)
    }
    else {
      out <- NA
    }
    return(out)
  }
  res$nmut <- sapply(1:res$n, function(i) findNmut(i))
  res$ngen <- rep(1, length(res$ances))
  res$call <- match.call()
  class(res) <- "simOutbreak"
  return(res)
}


# Vacations and quarnateens -----------------------------------------------

#' getGroupColsures
#' this does  this
#' 
#' @return 
#' @export
#'
getGroupColsures <- function() {
  # Masked for data protection
}


#' getVacations
#' Get vacation dates where schools were closed in Geneva.
#'
#' @return
#' @export
#'
#' @examples
getVacations <- function() {
  # Masked for data protection
}


getQuarantinePolicy <- function() {
  
  quarantine_duration <- readxl::read_xlsx("data/03_study_info/SCS_measures_all.xlsx", sheet = 2) %>%
    janitor::clean_names() %>% 
    select(contains("quarantaine")) %>% 
    rename(duration = duree_quarantaine_jours,
           TL = quarantaine_debut,
           TR = quarantaine_fin) %>% 
    mutate(across(c("TL", "TR"), ~ as.Date(., format = "%d.%m.%Y"))) %>% 
    mutate(duration = str_extract(duration, "[0-9]+") %>% as.numeric())
  
  quarantine_duration
  
}


# Plots -------------------------------------------------------------------

getSchoolColors <- function() {
  #  New attempt
  c("#825726", "#677CBF", "#FFCE4A", "#AA996A", "red")
  # paletteer::paletteer_d(5, palette = "ggsci::category20_d3")
}

# getSchoolColors <- function() {
#   c("#91826D", "#E6AC60", "#BA996E", "#66615A", "black")
# }

getSettingColors <- function() {
  c("#050404", "#6E0404", "#C78B8B")
}


getVariantColors <- function() {
  c("#3D8202", "#ED6403", "#5302C4")
  # c("#C2E4AD", "#989E94", "#8D648E")
  # c("orange", "green", "purple")
}



# Helpers for final figures -----------------------------------------------

#' makeSetting
#' Add new column with setting school/pre-school
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
makeSetting <- function(df) {
  df %>% 
    mutate(setting = case_when(school %in% c("SA", "SB") ~ "school", 
                               T ~ "preschool"))
}

#' computeTestStats
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
computeTestStats <- function(df) {
  df %>% 
    summarise(dat_date = median(date),
              n = n(),
              n_pos = sum(result),
              frac = n_pos/n) %>% 
    # Add binomial confidence intervals
    rowwise() %>% 
    mutate(lo = Hmisc::binconf(n_pos, n)[2],
           hi = Hmisc::binconf(n_pos, n)[3]) %>% 
    ungroup()
}

# Statistical model -------------------------------------------------------

#' makeIntervals
#' make observation intervals from long format event data
#' 
#' @param df 
#' @param group_var 
#'
#' @return
#' @export
#'
#' @examples
makeIntervals <- function(df,
                          group_var = "sugar_id") {
  
  # Compute intervals. By default we take interval information based on the
  # left time bound for outbreak id etc...
  intervals <- df %>% 
    group_by_at(group_var) %>% 
    group_modify(function(x, y) {
      x %>% 
        mutate(next_date = lead(date, 1),
               next_time_id = lead(time_id, 1),
               next_school = lead(school, 1),
               next_group = lead(group, 1)) %>% 
        select(outbreak_id, school, group, 
               TL = date, 
               TR = next_date, 
               tid_L = time_id, 
               tid_R = next_time_id) %>% 
        # Only keep intervals that have an end
        filter(!is.na(tid_R)) %>% 
        # Define the id of the interval for this participant
        mutate(part_int_id = row_number())    
    }) %>% 
    ungroup() %>% 
    arrange(!!rlang::sym(group_var), TL) %>% 
    mutate(int_id = row_number())
  
  
  # Covariates for computation of community rates during intervals
  intervals <- intervals %>% 
    mutate(
      # What variant was dominant during the interval
      variant = map_chr(TL, ~ getVariant(.)),
      # What vacations period was it
      vacations = str_extract(outbreak_id, str_c(getVacations()$name, collapse = "|"))) %>%
    # Set in_school if not in vacations
    replace_na(list(vacations = "in_school")) %>% 
    # Set references for covariate matrices
    mutate(variant = factor(variant) %>% forcats::fct_relevel("alpha"),
           vacations = factor(vacations) %>% forcats::fct_relevel("in_school"))
  intervals
}


# Math --------------------------------------------------------------------

inv_logit <- function(x) {
  1/(1+exp(-x))
}

logit <- function(x) {
  log(x/(1-x))
}

log_sum_exp <- function(x) {
  max(x) + log(sum(exp(x-max(x))))
}


# Seroprevalence data -----------------------------------------------------


getSeroprevalenceData <- function() {
  dir("data/09_other_data/", pattern = "HH", full.names = T) %>%
    map_df(function(x) {
      read_csv(x) %>%
        janitor::clean_names() %>%
        filter(str_detect(category, "[0-9]")) %>%
        {
          df <- .
          if (!any(str_detect(colnames(df), "any_antibody_response"))) {
            df %>% rename(c("any_antibody_response" = "seroprevalence_95_percent_ci"))
          } else {
            df
          }
        } %>%
        select(age_group = category, txt = any_antibody_response) %>%
        mutate(mean = str_extract(txt , "(.)*(?=\\()"),
               q5 = str_extract(txt , "(?<=\\()(.)*(?=-)"),
               q95 = str_extract(txt , "(?<=-)(.)*(?=\\))"),
               across(c("mean", "q5", "q95"), function(x) as.numeric(x)/100),
               date = case_when(str_detect(x, "2022") ~ "2022-06-09",
                                str_detect(x, "7352") ~ "2020-12-13",
                                T ~ "2021-07-07") %>%
                 as.Date(),
               school = "general_population") %>%
        select(-txt)
    })  %>%
    mutate(school = factor(school, levels = c("general_population", getAllSchools())),
           cohort = age_group)
}

getSeroprevalenceDates <- function(){
  c(
    "baseline" = "2020-12-13",
    "middle" = "2021-07-07",
    "endline" = "2022-06-09"
  )  
}


getSeroprevalenceFakeDates <- function(){
  c(
    "baseline" = getRunDateBounds("20/21")["t_min"] %>% as.character(),
    "middle" = "2021-07-07",
    "endline" = getRunDateBounds("21/22")["t_max"] %>% as.character()
  )  
}

getSeroprevDraws <- function(n_draws = 100) {
  seroprev_stats <- getSeroprevalenceData()
  
  dates <- getSeroprevalenceDates()
  
  res <- purrr::map(dates, function(x) {
    dat <- filter(seroprev_stats, 
                  date == x,
                  age_group %in% c("[0,6)", "[25,35)")) %>% 
      mutate(sd = (mean - q5)/1.65)
    
    draws <- matrix(NA, nrow = n_draws, ncol = 2)
    
    for (i in 1:nrow(dat)) {
      draws[, i] <- rnorm(n_draws, mean = dat$mean[i], sd = dat$sd[i])
    }
    
    draws
  })
  
  names(res) <- str_c("seroprev_", names(res))
  res
}


# Simulations -------------------------------------------------------------

#' Title
#'
#' @param group_specs 
#'
#' @return
#' @export
#'
#' @examples
makeSimulationIndividuals <- function(group_specs) {
  
  # Loop over individuals and create
  individuals <- map_df(
    seq_len(group_specs$n_groups),
    function(g) {
      
      # Create individuals
      indiv <- map_df(
        c("child", "adult"), 
        function(a) {
          # Define number of individuals to create
          n <- ifelse(a == "child", group_specs$pupils_per_group, group_specs$teachers_per_group)
          
          map_df(
            seq_len(n),
            function(p) {
              tibble(
                # Individual id
                sugar_id = str_c("g", g, "-", a, p),
                # Group
                group = g,
                # Household id
                hh_id = digest::digest(sugar_id) %>% str_sub(1, 5),
                # Age category
                age_cat = a
              )
            })
          
        })
    })
  
  # Reassing households to create families
  is_child <- which(individuals$age_cat == "child")
  new_hh_ids <- sample(is_child, 
                       round(group_specs$household_prob * length(is_child)))
  
  for (id in new_hh_ids) {
    individuals$hh_id[id] <- sample(individuals$hh_id[setdiff(is_child, new_hh_ids)], 1)
  }
  
  
  # Add fake term and school
  individuals <- individuals %>% 
    mutate(term = "20/21",
           school = "S")
  
  individuals
}

#' makeSimRunData
#'
#' @param config 
#'
#' @return
#' @export
#'
#' @examples
makeSimRunData <- function(config) {
  
  # A. Simulation individuals
  sim_individuals <- makeSimulationIndividuals(group_specs = config$group_specs)
  
  # Unique set of sugar ids to use in pomp Model
  u_sugar_ids <- sim_individuals %>% distinct(sugar_id) %>% pull(sugar_id)
  U <- length(u_sugar_ids)    # number of distinct participants
  u_sugar_ids_num <- makeUIdsNum(u_sugar_ids)
  
  sugar_ids_mapping <- tibble(
    sugar_id = u_sugar_ids,
    id = u_sugar_ids_num
  )
  
  
  # B. Adjacency
  adj_list_df <- buildAdjacencyDF(full_sugar_ids = sim_individuals,
                                  u_sugar_ids = u_sugar_ids,
                                  do_checks = FALSE) %>% 
    rename(year = term) %>% 
    arrange(year, from, to, type) %>% 
    distinct() %>% 
    # !! Keep only closest connection
    group_by(from, to) %>% 
    slice_min(type) %>% 
    ungroup() %>% 
    filter(from != to) %>% 
    arrange(year, from, to) %>% 
    mutate(from_chr = from,
           to_chr = to,
           from = map_dbl(from, ~ which(u_sugar_ids_num == .) - 1),
           to = map_dbl(to, ~ which(u_sugar_ids_num == .) - 1))
  
  adj_epi <- list(
    adj_list = select(adj_list_df, from, to),
    adj_type = adj_list_df$type,
    start_list =  getStartEndAdj(adj_list_df, "start") - 1,
    end_list = getStartEndAdj(adj_list_df, "end") - 1
  )
  
  # Check if everyone is connected
  testthat::expect_equal(sum(!(u_sugar_ids_num %in% adj_list_df$from_chr)), 0)
  testthat::expect_equal(sum(!(u_sugar_ids_num %in% adj_list_df$to_chr)), 0)
  
  
  # C. Participant masking
  # Fake information
  masks_info <- list(
    n_masks = rep(0, U),
    start_mask = matrix(0, nrow = U, ncol = 9),
    end_mask = matrix(0, nrow = U, ncol = 9)
  )
  
  
  # D. Phylo data
  gen <-  list(
    do_gen = TRUE,
    adj_list_ancestors = array(dim = 0), 
    adj_list_infectors = array(dim = 0), 
    adj_type_infectors = array(dim = 0), 
    dist_closest_comm_seq = array(dim = 0), 
    closest_comm_seq_len = array(dim = 0), 
    seq_len_mat = array(dim = 0), 
    dist_seq = array(dim = 0), 
    start_list_seq = array(dim = 0), 
    end_list_seq = array(dim = 0), 
    n_max_ancestors = array(dim = 0), 
    temp_diff_seq = array(dim = 0),
    closest_comm_seq_temp_diff = array(dim = 0),
    u_seq_all = array(dim = 0),
    w_mat = array(dim = 0),
    mu = array(dim = 0),
    max_kappa = array(dim = 0),
    seq_covar = array(dim = 0),
    time_seq =  rep(999, U),
    seq_delay =  array(dim = 0)
  )
  
  
  # E. Vaccination
  vaccination <- list(
    t_vacc_fake = dateToTime(as.Date("2021-01-30")),
    t_vacc = rep(-1, U)
  )
  
  
  # F. Test data
  # One time point at the start 
  data <- sim_individuals %>% 
    select(sugar_id) %>% 
    mutate(participant = map_chr(sugar_id, ~ u_sugar_ids_num[which(u_sugar_ids == .)]),
           time = dateToTime(config$sim_specs$date_start),
           n_pcr = NA,
           pcr_pos = NA,
           seq = NA,
           n_block = NA,
           age_cat = NA,
           sero = NA) %>% 
    select(-sugar_id)%>% 
    arrange(participant, time) 
  
  
  # G. Add age cat as data
  age_cat <- sim_individuals %>% 
    mutate(participant = map_chr(sugar_id, ~ u_sugar_ids_num[which(u_sugar_ids == .)])) %>% 
    select(participant, age_cat) %>% 
    mutate(age_cat = age_cat == "adult")
  
  
  # H. Covar
  # Load predicted proportion of SARS-CoV-2 variants in Geneva from SP4 analysis
  prob_variants <- readRDS("data/09_other_data/GE_variant_probs.rds") %>% 
    # Combine all Omicron variants
    mutate(variant = str_extract(variant, "Alpha|Delta|Omicron") %>% str_to_lower()) %>% 
    filter(!is.na(variant)) %>% 
    group_by(week, variant) %>% 
    summarise(prob = sum(mean)) %>% 
    ungroup() %>% 
    filter(week > "2021-01-01") %>% 
    mutate(variant = str_c("prob_", variant)) %>% 
    #!! shift forward by 10 days to capture change in transmission
    mutate(week = week + 10)
  
  covar <- prob_variants %>% 
    pivot_wider(values_from = "prob",
                names_from = "variant") %>% 
    mutate(time = dateToTime(week)) %>% 
    select(-week) %>% 
    arrange(time) %>% 
    mutate(participant = "001") %>% 
    select(participant, time, unique(prob_variants$variant))
  
  # Expand to cover all times 
  max_time_data <- max(data$time)
  max_time_covar <- max(covar$time)
  
  covar <- bind_rows(
    covar,
    covar %>% 
      slice_max(time) %>% 
      mutate(time = max_time_data + getDeltaT())
  )
  
  
  testthat::expect_lt(max_time_data, max(covar$time))
  
  # expand to all unit names
  covar <- map_df(sugar_ids_mapping$id,
                  function(x) {
                    covar %>%
                      mutate(participant = x)
                  }) %>%
    select(participant, time , contains("prob")) %>%
    arrange(time, participant)
  
  
  # G. Parameters for sensitivity
  sero_sens <- list(
    t_log_mean = array(dim = 0),
    t_log_sd = array(dim = 0),
    time_sens_coefs = array(dim = 0)
  )
  
  
  # H. Other globals
  test_perf <- list(
    sens_pcr = array(dim = 0),
    spec_pcr = array(dim = 0),
    sens_sero = array(dim = 0),
    spec_sero = array(dim = 0)
  )
  
  liks <- list(
    use_test_lik = 0,
    use_sero_lik = 0,
    use_gen_lik = 0
  )
  
  # Get natural history parameters
  nat_hist <- getNatHistParams()
  
  
  # Seroprev at baseline and endline
  seroprev <- list(
    seroprev_baseline = array(dim = 0),
    seroprev_middle = array(dim = 0),
    seroprev_endline = array(dim = 0)
  )
  
  globals_parlist <- list(
    test_perf = test_perf,
    adj_epi = adj_epi,
    gen = gen,
    liks = liks,
    nat_hist = nat_hist,
    masks_info = masks_info,
    vaccination = vaccination,
    sero_sens = sero_sens,
    seroprev = seroprev
  )
  
  # I. Block list
  block_list <- makeBlockList(full_sugar_ids = sim_individuals,
                              u_sugar_ids = u_sugar_ids,
                              u_sugar_ids_num = u_sugar_ids_num)
  
  
  
  # Map from units to blocks
  map_participant_block <- makeMapParticipantBlock(block_list = block_list)
  
  # Update
  globals_parlist$seroprev$map_participant_block <- map_participant_block
  globals_parlist$seroprev$map_participant_age_cat <- age_cat$age_cat
  
  
  
  # Compile
  run_data <- list(
    adj_list_df = adj_list_df,
    covar = covar,
    sugar_ids_mapping = sugar_ids_mapping,
    full_sugar_ids = sim_individuals,
    globals_parlist = globals_parlist,
    data = data,
    U = U,
    times = data$time,
    block_list = block_list
  )
  
  run_data
}



replicateParams <- function(params) {
  purrr::map(c("wt", "a", "d", "o"), 
             function(x) {
               vec <- params
               names(vec) <- str_c(names(params), "_", x)
               vec
             }) %>% 
    unlist()
}

# Define parameters
applyInterventions <- function(params, 
                               interventions) {
  
  params$lambda <- params$lambda * interventions$community
  params$lambdaA <- params$lambdaA * interventions$within_group
  params$lambdaB <- params$lambdaB * interventions$between_group
  
  unlist(params)
}


#' computeScenarioStats
#'
#' @param scenario_dirs 
#'
#' @return
#' @export
#'
#' @examples
computeScenarioStats <- function(scenario_dirs,
                                 out_file = "generated_data/scenario_effect_stats.rds",
                                 redo = FALSE) {
  
  # File of all simulation results to compute differences from
  sims_file <- "generated_data/all_scenario_sims.rds"
  
  if (!file.exists(out_file) | redo) {
    
    cat("-- Extracting scenario simulations. \n")
    
    all_sims <- map_dfr(scenario_dirs, function(x) {
      
      # Get simulation files
      sim_files <- dir(x, full.names = T)
      
      # Compile simulations
      sims <- map_df(sim_files, function(y) {
        
        res <- readRDS(y)
        
        res$agg_sim %>% 
          filter(state == "C") %>% 
          select(date, sim = .id, value, n_tot) %>% 
          mutate(scenario = res$config$name,
                 setting = res$config$group_specs$institution,
                 w = res$config$intervention_specs$within_group,
                 b = res$config$intervention_specs$between_group,
                 is_ref = w == 1 & b == 1)
      })
    })
    
    saveRDS(all_sims, file = sims_file)
    
  } else {
    
    cat("-- Loading pre-computed scenario simulations. \n")
    
    all_sims <- readRDS(sims_file)
  }
  
  
  if (!file.exists(out_file) | redo) {
    
    cat("-- Computing scenario effects. \n")
    
    effect_stats <- map_dfr(scenario_dirs, function(x) {
      
      # Get the simulations for this scenario type
      sims <- all_sims %>% 
        filter(
          str_detect(scenario, 
                     x %>% 
                       str_remove("generated_data/simulations") %>% 
                       str_remove_all("/") %>% 
                       str_c("^", .))
        )
      
      # Compute changes
      ref_sims <- filter(sims, is_ref)
      
      # Loop over scenarios that are not the reference
      map_df(unique(sims$scenario[!sims$is_ref]), function(y) {
        
        this_sims <- filter(sims, scenario == y)
        
        # Loop over simulation dates
        map_df(unique(this_sims$date), function(d) {
          
          # Get the reference and scenario simulations
          ref <- filter(ref_sims, date == d, value > 0)
          this <- filter(this_sims, date == d, value > 0)
          
          # If empty return NA
          if (nrow(ref) == 0 | nrow(this) == 0) {
            return(
              tibble(
                date = d,
                mean = NA,
                median = NA,
                q025 = NA,
                q975 = NA
              )
            )
          }
          
          # Compute pair-wise differences
          changes <- matrix(NA, nrow = nrow(ref), ncol = nrow(this)) 
          
          for (i in 1:nrow(changes)) {
            changes[i, ] <- this$value/ref$value[i]
          }
          
          # Compute stats
          tibble(
            date = d,
            mean = mean(changes, na.rm = T),
            median = median(changes, na.rm = T),
            q025 = quantile(changes, na.rm = T, 0.025),
            q975 = quantile(changes, na.rm = T, 0.975)
          )
        }) %>%
          # Add back scenario info
          bind_cols(
            this_sims %>% 
              select(scenario, setting, w, b) %>% 
              slice(1)
          )
      })
    }) %>% 
      mutate(variant = str_extract(scenario, "alpha_like|delta_like|omicron_like"))
    
    saveRDS(effect_stats, file = out_file)
    
  } else {
    cat("-- Loading pre-computed scenario effects. \n")
    
    effect_stats <- readRDS(out_file)
  }
  
  effect_stats
}


#' formatPct
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
formatPct <- function(x) {
  str_c(formatC(x * 100, digits = 0, format = "f"), "%")
}


# Outbreak fidelity -------------------------------------------------------


getCasesOutbreak <- function(outbreak_id,
                             cumulative = FALSE,
                             by_group = TRUE) {
  
  if (by_group) {
    group_vars <- c("date", "group")
  } else {
    group_vars <- c("date", "school")
  }
  
  if (!cumulative) {
    # Load outbreak data
    res <- getDataOutbreak(ob_ids = outbreak_id) %>% 
      filter(what == "pcr/antigen") %>% 
      group_by_at(group_vars) %>% 
      summarise(n_pos = sum(result),
                n_test = n()) %>% 
      mutate(what = "n_inf") %>% 
      ungroup()
    
  } else {
    
    
    res <- getDataOutbreak(ob_ids = outbreak_id) %>% 
      filter(what == "pcr/antigen") %>% 
      bind_rows(getDataOutbreak(ob_ids = 2) %>% 
                  filter(what == "pcr/antigen") %>% 
                  distinct(sugar_id, group) %>% 
                  mutate(date = as.Date("2021-04-15"),
                         result = 0)) %>% 
      group_by(sugar_id) %>% 
      arrange(sugar_id, date) %>% 
      filter(!is.na(result)) %>% 
      mutate(new_inf = c(0, diff(cummax(result)))) %>% 
      group_by_at(group_vars) %>% 
      summarise(new_inf = sum(new_inf)) %>% 
      group_by_at(group_vars) %>% 
      mutate(cumul_inf = cumsum(new_inf)) %>% 
      mutate(what = "cumul_inf")
    
  }
  
  res
}

computeOutbreakStatsByDate <- function(traj_df,
                                       dates,
                                       ids) {
  
  if (length(dates) != length(ids)) {
    stop("Dates should be the same length of ids")
  }
  
  map_df(seq_along(dates), function(x) {
    traj_df %>% 
      filter(sugar_id %in% ids[[x]]) %>% 
      mutate(date = timeToDate(time)) %>% 
      filter(date %in% dates[x]) %>% 
      computeOutbreakSimStats()
  })
  
}

computeOutbreakSimStats <- function(traj_df) {
  
  # Compute statistics
  traj_stats <- traj_df %>% 
    group_by(draw, time, group, sugar_id) %>% 
    summarise(infected = sum(value[compartment %in% c("I", "Pone", "Ptwo", "Pthree")]) == 1,
              past_infected = value[compartment == "C"] == 1) %>%
    group_by(draw, time, group) %>% 
    summarise(n = n(),
              infected = sum(infected),
              sim_infected = rbinom(1, infected, .9) + rbinom(1, n - infected, 0.01),
              past_infected = sum(past_infected),
              sim_past_infected = rbinom(1, past_infected, .85) + rbinom(1, n - past_infected, 0.01),
    ) %>% 
    select(-infected, -past_infected) %>% 
    pivot_longer(cols = c("sim_infected", "sim_past_infected"),
                 names_to = "what") %>% 
    group_by(time, group, n, what) %>% 
    summarise(mean = mean(value),
              median = quantile(value, 0.5),
              q025 = quantile(value, 0.025),
              q975 = quantile(value, .975)) %>% 
    ungroup() %>% 
    mutate(date = timeToDate(time))
  
  traj_stats
}

computeOutbreakTrajStats <- function(traj_df, 
                                     by_group = TRUE) {
  
  if (by_group) {
    group_vars1 <- c("draw", "group")
    group_vars2 <- c("draw", "group", "time")
    group_vars3 <- c("group", "time", "what")
  } else {
    group_vars1 <- c("draw", "school")
    group_vars2 <- c("draw", "school", "time")
    group_vars3 <- c("school", "time", "what")
  }
  
  # Compute statistics
  traj_stats <- traj_df %>% 
    group_by_at(group_vars1) %>% 
    mutate(start_cumul_inf = sum(value[compartment == "C" & time == min(time)])) %>% 
    group_by_at(group_vars2) %>% 
    summarise(n_inf = sum(value[compartment %in% c("I", "Pone", "Ptwo", "Pthree")]),
              cumul_inf = sum(value[compartment == "C"]) - start_cumul_inf,
              n = length(unique(sugar_id)),
              frac = n_inf/n) %>%
    pivot_longer(cols = c("n_inf", "frac", "cumul_inf"),
                 names_to = "what") %>% 
    group_by_at(group_vars3, n) %>% 
    summarise(mean = mean(value),
              median = quantile(value, 0.5),
              q05 = quantile(value, 0.05),
              q95 = quantile(value, .95),
              q025 = quantile(value, 0.025),
              q975 = quantile(value, .975)) %>% 
    ungroup() %>% 
    mutate(date = timeToDate(time))
  
  traj_stats
}

