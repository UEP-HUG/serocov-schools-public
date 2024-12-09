# This script contains functions used in the processing of data for
# the stan model



# Data statistics ---------------------------------------------------------

computeProp <- function(df,
                        res_col = "result") {
  df %>% 
    filter(!is.na(!!rlang::sym(res_col))) %>% 
    summarise(n = sum(!is.na(!!rlang::sym(res_col))),
              n_pos = sum(!!rlang::sym(res_col), na.rm = T),
              prop = n_pos/n,
              lo = Hmisc::binconf(n_pos, n)[2],
              hi = Hmisc::binconf(n_pos, n)[3]
    )
}



# Data processing ---------------------------------------------------------


makeEventType <- function(x, 
                          simplify = TRUE,
                          add_terms = TRUE) {
  
  if (!str_detect(x, "-")) {
    id <- ifelse(simplify, "outbreak", str_c("outbreak-", x))
  } else {
    id <- str_extract(x, "baseline|last_visit|debut_school_year|vacations|variant_change|vaccination")
    id <- str_replace(id, "debut", "start")
    if (add_terms) {
      id <- str_c(id, str_extract(x, "20/21|21/22"), sep = "-")
    }
  }
  id
}

makeMap <- function(df, col) {
  s <- rlang::sym(col)
  u_elem <- df %>% distinct(!!s) %>% pull(!!s) %>% sort()
  
  map_dbl(df %>% pull(!!s), 
          function(x){
            which(u_elem == x)
          })
}



# Plots -------------------------------------------------------------------

getEventTypeColors <- function(with_terms = FALSE) {
  if (with_terms) {
    c(
      RColorBrewer::brewer.pal(10, "Paired")[1:6],
      paletteer::paletteer_c(n = length(getAllOutbreaks()), "grDevices::Inferno", -1),
      RColorBrewer::brewer.pal(10, "Paired")[7:10],
      RColorBrewer::brewer.pal(12, "Paired")[11:12]
    ) %>% 
      magrittr::set_names(
        c(
          str_c(rep(c("baseline", "start_school_year","last_visit"), each = 2),
                rep(c("20/21", "21/22"), times = 3),
                sep = "-"),
          str_c("outbreak", getAllOutbreaks(), sep = "-"),
          str_c(rep(c("vacations", "variant_change"), each = 2),
                rep(c("20/21", "21/22"), times = 2),
                sep = "-"),
          str_c(rep(c("vaccination"), each = 2),
                rep(c("20/21", "21/22"), times = 1),
                sep = "-")
        )
      )
  } else {
    c(
      RColorBrewer::brewer.pal(10, "Paired")[c(2, 4, 6)],
      paletteer::paletteer_c(n = length(getAllOutbreaks()), "grDevices::Inferno", -1),
      RColorBrewer::brewer.pal(10, "Paired")[c(8, 10)],
      RColorBrewer::brewer.pal(12, "Paired")[c(11)]
    ) %>% 
      magrittr::set_names(
        c(
          c("baseline", "start_school_year","last_visit"),
          str_c("outbreak", getAllOutbreaks(), sep = "-"),
          c("vacations", "variant_change", "vaccination")
        )
      )
  }
}


# Helpers -----------------------------------------------------------------

make_sens_covar <- function(x, t_log_mean, t_log_sd) {
  
  n <- length(x)
  t_log_mat <- matrix(NA, nrow = n, ncol = 4)
  
  for(i in 1:n) {
    t_log_mat[i, 1] = 1;
    t_log_mat[i, 2] = (log(x[i]) - t_log_mean)/t_log_sd;
    t_log_mat[i, 3] = t_log_mat[i, 2]^2;
    t_log_mat[i, 4] = t_log_mat[i, 2]^3;
  }
  
  t_log_mat
}

#' custom_summaries
#' Custom summaries to get the 95% CrI
#'
#' @return
#' @export
#'
#' @examples
custom_summaries <- function() {
  
  c(
    "mean", "median", "custom_quantile2",
    posterior::default_convergence_measures(),
    posterior::default_mcse_measures()
  )
}

#' cri_interval
#' The Credible interval to report in summaries
#' @return
#' @export
#'
#' @examples
cri_interval <- function() {
  c(0.025, 0.975)
}

#' custom_quantile2
#' quantile functoin with custom cri
#'
#' @param x 
#' @param cri 
#'
#' @return
#' @export
#'
#' @examples
custom_quantile2 <- function(x, cri = cri_interval()) {
  posterior::quantile2(x, probs = cri)
}


