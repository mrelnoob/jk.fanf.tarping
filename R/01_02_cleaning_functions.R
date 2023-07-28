### Sourcing previous scripts:
source(file = here::here("R/01_01_prep_functions.R")) # I need to do that because one of the functions
# below uses one of the functions from the sourced R file.

library(magrittr)



### _______________
#' Clean data types
#'
#' @description The function `clean_my_data` automatically imports the \emph{raw_data} of the knotweed
#' tarping survey and cleans by transforming character variables into factors, ordinal variables into
#' ordered factors, and boolean/binary variables into factors.
#' @return A cleaned tibble.
#' @note The variables "plantation" and "age" are ordinal variables and I have thus coded them
#' as such (as an ordered factor). However, it might be preferable, from a statistical point of view,
#' to consider it as a numeric variable. The 2nd solution would be more parsimonious (less levels and
#' thus lighter models) but would assume that intervals between each level (between 0 and 1,
#' between 1 and 2, etc.) are equals when it's not a necessarily true assumption (for "plantation",
#' it's probably not). The 1st solution will cause models to add polynomial terms to my factor levels. \cr
#' I also coded "planned_duration" as an ordinal variable but it there's no problem here because I will
#' probably not use it in statistical analyses.
#'
#' @export
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' my_cleaned_data <- clean_my_data()
#' }
clean_my_data <- function(){
  raw_data <- import_raw_data()
  # Transform character variables into factors, ordinal variables into ordered factors, and boolean/binary
  # variables into factors (+ creation of new variables):
  raw_data$plantation[raw_data$plantation == 1] <- 0
  raw_data$plantation[raw_data$plantation == 2] <- 1

  raw_data %>%
    dplyr::mutate(manager_id = factor(x = manager_id),
                  planned_duration = factor(x = planned_duration, ordered = TRUE),
                  age = factor(x = age, ordered = TRUE)) %>%
    dplyr::mutate(geomem = ifelse(grepl("geomem$", fabric_type) | grepl("tarp$", fabric_type) |
                                    grepl("mixed", fabric_type) | grepl("unknown", fabric_type), 1, 0)) %>%
    dplyr::mutate(geotex = ifelse(grepl("geotex$", fabric_type) | grepl("mixed", fabric_type), 1, 0)) %>%
    dplyr::mutate(tarpfix_pierced = ifelse(c(grepl("*staples*", fabric_fixation) | fabric_fixation == "wired_stakes") &
                                             fabric_fixation != "staples_and_taped_patches", 1, 0)) %>%
    dplyr::mutate(stripsoverlap_ok = ifelse(multi_strips == 0 | strips_overlap > 42.5, 1, 0)) %>%
    dplyr::rename(reg_stripsoverlap = reg_stripoverlaps) %>%
    dplyr::mutate_if(is.character, factor) %>%
    dplyr::mutate_if(is_binary, factor) -> cleaned_data

  return(cleaned_data)
}





# _________________________________________
### Creation of a function that will create separate datasets for each response variable (depending on the
# explanatory variables their models will potentially include):

#' Datasets for models building
#'
#' @description The `model_datasets` function creates the respective reduced datasets that should be used
#' to model each response variables. For instance, if `response.var = "overlaps"`, the function will
#' produce a dataset containing \emph{reg_stripsoverlap} as the \strong{response} variable and a subset of
#' variables to be used as \strong{predictors/explanatory variables} (removing all the variables that
#' should not be used to model \emph{reg_stripsoverlap}'s variations).
#' @param response.var A character string specifying which modelling dataset should be produced (either:
#' effectiveness", "edges", "overlaps", or "tarpedarea"):
#' * "effectiveness" will produce the dataset for the 3 effectiveness evaluation variables (namely \emph{
#' eff_eradication}, \emph{eff_expansion} and \emph{eff_vigour});
#' * "edges" will produce the dataset having for response variable the tarping operations that observed
#' regrowth at the edge of the tarped area (at the end of the operation or during the latest visit to the
#' site);
#' * "overlaps" produce the dataset having for response variable the tarping operations that observed
#' regrowth at strip overlaps;
#' * "tarpedarea" will produce the dataset having for response variable the tarping operations that observed
#' regrowth at somewhere on the area covered by the fabric.
#'
#' @return A tibble.
#' @export
#' @import dplyr stringr
#'
#' @examples
#' \dontrun{
#' eff_data <- model_datasets(response.var = "effectiveness")
#' }
model_datasets <- function(response.var = c("effectiveness", "edges", "overlaps",
                                            "tarpedarea")){

  ppp <- clean_my_data()
  ppp %>%
    dplyr::filter(!stringr::str_detect(operation_type, "crushing_tarping_trial")) %>%
    dplyr::mutate(uprootexcav = ifelse(preparation == "excavation" |
                                         preparation == "uprooting", 1, 0)) -> qqq # Exclude rows which
  # belong to crushing-tarping trials, and creates a new variable called "uprootexcav"!

  ### For the 3 "effectiveness" models
  if (response.var == "effectiveness") {
    qqq <- dplyr::mutate(.data = qqq,
                         high_eff = ifelse(stringr::str_detect(latest_regrowth, stringr::regex("few.", dotall = TRUE)) |
                                             eff_eradication == "1", 1, 0)) %>%
      dplyr::mutate_if(is_binary, factor)
    # With str_detect and regex, I matched all obs. with latest_regrowth that begins with "few"
    # (quite similar to grep)!

    qqq <- dplyr::mutate(.data = qqq,
                         effectiveness = rowMeans(x = qqq[,80:82], na.rm = FALSE))


    tapioca <- qqq[,c("manager_id", "xp_id", "effectiveness", "eff_eradication", "high_eff",
                      "latitude", "longitude", "elevation", "goals",
                      "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness",
                      "granulometry", "obstacles", "flood",
                      "geomem", "geotex", "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex",
                      "pla_geotex", "weedsp_geotex", "other_unknown", "grammage", "thickness",
                      "maxveg", "uprootexcav", "stand_surface", "age", "fully_tarped", "distance", "tarping_duration",
                      "stripsoverlap_ok", "tarpfix_multimethod", "tarpfix_pierced", "sedicover_height", "trench_depth",
                      "pierced_tarpinstall", "plantation", "repairs", "add_control", "add_control_type",
                      "degradation", "pb_fixation", "pb_durability")]
  }

  ### For the "lreg_edges" model
  if (response.var == "edges") {
    qqq <- dplyr::mutate(.data = qqq, lreg_edges =
                           ifelse(c(grepl("*edges*", latest_regrowth) | grepl("*edges*", untarped_regrowth)), 1, 0)) %>%
      dplyr::mutate(.data = qqq, reg_elsewhere =
                      ifelse(reg_staples == 1 | reg_stripsoverlap == 1 | reg_obstacles == 1 |
                               reg_holes == 1 | reg_plantations == 1 | reg_pierced == 1, yes = 1, no = 0)) %>%
      dplyr::filter(fully_tarped == 1) # Only for fully tarped operations!

    tapioca <- qqq[,c("manager_id", "xp_id", "lreg_edges",
                      "latitude", "longitude", "elevation", "slope", "flood", "difficulty_access", "shade", "forest",
                      "ruggedness", "granulometry", "obstacles", "stand_surface",
                      "geomem", "geotex",
                      "maxveg", "uprootexcav",
                      "distance", "stripsoverlap_ok",
                      "tarpfix_multimethod", "tarpfix_pierced", "sedicover_height", "trench_depth", "plantation",
                      "repairs", "add_control", "freq_monitoring",
                      "degradation", "pb_fixation", "pb_durability", "tarping_duration", "reg_elsewhere")]
  }

  ### For the "overlaps" model
  if (response.var == "overlaps") {
    qqq <- dplyr::mutate(.data = qqq, reg_elsewhere =
                           ifelse(reg_staples == 1 | reg_edges == 1 | reg_obstacles == 1 |
                                    reg_holes == 1 | reg_plantations == 1 | reg_pierced == 1, yes = 1, no = 0)) %>%
      dplyr::mutate(stripfix_pierced =
                      ifelse(c(grepl("*staples*", strips_fixation) | grepl("*stakes*", strips_fixation)), 1, 0)) %>%
      dplyr::mutate(stripfix_taped =
                      ifelse(c(strips_fixation == "tape_and_heavyobj" | strips_fixation == "double_tape" |
                                 strips_fixation == "tape" | strips_fixation == "tape_and_staples"), 1, 0)) %>%
      dplyr::filter(multi_strips == 1) %>%
      dplyr::mutate_if(is_binary, factor)

    tapioca <- qqq[,c("manager_id", "xp_id", "reg_stripsoverlap",
                      "latitude", "longitude", "elevation", "slope", "flood", "difficulty_access","shade",
                      "forest", "ruggedness",
                      "granulometry", "obstacles", "stand_surface",
                      "geomem", "geotex",
                      "maxveg", "uprootexcav", "levelling", "fully_tarped", "distance",
                      "strips_overlap", "tarpfix_multimethod", "tarpfix_pierced", "stripfix_pierced", "stripfix_taped",
                      "sedicover_height", "trench_depth", "plantation",
                      "repairs", "add_control", "freq_monitoring", "pb_fixation", "pb_durability",
                      "reg_elsewhere", "tarping_duration")]
  }

  ### For the "tarpedarea" model
  if (response.var == "tarpedarea") {
    qqq <- dplyr::mutate(.data = qqq, lreg_tarpedarea =
                           ifelse(c(grepl("*tarpedarea*", latest_regrowth) |
                                      grepl("*tarpedarea*", untarped_regrowth)), 1, 0))

    tapioca <- qqq[,c("manager_id", "xp_id", "lreg_tarpedarea",
                      "latitude", "longitude", "elevation", "slope", "flood",
                      "difficulty_access", "shade", "forest", "ruggedness", "granulometry", "obstacles",
                      "stand_surface",
                      "geomem", "geotex", "woven_geotex",
                      "maxveg", "uprootexcav", "levelling", "fully_tarped", "distance",
                      "stripsoverlap_ok", "tarpfix_multimethod", "tarpfix_pierced", "pierced_tarpinstall",
                      "sedicover_height", "trench_depth", "plantation",
                      "repairs", "add_control", "freq_monitoring",
                      "pb_fixation", "pb_durability", "pb_trampiercing", "reg_edges", "tarping_duration")]
  }
  return(tapioca)
}
