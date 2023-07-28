### _____________________
# Import the raw dataset
#'
#' Import the raw .csv file of the 2020 knotweed's tarping survey data. \cr To avoid errors, please ensure
#' that there are no ";" in the cells of the CSV (e.g. if comments have been made in the table cells) and
#' that strings within the table are written in lower cases, without spaces or any special characters
#' (you can write in English for instance).
#'
#' @return A tibble (i.e. a kind of improved data.frame). For further information on tibbles, please refer to
#' the `tidyverse` or \link[readr]{readr} documentation.
#' @export
#' @importFrom readr read_delim
#' @importFrom here here
#'
#' @examples
#' \dontrun{
#' mydata <- import_raw_data()
#' }
import_raw_data <- function(){
  x <- readr::read_delim(here::here("data", "data_tarping_x.csv"), delim = ";", col_names = TRUE)
  return(x)
}





### ____________________
#' Find binary variables
#'
#' @description Automatically detect binary (boolean) variables, even if there is NAs in the data.
#' @param x Any variable (vector, columns etc.) of any length.
#'
#' @return Logical (i.e. a TRUE or FALSE answer).
#' @export
#' @importFrom stats na.omit
#'
#' @examples
#' xx <- c(0, 0, 1, 0, NA)
#' is_binary(xx) # TRUE
is_binary <- function(x) {
  x0 <- na.omit(x)
  length(unique(x0)) %in% 1:2 && all(x0 %in% 0:1)
} # Asks: is the length of the unique elements of x0 (excluding NAs) is 1 or 2 AND all either 0 or 1?





### ____________________________________________
#' Convert factor values as exact numeric values
#'
#' @description Automatically convert the input `x (factor)` into a `numeric` vector while keeping the
#' exact value of `x`.
#' @details  In \strong{R}, \emph{numeric factors} should be coded as `character` but it is rarely the
#' case because most functions actually prefer `factors`. It's one of the shortcomings of \strong{R}.
#' Consequently, there is no \strong{R} function that converts a \emph{numeric factor} (e.g. a binary
#' variable with only zeros and ones) into a `numeric` vector while keeping the correct value (if you use
#' `as.numeric(myfactor)`, you will get the underlying level codes, not the values as numbers), hence the
#' usefulness of `as.numfactor`.
#'
#' @param x A `factor` variable (may contain NAs).
#'
#' @return A `numeric` vector.
#' @export
#'
#' @examples
#' f <- as.factor(cars$speed)
#' as.numeric(f) # Converts the values into level ranks.
#' as.numfactor(f) # Works!
as.numfactor <- function(x){as.numeric(levels(x))[x]}
