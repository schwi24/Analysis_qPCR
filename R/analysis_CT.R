#' Extract the detectors from a csv files with melting data.
#'
#' This function reads a csv file with qPCR data as exported from ABI 7500.
#' @param file The file to load.
#' @param format The file format that is imported. Currently only csv files
#' from ABI 7500 and ABI 7500 Fast are supported.
#' @param repair_CT_limits A bolean value. If `TRUE`, it will replace
#' `Undetermined` or `NA` values with the next highest out of limit value.
#' @return A tibble with columns named `Well`, `Sample_Name`, `Detector`,
#' `Task`, `CT`, `Std_CT`, `Quantity`, `Mean_Quantity`, `Std_Quantity`,
#' `Filtered`, and  `Melting_Temp`. The columns can be empty.
#' @import dplyr
#' @import magrittr
#' @import stringr
#' @export
read_qPCR <- function(file, format = c("ABI_7500_csv", "Rotorgene_6000_csv", "Rotorgene_6000_rex"), repair_CT_limits = TRUE) {

  data <- read_csv(
    file = file,
    ## Skip over header of file (15 lines) and of the data table (1 line).
    skip = 15L,
    ## Keep filename as "Run".
    id = "Run",
    ## Rename columns.
    col_names = c(
      "Well", "Sample_Name", "Detector", "Task", "CT", "Std_CT", "Quantity",
      "Mean_Quantity", "Std_Quantity", "Filtered", "Melting_Temp"
    ),
    ## CT values as character because of sporadic "Undetermined" value.
    col_types = "ccffcnnnnnn"
  ) %>%
    mutate(., CT = na_if(CT, "Undetermined")) %>%
    #mutate(., CT = recode(.x = CT, !!!c("Undetermined" = "NA"))) %>%
    mutate(., CT = as.numeric(CT))

  if (repair_CT_limits) {

    ## Repair the non-detected CT values

    CT_limit <- c(40, 45, 50, 55, 60)
    test <- data %>%
      select(., CT) %>%
      pull(.) %>%
      max(., na.rm = TRUE) < CT_limit
    CT_limit <- CT_limit[test] %>%
      min(.) + 0.01
    data <- data %>%
      mutate(., CT = replace_na(data = CT, replace = CT_limit))
  }

  return(data)
}

#' Analyze CT values
#'
#' This function applies a function in a subset of the data (e.g. the references
#' assays or the calibrator sample) and performs a calculation with this value
#' on the remaining dataset and stores it as a new column. It is sensitive to
#' groups. This is useful for many calculations in the analysis of qPCR data.
#' @param data A tibble with qPCR data.
#' @param col A character string. The column to apply the function to.
#' @param subset A character string. Select a subset, either `Calibrator` or
#' `Reference`.
#' @param group A character vector. The groups to consider when applying the
#' function on the column.
#' @param funs A character string. Select a function to apply to the column,
#' either `mean`, `mean_centered`, `mean_substract`, `sd`, or `error`.
#' @param new_col_name A character string. Select a name for the new column. If not
#' specified the new column will be named according to the original column `col`
#' and the function `fun`.
#' @return A tibble with a column containing the calculated values.
#' @import dplyr
#' @import magrittr
#' @import stringr
#' @export
scale_group_by_subset <- function(
  data,
  col = "CT",
  subset = c("Calibrator", "Reference"),
  group = c("Detector"),
  funs = c("mean", "mean_centered", "mean_substract", "sd", "error"),
  new_col_name = NULL
) {

  a <- data %>%
    ungroup(.) %>%
    filter(., select(., !!rlang::sym(subset))) %>%
    group_by(select(., !!rlang::sym(group)))

  b <- switch(
    funs,
    "mean" = {
      summarise(a, new_col = mean(!!rlang::sym(col)))
    },
    "mean_centered" = {
      summarise(a, new_col = mean(!!rlang::sym(col)))
    },
    "mean_substract" = {
      summarise(a, new_col = mean(!!rlang::sym(col)))
    },
    "sd" = {
      summarise(a, new_col = sd(!!rlang::sym(col)))
    },
    "error" = {
      summarise(a, new_col = !!rlang::sym(col)) %>%
        unique(.)
    },
    stop("No valid function selected")
  )

  c <- left_join(ungroup(data), ungroup(b), by = group) %>%
    group_by(., select(., !!rlang::sym(group)))

  d <- switch(
    funs,
    "mean" = {mutate(c, new_col = mean(!!rlang::sym(col)))},
    "mean_centered" = {mutate(c, new_col = new_col - (!!rlang::sym(col)))},
    "mean_substract" = {mutate(c, new_col = (!!rlang::sym(col)) - new_col)},
    "sd" = {mutate(c, new_col = sd(!!rlang::sym(col)))},
    "error" = {mutate(c, new_col = sqrt(((!!rlang::sym(col))^2) + new_col^2))},
    stop("No valid function selected")
  )

  if (is_null(new_col_name)) {
    new_col_name <- paste(funs, col, sep = "_")
  }
  e <- d %>%
    ungroup(.) %>%
    rename(., !!new_col_name:=new_col)

  return(e)
}


#' A helper function to create subsets for calibrator sample and reference assay
#'
#' This function adds a column called `name` with boolean values. It is `TRUE`
#' if the column `col` equals the value `value` in the dataset `data`.
#' @param data A tibble with qPCR data.
#' @param col A character string. The column to look for matching values.
#' @param value. A character string. The value to look for matches in column.
#' @param name A character vector. The name of the created column with boolean
#' values.
#' @return A tibble with an additional column containing boolean values.
#' @import dplyr
#' @import magrittr
.set_subset <- function(data, col, value, name){
  data %>%
    ungroup(.) %>%
    mutate(., !!rlang::sym(name):=FALSE) %>%
    mutate(., !!rlang::sym(name):=replace(!!rlang::sym(name), !!rlang::sym(col) == value, TRUE))
}

#' Analyze qPCR data with ddCT method
#'
#' This function analyzes a qPCR dataset in `data` according to the ddCT method.
#' @param data A tibble with qPCR data.
#' @param sample A character string that identifies the sample column.
#' @param detector A character string that identifies the assay/detector column.
#' @param ct A character string that identifies the column with the CT values.
#' @param calibrator A character string with the name of the calibrator sample.
#' @param reference A character string with the name of the reference assay.
#' @param simplify A boolean value that specifies the returned data tibble.
#' If `FALSE` it will return intermediate dCT values.
#' If `TRUE` it will return only ddCT values.
#' @return A tibble with ddCT values.
#' @import dplyr
#' @import magrittr
analyze_ddCT <- function(

  data,
  sample = "Sample",
  detector = "Detector",
  ct = "CT",
  calibrator,
  reference,
  simplify = TRUE

) {

  res <- data %>%

    # Set calibrator sample and reference assay
    .set_subset(., "Sample", calibrator, "Calibrator") %>%
    .set_subset(., "Detector", reference, "Reference") %>%

    # Calculate mean-centered dCT
    scale_group_by_subset(data = ., col = "CT", subset = "Calibrator", group = "Detector", funs = "mean_centered", new_col_name = "dCT") %>%
    group_by(., Sample, Detector) %>%

    # Calculate mean_dCT and sd_dCT (equals sd_CT)
    mutate(.data = ., mean_dCT = mean(dCT)) %>%
    mutate(.data = ., sd_dCT = sd(dCT)) %>%

    # Calculate mean_ddCT and erro_ddCT
    scale_group_by_subset(data = ., col = "mean_dCT", subset = "Reference", group = "Sample", funs = "mean_substract", new_col_name = "mean_ddCT") %>%
    scale_group_by_subset(data = ., col = "sd_dCT", subset = "Reference", group = "Sample", funs = "error", new_col_name = "error_ddCT")

  # Simplify tibble
  if (simplify){
    res <- res %>%
      select(., -dCT, -mean_dCT, -sd_dCT, -Calibrator, -Reference)
  }

  # Return tibble
  return(res)

}
