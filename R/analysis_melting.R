#' Extract the detectors from a csv files with melting data.
#'
#' This function is a helper for the `read_melt` function. It extracts the
#' detectors and their line ranges in the csv file.
#' @param file_name The file to load from.
#' @return A tibble with three columns: `Detector`, `Position`, and `Length`.
#' @import dplyr
#' @import magrittr
#' @import stringr
#' @export
assay_positions <- function(file_name){

  # Extract the detectors and their positions.

  first <- readLines(
    con = file.path("data", file_name),
    n = 1000
  )
  assays_pos <- str_which(first, "^Detector")
  assays_name <- first[assays_pos] %>%
    str_replace(., "^Detector = ", "") %>%
    str_sub(
      .,
      start = 1,
      end = str_locate(., "Reporter = ")[,1] - 2
    )
  bind_rows(
    tibble(Detector = assays_name, Position = assays_pos),
    tibble(Detector = NULL, Position = length(first) + 1)
  ) %>%
    mutate(
      .,
      Length = (lead(Position)-Position - 4)/3
    ) %>%
    drop_na(.)
}

#' Load melting curves.
#'
#' This function reads a csv file with melting curves exported from Applied
#' Biosystems 7500 and 7500 Fast Real-Time PCR System.
#' @param file_name The file to load.
#' @param plate The plate format. The default is `"Plate_96"` but `"Plate_384"` is
#' also acceptable.
#' @return A tibble with the melting data.
#' @import dplyr
#' @import magrittr
#' @export
read_melt <- function(file_name, plate = "Plate_96") {
  # read a '.csv' file with melting data.

  # Initializations
  plt <- tibble(Plate = c("Plate_96", "Plate_384"), Value = 1:2) %>%
    filter(., Plate == plate) %>%
    select(., Value) %>%
    pull(.)
  variables <- c("Temperature", "Raw_Fluorescence", "Derivative")
  assay_pos <- assay_positions(file_name = file_name)
  melt_tbl <- tibble(NULL)

  # Loop over assays (outer) and variables (inner)
  for (i in 1:dim(assay_pos)[1]) {

    assay_i <- assay_pos[i,]
    detector_i <- assay_i %>% select(., Detector) %>% pull(.)

    for (j in 1:3) {

      n_max <- assay_i %>%
        select(., Length) %>%
        pull(.)

      skip <- assay_i %>%
        select(., Position) %>%
        pull(.) %>%
        .[] + (j-1) * (n_max + 1) + 1

      melt_tbl <- bind_rows(
        melt_tbl,
        read_csv(
          file = file.path("data", file_name),
          skip = skip ,
          n_max = n_max,
          id = "Run",
          col_names = c("Well_Number", paste("x", 1:100, sep = "_")),
          col_types = "idddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd",
          show_col_types = FALSE
        ) %>%
          mutate(., Detector = detector_i, Variable = variables[j]) %>%
          #mutate(., Run = str_replace(Run, "^data/", "")) %>%
          pivot_longer(
            data = .,
            cols = starts_with("x_"),
            names_to = "Step",
            names_prefix = "x_",
            values_to = "Value"
          )
      )

    }
  }

  melt_tbl <- melt_tbl %>%
    pivot_wider(
      data = .,
      names_from = "Variable",
      values_from = "Value"
    ) %>%
    left_join(
      .,
      switch(plt, Plate_96, Plate_384),
      by = c("Well_Number"="Well_Number"),
      keep = FALSE
    ) %>%
    select(
      .,
      c(
        "Run", "Well", "Row", "Column", "Detector", "Step", "Temperature",
        "Raw_Fluorescence", "Derivative"
      )
    )

  return(melt_tbl)

}
