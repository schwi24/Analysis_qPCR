#' Layout for 96-well plate
#'
#' This is a tibble with the coordinates for a 96-well plate.
#'
#' @import dplyr
#' @import magrittr
#' @import stringr
#' @export
Plate_96 <- tibble(
  Row = replicate(12, LETTERS[1:8]) %>% as.vector(.),
  Column = t(replicate(8, 1:12)) %>% as.vector(.)
) %>%
  mutate(., Well = str_c(Row, Column, sep = "")) %>%
  arrange(., Row, Column) %>%
  mutate(., Well_Number = 1:96)

#' Layout for 384-well plate
#'
#' This is a tibble with the coordinates for a 384-well plate.
#'
#' @import dplyr
#' @import magrittr
#' @import stringr
#' @export
Plate_384 <- tibble(
  Row = replicate(24, LETTERS[1:16]) %>% as.vector(.),
  Column = t(replicate(16, 1:24)) %>% as.vector(.)
) %>%
  mutate(., Well = str_c(Row, Column, sep = "")) %>%
  arrange(., Row, Column) %>%
  mutate(., Well_Number = 1:384)
