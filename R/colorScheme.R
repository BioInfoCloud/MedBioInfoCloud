#' Predefined Color Schemes for Data Visualization
#'
#' A list of pre-defined color palettes for consistent visualization across analyses.
#' Each element is a vector of hex color codes, tailored for different numbers of groups.
#'
#' @format A list with 13 elements:
#' \describe{
#'   \item{color_set1}{4-color palette for general use.}
#'   \item{color_set14}{13-color palette (extended for multiple groups).}
#'   \item{color_set5}{5-color palette for medium-sized group comparisons.}
#'   \item{color_set2}{2-color palette for binary comparisons.}
#'   \item{color_set6}{7-color palette (includes diverse hues).}
#'   \item{color_set2_2}{Alternative 2-color palette.}
#'   \item{color_set2_3}{Another alternative 2-color palette.}
#'   \item{color_set8}{8-color palette for multi-group analyses.}
#'   \item{color_set8_2}{Alternative 8-color palette (softer tones).}
#'   \item{color_set14_2}{15-color palette (extended with varied hues).}
#'   \item{color_set4}{4-color palette (high contrast).}
#'   \item{color_set3}{3-color palette for tri-group comparisons.}
#'   \item{color_set7}{7-color palette (vibrant tones).}
#' }
#' @export
colorScheme = list(
  color_set1 = c("#F3AE63", "#93CBAE", "#73558B", "#E2BECB"),
  # 4-color palette for general visualization
  color_set14 = c(
    "#D2EBC8",
    "#3C77AF",
    "#7DBFA7",
    "#AECDE1",
    "#EE934E",
    "#D1352B",
    "#9B5B33",
    "#F5CFE4",
    "#B383B9",
    "#8FA4AE",
    "#FCED82",
    "#F5D2A8",
    "#BBDD78"
  ),

  # 5-color palette for medium group counts
  color_set5 = c("#D55640", "#479D88", "#6CB8D2", "#415284", "#E69F84"),

  # 2-color palette for binary comparisons
  color_set2 = c("#AA3538", "#5891BF"),

  # 7-color palette with diverse hues
  color_set6 = c(
    "#C6307C",
    "#D0AFC4",
    "#4991C1",
    "#89558D",
    "#AFC2D9",
    "#435B95",
    "#79B99D"
  ),

  # Alternative 2-color palette
  color_set2_2 = c("#E9E55A", "#88558D"),

  # Another 2-color palette (high contrast)
  color_set2_3 = c("#C6367A", "#B3B5B1"),

  # 8-color palette for multi-group analyses
  color_set8 = c(
    "#57AB7E",
    "#46976D",
    "#4A92A4",
    "#95538B",
    "#8BAAACD",
    "#1A1919",
    "#8E9296",
    "#3F7689"
  ),

  # Alternative 8-color palette (softer tones)
  color_set8_2 = c(
    "#DBC9B3",
    "#EED0E0",
    "#EBAEA9",
    "#CB95BB",
    "#EBCC96",
    "#AED0DF",
    "#CBE5DE"
  ),

  # 15-color extended palette
  color_set14_2 = c(
    "#EFE2AA",
    "#95A6DA",
    "#83B4EF",
    "#BFA6C9",
    "#DBC9B3",
    "#F5E0BA",
    "#8ECFF8",
    "#AED0DF",
    "#7EC0C3",
    "#89B780",
    "#EED0E0",
    "#F5D8D0",
    "#EBAEA9",
    "#CB95BB",
    "#AAD0AC"
  ),

  # 4-color high-contrast palette
  color_set4 = c("#651222", "#173565", "#EA3F35", "#9FC4DB"),

  # 3-color palette for tri-group comparisons
  color_set3 = c("#3C234A", "#ECDC52", "#70BACC"),

  # 7-color vibrant palette
  color_set7 = c(
    "#C4452E",
    "#D4342C",
    "#D99943",
    "#1B3A8F",
    "#B0A568",
    "#3B3B7E",
    "#020001"
  )
)
