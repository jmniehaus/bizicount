#' Nigeria Terrorism Data
#'
#' Data on terrorist attacks by Fulani Extremists and Boko Haram in Nigeria,
#' from the year 2014. Attacks data from Global Terrorism Database, other variables
#' from UCDP PRIO-Grid data.
#'
#' @name terror
#' @docType data
#' @format ## `terror`
#' A data frame with 312 rows and 6 columns:
#' \describe{
#'   \item{att.ful, att.bok}{Integer number of attacks by Fulani Extremists and Boko Haram in 2014.}
#'   \item{xcoord, ycoord}{Longitude and Latidude of grid-cell centroid where attack occurs.}
#'   \item{pop}{Population in grid cell.}
#'   \item{mtns}{Proportion of terrain in grid cell that is considered mountainous.}
#' }
"terror"
