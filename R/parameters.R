#' Pollination model parameter values
#'
#' Expert-based parameters for the pollination model
#' 
#' @format A list with dataframes (or lists) containing:
#' \describe{
#'   \item{poll_names}{Basic information about the species}
#'   \item{distances}{Mean distances for foraging and dispersal to new nesting sites}
#'   \item{av}{Maximum nesting densities per hectare}
#'   \item{growth}{Growth function (growth.func) parameters}
#'   \item{florNestInfo}{Floral cover (first dataframe) and floral attractiveness and nesting quality (second dataframe)}
#'   \item{lfn}{Additional information on land-use types}
#'  } 
#' @references Haeussler J, Sahlin U, Baey C, Smith HG, Clough Y (2017) Predicting pollinator population size and pollination ecosystem service responses to enhancing floral and nesting resources. Ecology and Evolution, 7: 1898-1908.\url{http://dx.doi.org/10.1002/ece3.2765}
"parameters"