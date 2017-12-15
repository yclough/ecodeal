#' Multidiversity model 5
#'
#' A model estimating multidiversity from a range of landscape descriptors. Multidiversity was calculated as follows:
#' For each taxonomic group (plants, bees, butterflies, hoverflies, carabids, spiders and birds), we calculated species richness at the landscape level, i.e. across all three sampling sites and across all visits when multiple survey visits were conducted. We standardized species richness for each taxonomic group by centering species richness on its mean and scaling it based on its standard deviation across all regions.
#' We then calculated the average standardized species richness across all seven taxonomic groups (Sirami et al. in prep).
#'
#' @format A model with a response and six explanatory variables:
#' \describe{
#'   \item{Response}{Multidiversity, defined as the average standardized species richness across seven taxonomic groups}
#'   \item{Crop_SHDI}{Crop diversity, using the Shannon diversity index}
#'   \item{Crop_MFS}{Crop mean field size, in hectares}
#'   \item{Seminat_Cover}{proportion seminatural habitat cover in the landscape}
#'   \item{sampled.crop.nb}{number of crops sampled}
#'   \item{Region}{Factor indicating the study region}
#'   \item{Year}{Sampling year, 2013 or 2014}
#'   ...
#' }
#' @source \url{http://www.farmland-biodiversity.org/}
"model5"
