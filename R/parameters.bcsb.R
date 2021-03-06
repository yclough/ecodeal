#' Spring barley aphid biocontrol model parameter values
#'
#' Parameter values as used in the model. There are three stages: colonizing, growing ground-living and growing flying.
#' PRED_PAR contains the predation rate during aphid colonization (1st column),
#' predator units (columns 2-4), predator initial populations (columns 5-7), predator populations as influenced by the
#' land-use (columns 8:10), coefficients a0 (columns 11-13), b (columns 14-16), c (columns 17-19).
#' OTHER_PAR contains Wo, N00, m2, m3, f, maloss (all in their individual columns).
#' The number -99 means that the landscape does not influence population size during that stage.
#' The predator groups are, in order of appearence in the data tables Wolfspiders,Sheet-web spiders,
#' Large ground beetles,Carabidae, large, small ground beetles,Rove beetles, Adult lady beetles,
#' Juvenile lady beetles, Lacewing larvae,Syrphid larvae)

#' p0 are the attack rates calculated from U predatory units, and the othe parameters corresponding to 1 U;
#' p Attack rate accounting for the U of each predator group
#'
#' @format A list with dataframes (or lists) containing:
#' \describe{
#' \item{PRED_PAR}{Predator parameters, Predation rate during aphid colonization (see above)}
#'   \describe{
#'     \item{pw}{Predation rate during aphid colonization, calculated from mortality m=(N00-N0)/N00 in Östman et al. 2003}
#'     \item{U}{Predator units of the enemies during the different stages}
#'     \item{P}{Baseline predator densities during the different stages}
#'     \item{LU}{relevant land-use for the different groups and periods (1=PNAC135; 2=PNAC500; 3=PNAC1500; 4=PGL500)}
#'     \item{a0}{Used in the calculation of P0 (Phat in Jonsson et al) and Gompertz function coefficient a (Jonsson et al, eqn. 9)}
#'     \item{b}{Gompertz function coefficient (Jonsson et al, eqn. 9)}
#'     \item{c}{Gompertz function coefficient (Jonsson et al, eqn. 9)}
#'
#'   \item{OTHER_PAR}{Other parameters (see above) }
#'   \describe{
#'     \item{W0}{N0 from växtskyddscentralens barley fields (=3.6) }
#'     \item{N00}{N00 is the size of the colonizing population, before any predation. (=3.49) }
#'     \item{m2}{Mortality in the second phase (=0.027), estimated from the reduction in r }
#'     \item{m3}{Mortality in the third phase (=0.054), estimated from the reduction in r }
#'     \item{f}{Fecundity (per day), 1.9/7=.2714 from Dean (1974) }
#'     \item{T}{Duration of the growth period, in days (=14)}
#'     \item{th}{Threshold to apply insecticide, in number of aphids per tillar (=5)}
#'     \item{maloss}{Yield loss due to machinery driving over the crop when spraying in percent (=0.002)}
#'     \item{time_limit}{Time after which it is not worth to spray, in days (=10)}
#'     }
#'  }}
#' @references Jonsson, M., Bommarco, R., Ekbom, B., Smith, H. G., Bengtsson, J., Caballero-Lopez, B., Winqvist, C.
#' and Olsson, O. (2014), Ecological production functions for biological control services in agricultural landscapes.
#' Methods Ecol Evol, 5: 243-252. \url{http://dx.doi.org/10.1111/2041-210X.12149}
#' Östman, Ö., Ekbom, B., & Bengtsson, J. (2003). Yield increase attributable to aphid predation by ground-living polyphagous
#'  natural enemies in spring barley in Sweden. Ecological Economics, 45(1), 149-158. \url{https://doi.org/10.1016/S0921-8009(03)00007-7}
#' Dean, G. J. (1974). Effect of temperature on the cereal aphids Metopolophium dirhodum (Wlk.), Rhopalosiphum padi (L.) and
#' Macrosiphum avenae (F.)(Hem., Aphididae). Bulletin of Entomological Research, 63(3), 401-409. \url{https://doi.org/10.1017/S0007485300040888}
"parameters.bcsb"
