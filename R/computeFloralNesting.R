#' Compute floral and nesting resource values for flower-visiting insects
#' 
#' Function that computes the floral and nesting arrays according to the bees that are accounted for
#' @param landuseMap the map of the landscape containing the landuse categories, given as a raster
#' @param edgesMap the map containing the length of field edges in each cell, given as a raster or a rasterstack in the case of multiple edge types
#' @param landuseCodes the matching between landuse codes and landuse category names
#' @param bees the vector of bee species
#' @param num_floral the number of floral periods
#' @param florNestInfo a list containing the floral and nesting information, with 4 elements: 1) the floral coverage, 2) bumblee bee info, i.e. floral value and attractiveness, 3) honeybee info and 4) solitary bee info
#' @param codeEdges a vector containing the landuse codes for the edge types (e.g. grassy field edge or sown wildflower strip)
#' @param cell.size a number in meters corresponding to the widht and breadth of the input raster cell
#' @return a list containing two arrays: one of dimension nb of bees containing the nesting qualities, and one of dimension nb of bees * nb of periods containing the floral values for each species and each period

computeFloralNesting <- function(landuseMap, edgesMap, unitEdges="m", widthEdges, landuseCodes, bees, num_floral, florNestInfo,codeEdges,cell.size=25)
{

  require(raster)
  require(plyr)
  nr <- nrow(landuseMap)
  nc <- ncol(landuseMap)

  num_bees <- length(bees)

  # Creation of arrays
  emptyraster<-edgesMap[[1]]
  values(emptyraster)<-NA
  nest <- mget(rep("emptyraster",num_bees))
  floral_periods <- stack(mget(rep("emptyraster",num_floral)))
  names(floral_periods)=paste("floral period",1:num_floral,sep=" ")
  floral <- mget(rep("floral_periods",num_bees))
  names(floral)=bees
  names(nest)=bees
  attract=emptyraster
  # width of edges in the landscape
  if(unitEdges=="sqm"){
  propEdges <- edgesMap/(cell.size^2); # blrast contains the length of the edges
  }
  if(unitEdges=="m"){
  propEdges <- edgesMap*widthEdges/(cell.size^2); # blrast contains the length of the edges
  }
  # Add coverage due to edges
  floralCoverage <- floral_periods

	# to do -> soft code period number by introducing a loop:

 for(v in 1:num_floral){
  values(floralCoverage[[v]]) <- mapvalues(values(landuseMap),florNestInfo$floralCover[,'code'],florNestInfo$floralCover[,paste0('Flor_Cov_P',v,'_b')], warn_missing = FALSE)
  }

  # Derive nesting and floral values for the different species #
  for(s in 1:num_bees)
  {
	  values(nest[[s]]) <- mapvalues(values(landuseMap),florNestInfo$attract[florNestInfo$attract[,'species']==s,'code'],florNestInfo$attract[florNestInfo$attract[,'species']==s,'Nest_P1_b'], warn_missing = FALSE)
  	nest[[s]]<-nest[[s]] * (1-sum(propEdges)) + sum(propEdges * florNestInfo$attract[pmatch(codeEdges,florNestInfo$attract[,'code']),'Nest_P1_b'])
	values(nest[[s]]) <- pmax(values(nest[[s]]),0)

    for (v in 1:num_floral)
    {
	    values(attract) <- mapvalues(values(landuseMap),florNestInfo$attract[florNestInfo$attract[,'species']==s,'code'],florNestInfo$attract[florNestInfo$attract[,'species']==s,paste0("Flor_P",v,"_b")], warn_missing = FALSE)
      attractEdges <- florNestInfo$attract[pmatch(codeEdges,florNestInfo$attract[,'code']),paste0("Flor_P",v,"_b")]

      values(floral[[s]][[v]]) <- pmax(values(floralCoverage[[v]]* attract * (1 - sum(propEdges)) + sum(propEdges * florNestInfo$floralCover[pmatch(codeEdges,florNestInfo$floralCover[,'code']),paste0("Flor_Cov_P",v,"_b")] * attractEdges)),0)
    }
  }

  # check that nesting and floral values are different from NA
  testfunction=function(x,type){
	if(sum(is.na(values(x))>0)) warning(paste0("NA(s) detected in ",type," values"))
	}
  lapply(nest,testfunction,type="nesting")
  lapply(floral,testfunction,type="floral")

  return(list(nest=nest,floral=floral))
}
