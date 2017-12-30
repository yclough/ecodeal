#' Generate random crop rasterstacks
#'
#' Function to generate random crop raster sequences at plot or block level given an existing time series
#' The output us used to run dynamic ecosystem service models (e.g. Häussler et al. 2017)
#' Requires lu.data. The function takes as input a matrix with rows for cells and columns for cell, plot, block, land-use, plot and block edge information
#' The landuse info is for multiple years. The function generates crop frequencies for spatial units and samples from them. Allows to resample at different levels - e.g. plot or block. If random is FALSE than the observed pattern is simply replicated.
#' ls is the spatial scale for sampling (plots or blocks)

#' @param x takes as input a matrix with rows for cells and columns for cell, plot, block, land-use, plot and block edge information
#' @param rid a raster with plot ids
#' @param random If random is FALSE than the observed pattern is simply replicated.
#' @param exclude_field A vector of identifier values to be considered as non-agricultural areas. Defaults to 0.
#' @param exclude_crop A vector of crop codes to be considered as non-crop. Defaults to 14.
#' @param cell.size A number indicating the cell height/widht in meters. Defaults to 25.
#' @param widthgreenedges A number or vector of numbers indicating the assumed width of green field edges in metres. Defaults to 2.
#' @param model One of two final lmer models used in Sirami et al. (model3 or model5). Defaults to model3.
#' @param core logical (TRUE/FALSE), indicates whether the core of the raster should be used in the computation.
#' @param corewidth minimum distance from the centre point to the edge of the core, in units of the projection (usually metres).
#' @return A list containing the estimates, the model call, and a data frame of the explanatory variables
#' @note  if permanent grassland in the plot, then it is grassland in the complete time series. Same for permanent plot and non-agr areas. If we manage at the block-level then the permanent grassland disappears
#' @examples
#' load("landuse.data.RData") # list for each landscape (lu.data) with the landuse per year (lu15 to lu 10) per cell (column)
#' load("plot.raster.RData") # raster with plot codes (rplots)
#' load("lu.codes.RData") #
#' rblocks22<-raster("rblocks22.tif") # blocks - can also be retrieved from lu.data if necessary
#' outrandomplot<-rastergen(x=lu.data[["ls 22"]], rid=rplots[["ls 22"]], random=TRUE, ls="plots", lucols=4:9)
#' outobserved<-rastergen(x=lu.data[["ls 22"]], rid=rplots[["ls 22"]], random=FALSE, ls="plots", lucols=4:9)
#' outrandomblock<-rastergen(x=lu.data[["ls 22"]], rid=rblocks22, random=TRUE, ls="blocks", lucols=4:9)



rastergen <- function(x, rid, random=TRUE, ls, lucols,nyout=6 ,nlu=13){



  if(random==TRUE){

    colplots<-which(colnames(x)==ls)

    # split up the matrix into a list per spatial unit, get lu frequencies over the timeseries, sample a crop sequence for nyout years

      xsplit<-split(x[,lucols],x[,colplots])
      xsplit_names<-names(xsplit)
      lufreq<-do.call(rbind,lapply(xsplit,function(x) tabulate(x,nbins=nlu)/sum(tabulate(x))))
      lufreq<-lufreq[-which(xsplit_names=="0"),]
      cropseq=t(apply(lufreq,1,function(x) sample(x=1:nlu,size=nyout,prob=x,replace=TRUE)))

    # build a rasterstack from the crop sequence and rid

      for(y in 1:nyout){
          croptempy<-subs(x=rid,y=data.frame(as.numeric(rownames(cropseq)),cropseq[,y]),by=1,which=2)
          # previous line substitutes values in rid with the crops, using the rownames to link row to rid
          if(y==1){
            croptemp<-croptempy
          }
          if(y>1){
            croptemp<-stack(croptemp,croptempy)
          }
          rm(croptempy)
      }

      names(croptemp)=paste0("Year ",1:nyout)

  }

  if(random==FALSE){

      for(y in 1:nyout){
        croptempy<-rid
        croptempy[]<-x[,lucols[y]]
        values(croptempy)[values(croptempy)==14]<-NA
        if(y==1){
          croptemp<-croptempy
        }
        if(y>1){
          croptemp<-stack(croptemp,croptempy)
        }
        rm(croptempy)

      }

    names(croptemp)=paste0("Year ",1:nyout)
    warning("will only have produced correct output if x containts the vector of raster values in the right order")

  }
  return(croptemp)
}



