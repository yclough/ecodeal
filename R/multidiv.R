#' Computes multidiversity estimates for land-use crop-cover rasters
#'
#' @param rfield A raster containing identifier values for each agricultural field
#' @param rcrop A raster containing the crop codes
#' @param rgreenedges A raster or Rasterstack with values 0 or 1 indicating in which cells a green field adge is to be found.
#' @param exclude_field A vector of identifier values to be considered as non-agricultural areas. Defaults to 0.
#' @param exclude_crop A vector of crop codes to be considered as non-crop. Defaults to 14.
#' @param cell.size A number indicating the cell height/widht in meters. Defaults to 25.
#' @param widthgreenedges A number or vector of numbers indicating the assumed width of green field edges in metres. Defaults to 2.
#' @param model One of two final lmer models used in Sirami et al. (model3 or model5). Defaults to model3.
#' @param core logical (TRUE/FALSE), indicates whether the core of the raster should be used in the computation.
#' @param corewidth minimum distance from the centre point to the edge of the core, in units of the projection (usually metres).
#' @return A list containing the estimates, the model call, and a data frame of the explanatory variables
#' @examples
#' multidiv(rfield=lu.plot[[1]][[2]],rcrop=lu.plot[[1]][['lu10']],rgreenedges=rbled[[1]],exclude_field=c(0),exclude_crop=c(14),cell.size=25,widthgreenedges=2,model=model3,core=FALSE,corewidth=500)
#'
#' output_multidiv_multiple_landscapes<-do.call("rbind",lapply(as.list(1:13),function(x) multidiv(rfield=lu.plot[[x]][[2]],rcrop=lu.plot[[x]][['lu10']],rgreenedges=rbled[[x]], exclude_field=c(0),exclude_crop=c(14),cell.size=25,widthgreenedges=2,model=model3,core=FALSE,corewidth=500)$pred_multidiv))
#' 
#' output_explvars_multiple_landscapes<-do.call("rbind",lapply(as.list(1:13),function(x) multidiv(rfield=lu.plot[[x]][[2]],rcrop=lu.plot[[x]][['lu10']],rgreenedges=rbled[[x]], exclude_field=c(0),exclude_crop=c(14),cell.size=25,widthgreenedges=2,model=model3,core=FALSE,corewidth=500)$explvars))



multidiv<-function(rfield=lu.plot[[1]][[2]],rcrop=lu.plot[[1]][[3]],rgreenedges=rbled[[1]],
                          exclude_field=c(0),exclude_crop=c(14),cell.size=25,widthgreenedges=2,model=model3,core=TRUE,corewidth=500) {
  
  

  # the following function gets the coordinates for the centrepoints of the rasters to be used as explanaotry variables Lat Lon later on.
  
  
  getrastercentre<-function(x){
    library(spatial.tools)
    if(is.Raster(x)==FALSE) stop("x needs to be a a RasterLayer, RasterBrick, or a RasterStack")
    return(data.frame(Lon=mean(extent(x)[c(1,2)]),Lat=mean(extent(x)[c(3,4)])))
  }
  
  rastercentre<-getrastercentre(rfield)
  rastercentre.spdf= SpatialPointsDataFrame(rastercentre, rastercentre)
  crs(rastercentre.spdf)<-crs(rfield)
  
  if(core==TRUE & is.numeric(corewidth)){
    rfield<-crop(rfield,extent(rgeos::gBuffer(rastercentre.spdf,width=corewidth)))
    rcrop<-crop(rcrop,extent(rgeos::gBuffer(rastercentre.spdf,width=corewidth)))
    rgreenedges<-crop(rgreenedges,extent(rgeos::gBuffer(rastercentre.spdf,width=corewidth)))
  }

    # 4.2.1. Farmland compositional heterogeneity
        # We used crop diversity as a measure of farmland compositional heterogeneity.
        # We measured crop diversity using the Shannon diversity index, a widely used metric of landscape heterogeneity (e.g. Bertrand et al. 2015; Bosem Baillod et al. in press):
        # H^'= - sum over i (p[i] * ln(p[i]) where p[i] is the proportion of crop type i in the landscape.
        # Note that this metric assumes that all agricultural cover types (defined in 4.1) are considered equally different
        # Agricultural cover types included: cereal, fallow, alfalfa, clover, ryegrass, rice, corn, sunflower, sorghum, millet,
        # moha, oilseed rape, mustard, pea, bean, soybean, linseed, orchard, almond, olive, vineyard, mixed vegetables, sugar beet,
        # asparagus, carrot, onion, parsnip, potato, tomato, melon, strawberry, raspberry, wild bird cover, grassland (including temporary
        # and permanent grassland managed for production purpose) and other crops (unknown or rare crops)

        # 1-3: cereal,
        # ? set-aside?

          crop.vals<-as.data.frame(freq(rcrop))
          crop.vals$prop<-crop.vals$count/sum(crop.vals$count[(crop.vals$value %in% exclude_crop)==FALSE ])
          crop.vals$prop[(crop.vals$value %in% exclude_crop)==TRUE ] <- NA

          shdi<- -sum(crop.vals$prop * log(crop.vals$prop),na.rm=T)

    # 4.2.2. Farmland configurational heterogeneity
        # We used mean field size (ha) as a measure of farmland compositional heterogeneity. We chose this metric
        # over total field perimeter length per landscape (e.g. Bosem Baillod et al. in press) because it is directly
        # related to our hypotheses (see Appendix 1). Moreover it is easier to base practical recommendations for future
        # agricultural policies on mean field size rather than on total field perimeter length.

          field.vals<-as.data.frame(freq(rfield))
          mfs<-mean(field.vals$count[(field.vals$value %in% exclude_field)==FALSE ])*(cell.size^2)/10000  # CHECK do only non-agri lus have a plotId == 0?

    # 4.2.3. Semi-natural cover proportion
        # We calculated the sum of woodland (including woody linear elements), open land (e.g. shrubland, grassy margins)
        # and wetland cover (including lakes, rivers, ditches) in the landscape.
          if(class(rgreenedges)[1]=="RasterLayer"){
            propsnh<-
              (
                (sum(field.vals$count[(field.vals$value %in% exclude_field)==TRUE ]) * 25)+
                  (sum(values(rgreenedges)[values(rfield) %in% exclude_field ==FALSE]*widthgreenedges))
              )/(
                sum(field.vals$count)*25
              )
          }         
          if(class(rgreenedges)[1]=="RasterBrick" | class(rgreenedges)[1]=="RasterStack"){
            propsnh<-
              (
                (sum(field.vals$count[(field.vals$value %in% exclude_field)==TRUE ]) * 25)+
                  (sum(values(sum(rgreenedges*widthgreenedges))[values(rfield) %in% exclude_field ==FALSE]))
              )/(
                sum(field.vals$count)*25
              )
          }

    # 4.2.4. Total length of semi-natural linear elements
        # We assessed the total length of semi-natural linear elements between fields (SNL, in meters) by calculating half the sum
        # of all semi-natural linear elements located between two crops. Note that semi-natural linear elements located along roads
        # or urban areas were not included in the calculation of SNL.

        # sum all

        snl <- sum(values(rgreenedges)[values(rfield) %in% exclude_field ==FALSE]) * 25


    # Lat and Long


        rastercentre.spdf_EPSG3395 <- spTransform(rastercentre.spdf,  CRS("+init=epsg:3395"))
        Lat=rastercentre.spdf_EPSG3395@coords[,'Lat']
        Lon=rastercentre.spdf_EPSG3395@coords[,'Lon']

    # others
        # make a SpatialPointsDataFrame (if no data, just coordinates: use SpatialPoints() )


        Year<-2013
        sampled.crop.nb<-2
        Region<-factor(x="Goettingen",levels=levels(model5@frame$Region))


      explvars_raw <-  data.frame(Crop_SHDI = shdi, Crop_MFS = mfs, Seminat_Cover = propsnh,
                   snl = snl, Lat = as.numeric(Lat) ,
                   Lon = as.numeric(Lon) ,sampled.crop.nb = sampled.crop.nb,
                   Year = Year, Region = Region )

      pred_multidiv <- predict(model,newdata=explvars_raw)


    # use the predict function
    return(list(pred_multidiv = pred_multidiv , model=model@call, explvars=explvars_raw ))

}


