#' Springbarley aphid biocontrol model function
#'
#' Uses spatial data on floral and nesting quality to compute flower visitation and population dynamics for wild pollinators for a given year.
#' The function can be used for mutliple years by running a loop and using as input the output of the previous model (see example).
#' useraster=TRUE,cSMD,mask,parameters,landscvals,reclassin_PGL,reclassin_PAC,N00="usedefault"
#'
#' @param useraster whether to use a raster, if FALSE the function assumes a matrix.
#' @param cSMD land-use data, provided as a raster or a matrix.
#' @param mask same format as cSMD but containing binary values, 1 = produce a prediction, 0 = ignore
#' @param parameters list of two dataframes containing parameter values
#' @param reclassin_PGL matrix indicating for a sequence of land-use codes the grasslands (1) and the non-grassland uses (0)
#' @param reclassin_PAC matrix indicating for a sequence of land-use codes the annual land-uses (1) and the permanent ones (0)
#' @param N00 the density of colonizing aphids, does not vary in space, if set to usedefault its value is taken from OTHER_PAR.

#' @return a list containing :
#' \describe{
#'   \item{mc1}{mortality phase 1}
#'   \item{mc2}{mortality phase 2}
#'   \item{ad}{aphid days, calculated as N00*(1-mc1)*(exp((f-mc2)*T)-1)/(f-mc2)}
#'   \item{cd}{crop loss in kg.ha-1, calculated as 12.7*ad^.66}
#'  }

#' @references Jonsson, M., Bommarco, R., Ekbom, B., Smith, H. G., Bengtsson, J., Caballero-Lopez, B., Winqvist, C.
#' and Olsson, O. (2014), Ecological production functions for biological control services in agricultural landscapes.
#' Methods Ecol Evol, 5: 243-252. \url{http://dx.doi.org/10.1111/2041-210X.12149}
#' @examples
#' PRED_PAR <- read.table("C:/a1workmain/ES_synthesis/BC_springbarley/yann/PRED_PAR.txt")
#' OTHER_PAR <- read.table("C:/a1workmain/ES_synthesis/BC_springbarley/yann/OTHER_PAR.txt")
#' parameters.bcsb<-list("PRED_PAR"=PRED_PAR,"OTHER_PAR"=OTHER_PAR)
#' reclassin_PAC<-as.matrix(reclassin_info[,c(1,2,4)])
#' reclassin_PGL<-as.matrix(reclassin_info[,c(1,2,5)])
#' biocontrol_pred<-springbarley_biocontrol(cSMD=cSMD_Sk2020,mask=landsmask,
#' parameters=parameters.bcsb,reclassin_PGL=reclassin_PGL,reclassin_PAC=reclassin_PAC)

springbarley_biocontrol<-function(useraster=TRUE,cSMD,mask,parameters, reclassin_PGL,reclassin_PAC,N00="usedefault"){

if(useraster==TRUE){
	# reclass PAC, PGL
	reclassin_PNAC<-reclassin_PAC
	reclassin_PNAC[,3]<-abs(reclassin_PAC[,3]-1) #augment the values so that NAs can get 0 later
	raster_PNAC <- reclassify(cSMD, reclassin_PNAC, right = NA)
	raster_PGL <- reclassify(cSMD, reclassin_PGL, right = NA)
	raster_PAC <- reclassify(cSMD, reclassin_PAC, right = NA)

	# for edges, we could distribute them too with filter2. divide by 25*25, and substract from the other filter2 output
	# now we assume that uncultivated field amrgins are arable-like.

	# define buffer-brushes and implement filter2
	# PNAC 135, 500, 1500 and PGL 500



	brushsize<-c(135,500,1500)*2/25
	brushsize<-c(11,41,121)
	brushit<-function(x){
	temp<-EBImage::makeBrush(x,shape="disc",step=TRUE)
	return(temp/sum(temp))
	}

	raster_PNAC135<-raster(EBImage::filter2(raster::as.matrix(raster_PNAC),brushit(brushsize[1])))
	raster_PNAC500<-raster(EBImage::filter2(raster::as.matrix(raster_PNAC),brushit(brushsize[2])))
	raster_PNAC1500<-raster(EBImage::filter2(raster::as.matrix(raster_PNAC),brushit(brushsize[3])))
	raster_PGL500<-raster(EBImage::filter2(raster::as.matrix(raster_PGL),brushit(brushsize[2])))

	# reduce the cells to only those we want data for.
	temp<-c(values(raster_PNAC135)[values(mask)==1],
		values(raster_PNAC500)[values(mask)==1],
		values(raster_PNAC1500)[values(mask)==1],
		values(raster_PGL500)[values(mask)==1])
	tempm<-matrix(temp,nrow=length(temp)/4,ncol=4,byrow=F)
	rm(temp);gc()

}

if(useraster==FALSE){

	tempm<-matrix(landscvals,nrow=length(landscvals)/4,ncol=4,byrow=F)

}

  	PRED_PAR <- parameters$PRED_PAR
  	OTHER_PAR <- parameters$OTHER_PAR
	#   convert into specific paramaters
 	pw <- PRED_PAR[,1] # Predation rate during aphid colonization
 	U <- PRED_PAR[,2:4] # Predator units
	P <- PRED_PAR[,5:7] # Predator initial pop  # table S1
  	LU_ind <- PRED_PAR[,8:10] # Predator pop influenced by LU
  	a0 <- PRED_PAR[,11:13]
  	b <- PRED_PAR[,14:16]
  	c <- PRED_PAR[,17:19]
  	Wo <- OTHER_PAR[1,]
  	if(length(N00)==1 & N00=="usedefault") N00 <- OTHER_PAR[2,]
  	m2 <- OTHER_PAR[3,] # Aphid attack rate 2nd GS - Ground
  	m3 <- OTHER_PAR[4,] # Aphid attack rate 3rd GS - Flying
  	f <- OTHER_PAR[5,] # Aphid fecundity (baseline=0.2714, Jonsson et al also tested +/- 50%)
  	maloss <- OTHER_PAR[8,] # Loss due to wheel tracks
  	Wa <- Wo/(1-pw) #
  	N0 <- mean(Wo) #
  	m <- c(mean(pw),m2,m3) # attack rate
  	p0 <- -log(1-m)/colSums(U*P) # 'Attack rate' of 1 U
  	p <- matrix(rep(t(p0),each=9),nrow=nrow(P),ncol=ncol(P))*U # 'Attack rate' of different groups - eq 1 och 2 in the paper, repmat replicate and tile array
  	Pin <- P
  	P0 <- Pin/a0
  	a <- Pin*2-P0

  LU_temp<-LU_ind
  LU_temp[LU_ind==-99]<-1

 ## Relevant Land use for the different groups
  #   Returns predator populations as a function of land use around the focal point and the matrix Pin. Pin is the base-line
  #   densities of predator populations. The function uses a single Gomperz-form (positive or negative) for all enemy groups,
  #   but with different parameters and land use dependence. For the groups that do not have any known relation with land use Pout=Pin.
  #   LU.prop is PNAC (propotion not annualy tilled crops) within 135, 500, 1500 and PGL (proportion grassland) within 500 m radius,
  #   for field margin 0, 2.5 and 5 m

T<-14
	usePin<-1*(LU_ind==-99) # use the intercept when no landscape effect
	usePout<-1-usePin       # complement of usePin
	funmc1<-function(x)(1-exp(t(-p[,1])%*%(Pin[,1]*usePin[,1]+usePout[,1]*(P0[,1]+a[,1]*exp(b[,1]*exp(x[LU_temp[,1]]*c[,1])))))  )
	funmc2<-function(x){
			1-exp(-sum(
			p[,2]*(Pin[,2]*usePin[,2]+usePout[,2]*(P0[,2]+a[,2]*exp(b[,2]*exp(x[LU_temp[,2]]*c[,2]))))+
			p[,3]*(Pin[,3]*usePin[,3]+usePout[,3]*(P0[,3]+a[,3]*exp(b[,3]*exp(x[LU_temp[,3]]*c[,3]))))
			))
	}
# funmc1 and funmc2 correspond to equation 1, and equation 9
if(useraster==TRUE){

mc1<-apply(tempm,1,funmc1)

	raster_mc1<-raster_PNAC135
	values(raster_mc1)<-NA
	values(raster_mc1)[(values(cSMD)%in%c(1:8,12))&(values(GSSGMB_rast)==1)]<-mc1

mc2<-apply(tempm,1,funmc2)

	raster_mc2<-raster_PNAC135
	values(raster_mc2)<-NA
	values(raster_mc2)[(values(cSMD)%in%c(1:8,12))&(values(GSSGMB_rast)==1)]<-mc2
ADcalc<-function(x,y){N00*(1-x)*(exp((f-y)*T)-1)/(f-y)}
raster_ad<-N00*(1-raster_mc1)*(exp((f-raster_mc2)*T)-1)/(f-raster_mc2) #equations 5 (2 and 3)
raster_cd<-12.7*raster_ad^.66 # crop loss in kg.ha-1  # equation 6

output<-stack(raster_ad,raster_cd,raster_mc1,raster_mc2)
names(output)<-c("ad","cd","mc1","mc2")
return(output)
}

if(useraster==FALSE){
mc1<-apply(tempm,1,funmc1)
mc2<-apply(tempm,1,funmc2)
ADcalc<-function(x,y){N00*(1-x)*(exp((f-y)*T)-1)/(f-y)}
ad<-N00*(1-mc1)*(exp((f-mc2)*T)-1)/(f-mc2)
cd<-12.7*ad^.66 # crop loss in kg.ha-1

output<-list(ad,cd,mc1,mc2)
names(output)<-c("ad","cd","mc1","mc2")
return(output)
}
#save(raster_ad,raster_cd,raster_mc1,raster_mc2,file=paste("bcsb_out_",scenario,"_iseed",iseed,".RData",sep=""))

}
