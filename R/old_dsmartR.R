#dsmartR

###Based on the C5 tree modelling the function determines:
#1. The pixel based probability to each class
#2. The n most probable classes at each pixel
##All outputs are save to file in raster format (It is preferable is the working directory is set to the )
# Class probability maps
# n most probable maps
#R pbjects are saved:
#1. The probability rasters (rasterStack) is requested
#2. The n most probable classes (raster Stack)
#FUnction Requires:
# rLocs: A rasterStack of the rasters generated from the dsmart algorithm
# nprob: the n most probable class maps to be produced
# sepP: logical of whether class probability maps should be produced
# lookup: lookup table produced from dsmart that numerically links soil class codes to a number

dsmartR<- function(rLocs= NULL, nprob = 2, sepP=FALSE, lookup= NULL, cpus=1){
  pb <- txtProgressBar(min=0, max=100, style=3)
  beginCluster(cpus)
  #setwd(rLocs)
  dir.create("dsmartOuts/summaries/counts/",showWarnings = F)
  dir.create("dsmartOuts/summaries/probabilities/",showWarnings = F)
  dir.create("dsmartOuts/summaries/nProbable/",showWarnings = F)
  strc<- paste(getwd(),"/dsmartOuts/summaries/counts/",sep="")
  strp<- paste(getwd(),"/dsmartOuts/summaries/probabilities/",sep="")
  strn<- paste(getwd(),"/dsmartOuts/summaries/nProbable/",sep="")
  
  s1<- rLocs
  param<- nrow(lookup)
  
  
  #counts
  nme1<- paste(strc,"countOuts.tif" ,sep="") 
  f1<- function(x) {
    tabulate(x, nbins=param)}
  assign("param", param, envir=.GlobalEnv)
  counts<-clusterR(s1, calc, args=list(fun= f1), export = "param",filename=nme1,format="GTiff",overwrite=TRUE)
  setTxtProgressBar(pb, 25)
  
  
  #probabilities
  nme2<- paste(strp,"countOutsProbs.tif" ,sep="")
  param2<-nlayers(s1)
  f2<- function(x) (x/param2)
  assign("param2", param2, envir=.GlobalEnv)
  probs= clusterR(counts, calc,  args=list(fun=f2), export= "param2",filename=nme2,format="GTiff",overwrite=TRUE )
  setTxtProgressBar(pb, 50)
  
  if (sepP==TRUE) {s3<- stack()
                   for (np in 1:nlayers(probs)){
                     nme5<- paste(paste(strp, as.character(lookup[np,1]), sep=""), "_probs.tif", sep="")
                     names(probs[[np]])<- as.character(lookup[np,1])
                     s3<- stack(s3,probs[[np]])
                     writeRaster(probs[[np]],filename=nme5,format="GTiff",overwrite=TRUE)}}
  
  #Most probable
  nme3<- paste(strn,"nProbable.tif",sep="")
  f3<- function(x) order(x, decreasing=TRUE, na.last=TRUE)
  ordered.indices= clusterR(counts, calc,  args=list(fun=f3), filename=nme3,format="GTiff",overwrite=TRUE )
  s4<- stack()
  for (zz in 1:nprob){
    nme4<- paste(strn,paste(zz,"_probable.tif",sep=""),sep="")
    s4<- stack(s4,ordered.indices[[zz]])
    writeRaster(ordered.indices[[zz]],filename=nme4,format="GTiff",overwrite=TRUE)}
  setTxtProgressBar(pb, 75)
  #Most probable probabilities
  nme5<- paste(strp,"OrderedProbs.tif" ,sep="")
  f4<- function(x) sort(x, decreasing=TRUE, na.last=TRUE)
  ordered.probs= clusterR(probs, calc,  args=list(fun=f4), filename=nme5,format="GTiff",overwrite=TRUE )
  s5<- stack()
  for (zz in 1:nprob){
    nme6<- paste(strn,paste(zz,"_probableProbs.tif",sep=""),sep="")
    s5<- stack(s5,ordered.probs[[zz]])
    writeRaster(ordered.probs[[zz]],filename=nme6,format="GTiff",overwrite=T)}
  setTxtProgressBar(pb, 90)
  #Confusion Index
  nme7<- paste(strn,"confusionIndex.tif",sep="")
  f4<- function(x) (1-(x[[1]]-x[[2]]))
  confusInd<- clusterR(ordered.probs, fun=f4, filename=nme7, format="GTiff", overwrite=TRUE)
  setTxtProgressBar(pb, 100)
  endCluster()
  
  if (sepP==TRUE){retval<-list(s3,s4,s5,confusInd)} else{ retval<-list(s4)}
  close(pb)
  message(paste("DSMART outputs can be located at:",getwd(), sep=" "))
  return(retval)}

#END