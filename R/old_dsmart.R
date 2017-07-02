#Purpose: dsmart - disaggregation and harmonisation of soil map units through resampled classification trees
#Maintainer: Nathan Odgers (nathan.odgers@sydney.edu.au); Brendan Malone (brendan.malone@sydney.edu.au)
#Note: Algorithm described in [doi:10.1016/j.geoderma.2013.09.024]

#Basic description
# Randomly samples spatial polygons and builds a See5 classification tree
# with the sampling points.
#
# Args:
#   covariates: A RasterStack of covariate layers.
#   polygons: A SpatialPolygonsDataFrame containing the polygons to be
#      randomly sampled.
#   composition: A data frame containing the map unit composition for each
#      map unit polygon. First column is the polygon number (corresponding
#      to the first field of the SpatialPolygonsDataFrame attribute table);
#      second column is the map unit code; third column is the soil class
#      code; fourth column is the areal proportion of the map unit the soil
#      class is assumed to occupy.
#   n: The number of samples to draw from each polygon.
#   reals: Number of model resamplings to execute
#   cpus: Number of compute nodes to use

# Returns:
#   A number of items saved to file:
#          1. R date object that holds to model parameters from each fitted C5 model
#          2. Text files of the model summary output from each fitted C5 model
#          3. Rasters of each C5 model realisation. 
#          4. class and unique number lookup table
#   


dsmart<-function(covariates = NULL, polygons = NULL, composition = NULL, obsdat=NULL, n=NULL, reals = NULL, cpus=1){
  beginCluster(cpus)
  # Generate lookup table
  names(composition)<- c("poly", "mapunit", "soil_class", "proportion")
  lookup = as.data.frame(sort(unique(composition$soil_class)))
  lookup$code = seq(from=1, to=nrow(lookup), by=1)
  colnames(lookup) = c("name", "code")
  if(!is.null(obsdat)){names(obsdat)<- c("x", "y", "class")}
  
  #create output repositories
  model_lists<- vector("list", reals) #empty list 
  dir.create("dsmartOuts/",showWarnings = F)
  dir.create("dsmartOuts/rasters",showWarnings = F)
  dir.create("dsmartOuts/models",showWarnings = F)
  dir.create("dsmartOuts/summaries",showWarnings = F)
  strg<- paste(getwd(),"/dsmartOuts/rasters/",sep="")
  strm<- paste(getwd(),"/dsmartOuts/models/",sep="")
  strs<- paste(getwd(),"/dsmartOuts/summaries/",sep="")
  write.table(lookup, paste(strg,"classLookupTable.txt",sep=""),sep=",", col.names=T,row.names=F) 
  
  pb <- txtProgressBar(min=0, max=reals, style=3)
  for (j in 1:reals){
    # Empty data frame to store samples
    coordF<- matrix(NA, nrow=1000, ncol=3) #need to fix up the number of rows to better suit
    coordF<- data.frame(coordF)
    names(coordF)<- c("x", "y", "class")
    cf<- 1
    
    # take random samples within each polygon
    for(poly.id in polygons@data[,1]){
      #print(poly.id)
      # Subset a single polygon
      poly = subset(polygons, polygons@data[,1]==poly.id)
      #randomise points within polygon
      coordF[cf:(cf+(n-1)),1:2] = as.data.frame(spsample(poly, n , type="random", iter=10))
      
      # Allocate soil classes from within map unit
      poly.comp=subset(composition, composition$poly==poly.id)
      # Draw from Dirichlet distribution
      s=rdirichlet(1, poly.comp$proportion)
      
      # Weighted-random sample
      coordF$class[cf:(cf+(n-1))]=as.character(sample(poly.comp$soil_class, size=n, replace=TRUE, prob=s[1,]))
      cf<- cf+n}
    
      #spatial object
      locs<- as.data.frame(coordF[complete.cases(coordF),])
      locs<- rbind(locs,obsdat) # bind sampled data with observed data
      locs$num<- match(locs$class,lookup$name)
      coordinates(locs)<- ~ x + y  
      # Extract covariate values for the sampling locations
      values=extract(covariates,locs)
    
      #sample frame
      samples = cbind(as.data.frame(values),as.data.frame(locs)[,4])  
      names(samples)[ncol(samples)]<- "soil_class"
      samples$soil_class<- as.factor(samples$soil_class)
      samples<- samples[complete.cases(samples), ]
    
      #Fit model####
      res = C5.0(samples[,-ncol(samples)], y=samples$soil_class)
      model_lists[[j]]<- res
      #Capture output
      out<-capture.output(summary(res))
      f2<- paste(strm,paste("C5_model_", j,".txt",sep=""),sep="" )
      cat(out,file=f2,sep="\n",append=TRUE)
    
    nme<- paste(paste(paste(strg,"map",sep=""),"_",j,sep=""), ".tif", sep="")
    r1 <- clusterR(covariates, predict, args=list(res),filename=nme,format="GTiff",overwrite=T, datatype="INT2S")
    #plot(r1)
  setTxtProgressBar(pb, j)}
  
  #Save models to file
  save(model_lists, file = paste(paste(getwd(),"/dsmartOuts/",sep=""),"dsmartModels.RData", sep="") )
  endCluster()
  close(pb)
  message(paste(paste("DSMART outputs can be located at:",getwd(), sep=" "), "/dsmartOuts/",sep="") )}

#END