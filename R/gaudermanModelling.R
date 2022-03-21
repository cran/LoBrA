#' @title Create the Gouderman Data Arrangement.
#' 
#' @description  Using the Gouderman methodology to create the Gouderman-Data Arrangement.
#' @param selectedLDO Longitudinal Data Object, containing all selected metabolites to be used for the final Gouderman model.
#' @param breaks break points for the spline model
#' @param center Time point that corresponds to the center time t0. The algorithm will test whether there is a significant difference between the groups at this point.
#' @param timeperiod If the user defines the time period or segment, in the spline to be tested. Note, a 3 break point spline has 4 segments.
#' @param range If the user defines a range, the algorithm will test whether there is a significant difference between the groups in this range.
#' @return The function returns a 'GaudermanLDO' object. For more information @seealso 'GaudermanLDO' .
#' @importFrom methods new
#' @export 
#' @examples \dontrun{
#' } 
#'   data(LoBraExample)
#'   selectedLDO <- selectComponents(ldo, components)
#'   breaks<- c(8, 12)
#'   center<- 12
#'   timeperiod <- 2;
#'   gaudermanLDOexample <- createGoudermanData(selectedLDO, breaks, center, timeperiod)
#'   
createGoudermanData<-function(selectedLDO, breaks, center, timeperiod=NA, range=NA){

  ## Checking time period and range values.
  times<-selectedLDO@times
  bs<-c(min(times), breaks, max(times))
  if(!is.na(timeperiod)){
    range<-c(bs[timeperiod], bs[timeperiod+1])
  }else if( length(!is.na(range))>1 ){
    if(length(range)!=2){
      stop("Range must be a vector with two entries!")
    }
    range<- sort(range)
    i<- which(bs==range[1])
    j<- which(bs==range[2])
    # i<- which(bs==4)
    if(length(i)==1 & length(j)==1){
      timeperiod<-i;
    }else{
      stop("Each value in range must included once and only once in the breaks vector!")
    }
  }else{
    stop("Function requires either time period or range to be defined!")
  }
  if(!(center>=range[1] & center<=range[2])){
    stop("The center variable must be within the defined range!")
  }
  
  ## Read required parameters.
  sampleIds<-selectedLDO@ids
  peaknames<-selectedLDO@peaknames
  classes<- selectedLDO@labels
  
  ##Caluculate 'Gauderman' range and the new beginning and end of this range. If one of them is equal to the center, it will be set to 0, the remaining on will be either -1 or 1.
  gaudermanRange<- abs(range[1] - range[2]);

  ## Determine the new breaks for the spline.
  k<-c()
  for(i in 1:length(breaks)){
    b<-breaks[i]
    b<- (b-center)/gaudermanRange
    k<-c(k, b)
  }
  names(k)<-1:length(k)
  
  # Create new Gouderman data frames
  dataFrames<-list();
  newtimes<-c()
  # selectedLDO@times
  p<- peaknames[1]
  for(p in peaknames){
    peakmatrix<-selectedLDO@dataMatrices[[p]]
    
    myDataFrame<- getGeneralizedGaudermanDataFrame(peakmatrix, sampleIds, classes, center, timeperiod, gaudermanRange, k)
    dataFrames[[p]]<-myDataFrame;
    newtimes<-c(  newtimes, myDataFrame[,"time"])
    
  }
  newTimeVars<-setdiff(colnames(myDataFrame), colnames(peakmatrix));
  newtimes<- as.matrix(unique(newtimes));
  myGaudermanLDO<- methods::new("GaudermanLDO",
                       name=paste0("Gauderman-",selectedLDO@name), 
                       dataFrames=dataFrames, 
                       peaknames=peaknames, 
                       k= k,
                       times=newtimes, 
                       newTimeVars=newTimeVars,
                       ids=sampleIds, 
                       labels=classes)
  return(myGaudermanLDO);
}




#' @title Create Peak Matrices for Generalized 'Gauderman' linear mixed effect regression (LMER) Model with parameterized Times
#' 
#' @param peakmatrix Peak matrix to be converted.
#' @param sampleIds Ids of samples in the matrix
#' @param classes Classes of samples
#' @param center Time point that corresponds to the center time t0. The algorithm will test whether there is a significant difference between the groups at this point.
#' @param timeperiod defines the time period or segment, in the spline to be tested. Note, a 3 break point spline has 4 segments.
#' @param gaudermanRange range to be tested for a significant difference between the groups.
#' @param k break points for the generalized 'Gauderman' spline model.
#' @return Return the new peak matrix data frame for this peak.
#' 
getGeneralizedGaudermanDataFrame=function(peakmatrix, sampleIds, classes, center, timeperiod, gaudermanRange, k){
  ## Scale the old times to the new 'Gauderman' range, standardizing the range between rangeStart and rangeEnd to 1.
  oldMins<- as.numeric(peakmatrix[,'time']);
  newMins = (oldMins - center)/gaudermanRange;
  
  ## Define the new Gouderman time variables for (1) the time points before k1 (2) the time points between k1 and k2 (k1>=t<=k2) which is the time of interest and (3) the time points larger than k2.
  numberS<-length(k)+1
  timematrix<-matrix(rep(newMins, numberS), ncol = numberS, nrow = length(newMins))
  k0<- c(min(newMins), k, max(newMins))
  names(k0)<-1:length(k0)
  i<-4
  for(i in 1:numberS){
    t<- timematrix[,i]
    timematrix[t<=k0[i],i]<- k0[i]
    timematrix[t>k0[i+1],i]<- k0[i+1]
  }

  ## Create a data.frame with all previous fields.
  peakFrame <- data.frame(class=factor(peakmatrix[,"class"], levels=classes), 
                          id = factor(peakmatrix[,"id"], levels=sampleIds), 
                          time=as.numeric(newMins),
                          value=as.numeric(peakmatrix[,"value"]))

  ## Remove k[i] for each column of the time-matrix that is not equal to the specified time period according to the generalized 'Gauderman' algorithm.
  ## Add each column to the data.frame.
  i=1;
  j=1;
  while(i <=numberS){
    t<-as.numeric(timematrix[,i])
    n<-paste0("timep",i)
    if(i==timeperiod){
      peakFrame[,n]<-t
      i=i+1
      next;
    }
    peakFrame[,n] <- t-as.numeric(k[j])
    j=j+1;
    i=i+1;
  }
  
  ## Return the new data frame for this peak.
  return(peakFrame);
}




#' @title Fitting the Gouderman LME Model with using Gouderman-Data Arrangement.
#' 
#' @description Uses the linear mixed effects modeling to build the final 'Gauderman' model. The 'Gauderman' modification enables the exact calculation of the significance of a specified section of the spline model.
#' @param mygaudermanLDO GaudermanLDO data object, created by the generalized 'Gauderman' algorithm (GGA).
#' @param correctionMethod correction for p-values. Possible methods: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'
#' @return 'GaudermanModelEvaluation' Results of the evaluation of the Fitted linear mixed effect models for the defined time periods.
#' @importFrom methods new
#' @export 
#' @examples 
#'   
#'   data(LoBraExample)
#'   selectedLDO <- selectComponents(ldo, components)
#'   gaudermanLDOexample <- createGoudermanData(selectedLDO, breaks=c(8, 12), center=12, timeperiod=2)
#'   evalResult<- modelGoudermanLongitudinal(gaudermanLDOexample)
#'   message(evalResult@correctedpvalues<0.005)
#'   
modelGoudermanLongitudinal<-function(mygaudermanLDO, correctionMethod="bonferroni"){
  
  ## Get required parameter from 'mygaudermanLDO'. 
  peaknames<- mygaudermanLDO@peaknames
  inter.knot2<- mygaudermanLDO@newTimeVars
  k<- mygaudermanLDO@k;
  classes<- mygaudermanLDO@labels
  
  
  ## Arrays and Lists to gather the models, p-values that go into the 'GaudermanModelEvaluation'.
  ctrl <- nlme::lmeControl(opt='optim');
  pValues<- c();
  modelparameter<- c();
  modellist<- list();
  allIntercepts <- c()
  allSlopes<-c()
  
  ## Run over all selected peaks and fit the specified linear mixed effect model.
  for(currentP in peaknames){
    ## Get the generalized gauderman data frame.
    mat<-mygaudermanLDO@dataFrames[[currentP]]

    ## Fit the linear mixed effect model.
    fixedf <- stats::as.formula(paste(c("value ~ class ", inter.knot2 , paste("class:", inter.knot2, sep="")) ,collapse =" + "))
    model.spline <- nlme::lme(fixed = fixedf , random = ~ 1+time | id, data = mat, method = "ML", na.action = stats::na.exclude, control=ctrl);
    modellist[[currentP]]<-model.spline;
    
    ## Aquiring p-values for the Class intercept comparison and the slops for all modeled time periods.
    classname<-paste0("class",classes[2])
    pvaluesOfinterest<- c( classname, paste0(classname, ":", inter.knot2))
    pv<- summary(model.spline)$tTable[pvaluesOfinterest,"p-value"];
    pValues<- rbind(pValues, pv);
    
    ## Aquiring the model parameter.
    modelparameter<- rbind(modelparameter, nlme::fixef(model.spline));

  }
  
  correctedpValues<- c();
  mycol<-colnames(pValues)[1]
  for(mycol in colnames(pValues)){
    corrected<- stats::p.adjust(pValues[,mycol], method = correctionMethod);
    correctedpValues<- cbind(correctedpValues,corrected);
  }
  rownames(pValues) <- peaknames;
  rownames(correctedpValues) <- peaknames;
  colnames(correctedpValues) <- colnames(pValues);
  rownames(modelparameter) <- peaknames;
  
  
  
  #Recover labels
  b<-unique(paste(mygaudermanLDO@dataFrames[[1]][,'class'],mygaudermanLDO@dataFrames[[1]][,'id']))
  labels<- as.factor(sapply(b, function(x){unlist(strsplit(x," "))[1]}))
  names(labels)<- sapply(b, function(x){unlist(strsplit(x," "))[2]})

  
  gModelEvaluation<- methods::new("GaudermanModelEvaluation",
                       name=paste0("Evaluation-",mygaudermanLDO@name), 
                       gaudermanLDO=mygaudermanLDO, 
                       models=modellist,
                       labels=labels,
                       pvalues=pValues,
                       correctedpvalues=correctedpValues,
                       modelparameter=modelparameter)
  
  return(gModelEvaluation) 
}




#' @title Plotting the 'Gouderman' LME Model and Results.
#' 
#' @param evaluationresult 'GaudermanModelEvaluation' data object, created by the modelGoudermanLongitudinal function.
#' @param main title of the plot
#' @param ylab y axis label
#' @param xlab x axis label
#' @param peaknames selection of peaks to be plotted
#' @return No return value
#' @export 
#' @examples 
#'   
#'   wd <- tempdir()
#'   data(LoBraExample)
#'   selectedLDO <- selectComponents(ldo, components)
#'   gaudermanLDOexample <- createGoudermanData(selectedLDO, breaks=c(8, 12), center=12, timeperiod=2)
#'   evalResult<- modelGoudermanLongitudinal(gaudermanLDOexample)

#'   # Plot all peaks
#'   filename<- file.path(wd, "finalModelEvaluation.pdf") ;
#'   oldpar <- par("mfrow")
#'   grDevices::pdf(filename, width=16, height=8);
#'     graphics::par(mfrow=c(1,1));
#'     plotGoudermanLongitudinalResults(evalResult);
#'   par(mfrow = oldpar)
#'   grDevices::dev.off();
#'   
#'   #Plot a selection of Peaks
#'   peaknames<- evalResult@gaudermanLDO@peaknames;
#'   filename<- file.path(wd, "finalModelEvaluation-components.pdf") ;
#'   oldpar <- par("mfrow")
#'   grDevices::pdf(filename, width=20, height=8);
#'     graphics::par(mfrow=c(2,5));
#'     plotGoudermanLongitudinalResults(evalResult, main="", peaknames=peaknames);
#'   par(mfrow = oldpar)
#'   grDevices::dev.off();
#'  
plotGoudermanLongitudinalResults<-function(evaluationresult, main="Mixed Effect Spline Model Evaluation", ylab = "Value", xlab="Time", peaknames=NULL){
  
  labels<- evaluationresult@labels
  breaks<- evaluationresult@gaudermanLDO@k

  # Set color ranges for each class:
  ul<- unique(labels)
  mycol<-rep("", length(labels))
  maincol<-c()
  for(l in ul){
    mycol[labels==l]<- getColor(l, sum(labels==l))
    maincol<- c(maincol, getColor(l, 1))
  }
  names(mycol)<-names(labels)
  names(maincol)<-ul
  

  
  #plot data
  if(is.null(peaknames))
    peaknames<- evaluationresult@gaudermanLDO@peaknames
  
  for( p in peaknames ){
    data <- evaluationresult@gaudermanLDO@dataFrames[[p]]
    tempmodel<- evaluationresult@models[[p]]
    colores<- mycol[as.character(data[,'id'])]
    maxVal<-max(data[,'value'])
    plotGaudermanModel(data, labels, ul,  tempmodel, colores, maincol, breaks, main=paste(main," '", p, "'", sep = ""), ylab , xlab)

    # Add Pvalues for intercept and each segment
    minTime<-min(data[,'time'])
    maxTime<-max(data[,'time'])
    pvalues<- evaluationresult@correctedpvalues[p,]
    interceptPV<- pvalues[1]
    pvalues<-pvalues[-1]
    ts<-(c(minTime, breaks)+c(breaks, maxTime))/2
    pvcolors<-c("gray", "darkred")
    names(pvcolors)<- c(FALSE,TRUE)
    graphics::text(c(0), c(maxVal)-1.5, labels = c(paste(ifelse(interceptPV<0.005, "*", ""), "(", round(interceptPV,3), ")", sep = "" ) ), col = pvcolors[as.character(interceptPV<0.005)])
    graphics::text(ts, c(maxVal)-1.5, labels = c(paste( ifelse(pvalues<0.005, "* ", ""), "(", round(pvalues,3), ")", sep="")), col=pvcolors[as.character(pvalues<0.005)])
  }
    
}


#' @title Plotting helper function to plot a single generalized gouderman Model
#' 
#' @param data      data matrix used to fit the model
#' @param labels    class labels for all samples
#' @param ul        unique class labels
#' @param tempmodel model to be plotted
#' @param colores   predefined colors for the single samples
#' @param maincol   predefined colors for the fitted spline
#' @param breaks    break points of the spline to be plotted
#' @param main      main title of the plot
#' @param ylab      y label of the plot
#' @param xlab      x label of the plot
#' 
plotGaudermanModel=function(data, labels, ul,  tempmodel, colores, maincol, breaks, main, ylab , xlab){
  #Plot Original Data
  time<- data[,'time']
  value<- data[,'value']
  minVal<-min(value)
  maxVal<-max(value)
  graphics::plot(time, value, col=colores, pch=8, ylab=ylab, xlab = xlab, main=main)
  for (b in breaks) {
    graphics::lines(c(b,b),  c(minVal, maxVal), lty=3, col='darkgrey', lwd=3)
  }
  graphics::text(c(0), c(maxVal)-0.5, labels = c("Intercept"), col = "darkgrey")
  
  #Plot fitted Spline
  fittedvalues<- c(stats::fitted(tempmodel, 0))
  classes <- labels[names(fittedvalues)]
  l<-ul[1]
  for(l in ul){
    idx<- classes==l
    idx2 <- data[,"class"]==l
    utimes<- time[idx]
    names(utimes)<-1:length(utimes)
    utimes <- sort(utimes);
    ufitted<- fittedvalues[idx]
    ufitted<- ufitted[as.numeric(names(utimes))]
    graphics::lines(utimes, ufitted, col=maincol[l] , type = 'l', lwd=5)
  }
  
  #Add model parameter
  # modelfixed= tempmodel$coefficients$fixed
  # interceptC1= modelfixed["(Intercept)"]
  # interceptC2= interceptC1 + modelfixed["class2"]
  # text(c(0,0), c(interceptC1, interceptC2)+0.5, labels = c(round(interceptC1,2), round(interceptC2,2) ))
}

