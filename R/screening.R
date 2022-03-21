#' @title Screening of background or confounding components
#' 
#' @description Background noise signals originating from experimental settings or random events can hugely influence the signal pattern of the breath. Background data enables the detailed evaluation and differentiation of the compounds originating primarily from the background or confounding factors as compared to those from the sample itself. The method assumes that all compounds of interest show a larger variation in the sample as compared to the background noise.
#' @param ldo Longitudinal Data Object
#' @param method list of tests to perform, standard values: 'bf', 'levene' or 'bartlett'). 'bf' relates to "Brown-Forsythe" Levene-type procedure, 'levene' uses classical "Levene's" procedure and 'bartlett' applies Bartlett's test.
#' @param alpha A numeric value to defining the cutoff to select peaks.
#' @param criteria indicators which criteria to use for screening decision.
#' @return Returns an object of type 'LDOscreening' containing the original 'ldo' object and the results of the screening. The variable 'selectedPeaks' contains a matrix including the results (TRUE = Significant, FALSE = not Significant) of the specified tests ('bf', 'levene', 'bartlett').
#' @importFrom methods new
#' @export 
#' @examples \dontrun{
#' } 
#' 
#'   data(LoBraExample)
#'   method= c('bf', 'levene', 'bartlett')
#'   alpha =0.05
#'   criteria=c(1,1)
#'   ldos<-screening(ldo, method, alpha, criteria)
#'   components <- ldos@selectedPeaks[,"levene"]
#'   components <- names(components)[components]
#'   selectedLDO <- selectComponents(ldo, components)
#' 
#' 
screening<-function(ldo, method=c("bf", "levene", "bartlett"), alpha =0.05, criteria=c(1,1)){
  allalphapvals<-c();
  allrespvals<-c();
  
  sAlphasVar<- c();
  bgAlphasVar<- c();
  sResVar<- c();
  bgResVar<- c();
  experimentResiduals<-list();
  experimentIntercept <-list();
  
  currentP<- ldo@peaknames[71]
  oldpar <- graphics::par(no.readonly = TRUE);
  on.exit(graphics::par(oldpar));
  
  for(currentP in ldo@peaknames){
    graphics::par(mfrow=c(1,2));
    sPeaks<- c()
    ventilatorPeaks <- c()
    rmeans<- c()
    rvars<-c()
    vmeans<- c()
    vvars<-c()
    
    smodel<- stats::lm(value~id, ldo@dataMatrices[[currentP]])
    sAlphas<- smodel$coefficients[-1]
    sRes<- smodel$residuals
    sResVar<- c(sResVar, stats::var(sAlphas));
    sAlphasVar<- c(sAlphasVar, stats::var(sRes));
    
    bmodel<- stats::lm(value~id, ldo@backgroundMatrices[[currentP]])
    bgAlphas<- bmodel$coefficients[-1]
    bgRes<- bmodel$residuals
    bgAlphasVar<- c(bgAlphasVar, stats::var(bgAlphas));
    bgResVar<- c(bgResVar, stats::var(bgRes));
    
    alpy<- c(sAlphas, bgAlphas)
    alpgroup<-factor(c(rep("S", length(sAlphas)), rep("B", length(bgAlphas))))
    resy<- c(sRes, bgRes)
    resgroup<-factor(c(rep("S", length(sRes)), rep("B", length(bgRes))))
    alphapvals<-c();
    respvals<-c();
    
    for(m in method){
      alphapvals<-c( alphapvals, getPvalue(alpy, alpgroup, m));
      respvals<-c( respvals, getPvalue(resy, resgroup, m));
    }
    
    allalphapvals<-rbind(allalphapvals, alphapvals);
    allrespvals<-rbind(allrespvals,respvals);
    experimentIntercept[[currentP]] <- list(sAlphas, bgAlphas);
    experimentResiduals[[currentP]] <- list(sRes, bgRes);
  }
  rownames(allalphapvals)<-ldo@peaknames;
  rownames(allrespvals)<-ldo@peaknames;
  colnames(allalphapvals)<-method;
  colnames(allrespvals)<-method;
  
  
  selectedPeaks<-c(); 
  for(m in method){
    pv1<-allalphapvals[,m];
    pv2<-allrespvals[,m];
    cbind(pv1,pv2);
    tested <- cbind(pv1<alpha, pv2<alpha, bgAlphasVar<sAlphasVar, bgResVar<sResVar);
    if(sum(criteria)>1){
      accepted <- (tested[,1] & tested[,3]) | (tested[,2] & tested[,4]);
    }else if(criteria[1]==1){
      accepted <- (tested[,1]) | (tested[,2]);
    }else if(criteria[2]==1){
      accepted <- (tested[,3]) | (tested[,4]);
    }
    selectedPeaks<-cbind(selectedPeaks, accepted); 
    # Accepted <- rbind(Accepted, apply(Accepted, 2, FUN=sum));
    
  }
  colnames(selectedPeaks)<-method;
    
  
  ldoscreen<-methods::new("LDOscreening",ldo=ldo, experimentIntercept=experimentIntercept, experimentResiduals=experimentResiduals,  interceptPvalues=allalphapvals, residualPvalues=allrespvals, selectedPeaks=selectedPeaks)
  return(ldoscreen);
}




#' @title Plotting the screening results. 
#' @description For each peak two box plots are created. The first plot shows a boxplot of the Sample Intercept Comparison of the sample and the background, and the corresponding p-values. The second plot shows a boxplot of the Residual Comparison of the sample and the background, and the corresponding p-values.
#' @param ldoscreen LDO screening result
#' @param plotAll Select all components to be plotted. Default plots only the selected peaks using the correction method.
#' @param correctionmethod Version of correction method to be used to select the peaks. Valid values are 'bf', 'levene', and 'bartlett'.
#' @param decs decimal numbers of p-values to be plotted.
#' @param ask logical. Modifies the graphical parameter \code{ask} in \code{par} (If TRUE (and the R session is interactive) the user is asked for input, before a new figure is drawn. As this applies to the device, it also affects output by packages grid and lattice. It can be set even on non-screen devices but may have no effect there.) 
#' @param peaknames Defining a list of peaks to be plotted. By default all peaks will be plotted.
#' @return No return value
#' @export 
#' @examples \dontrun{
#' } 
#' 
#'   wd <- tempdir()
#'   data(LoBraExample)
#'   ldos<-screening(ldo, method= c('levene'), alpha =0.05, criteria=c(1,1))
#'   filename<- file.path(wd, "screeningresults.pdf") 
#'   grDevices::pdf(filename, width=16, height=8)
#'   plotLDOScreening(ldos)
#'   grDevices::dev.off();
#' 
plotLDOScreening=function(ldoscreen, plotAll=FALSE, correctionmethod="levene", decs=3, ask=FALSE, peaknames=rownames(ldoscreen@selectedPeaks)){
  if(!plotAll){
    peaknames<- names(which(ldoscreen@selectedPeaks[peaknames,correctionmethod] ))
  }
  
  p<-peaknames[30]
  oldpar <- graphics::par(no.readonly = TRUE);
  on.exit(graphics::par(oldpar));
  for(p in peaknames){
    graphics::par(mfrow=c(1,2), ask=ask)
    intensities <-ldoscreen@experimentIntercept[[p]]
    pv<-ldoscreen@interceptPvalues[p,]
    names(intensities)<- c("Sample", "Background")
    graphics::boxplot(intensities, main = paste("Sample Intersept Comparison", p), col="lightblue")
    graphics::legend("topright", c("P-values:", paste("Brown-Forsythe:", round(pv["bf"], decs)), paste("Levene Test" ,round(pv["levene"], decs)), paste("Bartlett's Test" ,round(pv["bartlett"], decs))))
    
    intensities <-ldoscreen@experimentResiduals[[p]]
    pv<-ldoscreen@residualPvalues[p,]
    names(intensities)<- c("Rat", "Ventilator")
    graphics::boxplot(intensities, main = paste("Rat Residual Comparison", p), col="lightblue")
    graphics::legend("topright", c("P-values:", paste("Brown-Forsythe:", round(pv["bf"], decs)), paste("Levene Test" ,round(pv["levene"], decs)), paste("Bartlett's Test" ,round(pv["bartlett"], decs))))
    
  }
}


#' @title Testing differences of groups with respect to a specific value and test. 
#' @param y Values to be tested
#' @param group corresponding groups whose difference we want to test
#' @param test specific test to be used. Can be each of the following 'bf', 'levene' or 'bartlett'.
getPvalue=function(y, group, test){
  if(test=="bf"){
    pv<- lawstat::levene.test(y, group)$p.value
  }else if(test=="levene"){
    pv<- lawstat::levene.test(y, group, location="mean")$p.value
  }else if(test=="bartlett"){
    pv<- stats::bartlett.test(x=y, g=group)$p.value
  }else{
    stop("Error: test needs to be either of bf, levene or bartlett"); 
  }
  return(pv)
}
