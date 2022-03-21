#' @title Evaluation of different spline variants.
#' 
#' @description The model selection method evaluates which spline models achieve the best quality among all tested metabolites. 
#' @param selectedLDO \code{\link{LDO}} containing all selected metabolites to be used for the model selection.
#' @param potentialBreaks Vector of all possible knots to be used for the spline modeling.
#' @param nknots Vector of number of spline knots to be used.  Therefore, 0 ~ no spline, 1 ~ spline with one knot, 2 ~ spline with two knots, etc.
#' @param splinetype spline type default is 'linear'. (Currently only linear is supported.)
#' @param qualityMeasure Vector of quality measures to be used. Possible options are 'AIC', 'BIC', and 'logLik'.
#' @return \code{\link{LDOmodelselection}} Object.
#' For each quality measure the model list contains a list of models for each spline tested. Additionally, the output contains a matrix of qualities for each Spline Component pair. And finally there is a list of breaks for each spline tested.
#' @importFrom methods new
#' @export 
#' @examples \dontrun{
#' } 
#'   data(LoBraExample)
#'   potentialBreaks <- c(8,12)
#'   selectedLDO <- selectComponents(ldo, components)
#'   ldoSelect<- lobraModelSelection(selectedLDO, potentialBreaks, nknots=c( 1, 2))
#'   length(ldoSelect@ldo@peaknames)
#'   
#'   
lobraModelSelection<-function(selectedLDO, potentialBreaks=c(), nknots=c(0, 1, 2), splinetype="linear", qualityMeasure=c("AIC", "BIC", "logLik")){
  modelList= list()
  quality= matrix(0,0,0)
  ctrl <- nlme::lmeControl(opt='optim');
  AICTable<-matrix(0,nrow=0,ncol=length(selectedLDO@peaknames));
  colnames(AICTable)<-selectedLDO@peaknames;
  BICTable<-AICTable;
  logLikTable<-AICTable;
  modelList<-list();
  breaks<-list();
  
  message("Start Model Selection ...");
  currentP <-selectedLDO@peaknames[1]
  if(0 %in% nknots || length(potentialBreaks)<1){
    AICS<-c();
    BICS<-c();
    logLikS<-c();
    mlist<-list();
    message("Evaluate Linear Models ...");
    for(currentP in selectedLDO@peaknames){
      mat<-selectedLDO@dataMatrices[[currentP]]
      model <- nlme::lme(fixed = value ~ class + time + class:time, random = ~ 1+time | id, data = mat, method = "ML", na.action = stats::na.exclude, control=ctrl);
      mlist[[currentP]]<-model;
      AICS<- c(AICS, summary(model)$AIC);
      BICS<-c(BICS, summary(model)$BIC);
      logLikS<-c(logLikS, summary(model)$logLik);
    }
    
    AICTable<-rbind(AICTable, AICS);
    BICTable<-rbind(BICTable, BICS);
    logLikTable<-rbind(logLikTable, logLikS);
    rownames(AICTable)<-rownames(BICTable)<-rownames(logLikTable)<-"No Spline"
    breaks[["No Spline"]] <-c();
    modelList[["No Spline"]]<-    mlist
    breaks[["No Spline"]] <-"No Spline";
  }

  
  ################################# Single knot spline models ##############################################################################################
  if(1 %in% nknots && length(potentialBreaks)>=1){
    message("Evaluate Single Knot Spline Models...");
    n<-rownames(AICTable);
    t<-8
    for(t in potentialBreaks){
      AICS<-c();
      BICS<-c();
      logLikS<-c();
      mlist<-list();
      
      for(currentP in selectedLDO@peaknames){
        mat<-selectedLDO@dataMatrices[[currentP]]
        inter.knot <- c(t)
        names(inter.knot) <- c("TimeV")
        E2<- outer(mat$time, inter.knot,"-")
        ls.mat <- E2*(E2>0)
        mat$timeV<- ls.mat[,1]
        model.spline <- nlme::lme(fixed = value ~ class + time + timeV + class:time + class:timeV , random = ~ 1+time | id, data = mat, method = "ML", na.action = stats::na.exclude, control=ctrl);
        mlist[[currentP]]<-model.spline;
        AICS<- c(AICS, summary(model.spline)$AIC);
        BICS<-c(BICS, summary(model.spline)$BIC);
        logLikS<-c(logLikS, summary(model.spline)$logLik);
      }
      
      AICTable<-rbind(AICTable, AICS);
      BICTable<-rbind(BICTable, BICS);
      logLikTable<-rbind(logLikTable, logLikS);
      breaks[[paste("T-",t, sep="")]] <-c(t);
      modelList[[paste("T-",t, sep="")]]<-    mlist
    }
    n<- c(n, paste("T-",potentialBreaks, sep=""))
    rownames(AICTable)<-rownames(BICTable)<-rownames(logLikTable)<-n
  }
  
  
  ################################# Two+ knots spline #################################################################################################
  nknots<- nknots[nknots>1]
  if(length(nknots)>0 && length(potentialBreaks)>= min(nknots)){
    
    pS<-powerSet(1:length(potentialBreaks))
    pSlengths<-lapply(pS, length)
    k<-2
    for(k in nknots){
      n<-rownames(AICTable);
      message(paste("Evaluate ", k , "-knot Spline Models..."));
      idsSets<- pS[pSlengths==k];
      #ids<-idsSets[[4]]
      for(ids in idsSets){
        inter.knot <- potentialBreaks[ids]
        names(inter.knot) <- paste("time",inter.knot, sep = "")
        tcomb<-paste("time",inter.knot, collapse = "-", sep = "")
        breaks[[tcomb]] <-c(inter.knot);

        AICS<-c();
        BICS<-c();
        logLikS<-c();
        mlist<-list();
        
        
        currentP <- selectedLDO@peaknames[[4]]
        for(currentP in selectedLDO@peaknames){
          mat<-selectedLDO@dataMatrices[[currentP]]
          E2<- outer(mat$time, inter.knot,"-")
          ls.mat <- E2*(E2>0)
          colnames(ls.mat)<- names(inter.knot)
          
          for(t in names(inter.knot)){
            mat[t]<- as.vector(ls.mat[,t])
          }
          
          fixedf <- stats::as.formula(paste(c("value ~ class + time + class:time", names(inter.knot) , paste("class:", names(inter.knot), sep="")) ,collapse =" + "))
          model.spline <- nlme::lme(fixed = fixedf , random = ~ 1+time | id, data = mat, method = "ML", na.action = stats::na.exclude, control=ctrl);
          mlist[[currentP]]<-model.spline;
          
          AICS<- c(AICS, summary(model.spline)$AIC);
          BICS<-c(BICS, summary(model.spline)$BIC);
          logLikS<-c(logLikS, summary(model.spline)$logLik);
        }
        
        AICTable<-rbind(AICTable, AICS);
        BICTable<-rbind(BICTable, BICS);
        logLikTable<-rbind(logLikTable, logLikS);
        modelList[[tcomb]]<-    mlist
        n<- c(n, tcomb)
      }
      rownames(AICTable)<-rownames(BICTable)<-rownames(logLikTable)<-n
    }
  }
  quality=list()
  if("AIC" %in% qualityMeasure){
    quality[["AIC"]]<- AICTable;
  }
  if("BIC" %in% qualityMeasure){
    quality[["BIC"]]<- BICTable;
  }
  if("logLik" %in% qualityMeasure){
    quality[["logLik"]]<- logLikTable;
  }
  
  lobraModelSelectionObject<- methods::new("LDOmodelselection", ldo=selectedLDO, 
                             potentialBreaks=potentialBreaks, 
                             splinetype=splinetype, 
                             qualityMeasure=qualityMeasure, 
                             modelList=modelList,
                             quality=quality,
                             breaks=breaks)
  return(lobraModelSelectionObject)
}


#' @title Plotting results of Model Evaluation and Selection.
#' 
#' @description Plotting the results of Model Evaluation and Selection. The plot shows a vertical boxplot for each spline tested starting with the best average fit according to the selected quality measure. The label of each spline can be found on the left, the median quality measure on the right. The x-axis denotes the selected quality measure. 
#' @param lobraModelSelectionObject Object of type LDOmodelselection that was created during the model evaluation. @seealso 'lobraModelSelection'
#' @param qualityMeasure List of quality measures to be visualized. 
#' @param title Title of the plot.
#' @return No return value
#' @export 
#' @examples \dontrun{
#' } 
#'
#' wd <- tempdir()
#' data(LoBraExample)
#' selectedLDO <- selectComponents(ldo, components)
#' ldoSelect<- lobraModelSelection(selectedLDO, potentialBreaks=c(8, 12), nknots=c(1, 2))
#' 
#' filename<- file.path(wd, "evaluateBestSplineAIC.pdf") ;
#' grDevices::pdf(filename, width=16, height=8);
#'   plotmodelSelectionEvaluation(ldoSelect, "AIC", "Best Spline Models");
#' grDevices::dev.off();
#'   
#' qualityMeasure=c("AIC", "BIC", "logLik")
#' filename<- file.path(wd, "evaluateBestSplineAllMeasures.pdf") ;
#' grDevices::pdf(filename, width=16, height=8);
#' oldpar <- par("mfrow")
#' par(mfrow=c(3,1))
#'   plotmodelSelectionEvaluation(ldoSelect, qualityMeasure);
#' par(mfrow = oldpar)
#' grDevices::dev.off();
#' 
plotmodelSelectionEvaluation<-function(lobraModelSelectionObject, qualityMeasure, title=NULL){
  q<-qualityMeasure[3]
  for(q in qualityMeasure){
    if(is.null(title)){
      mtitle=paste("", q, " Spline Model Comparison");
    }else{
      mtitle=title;
    }
    data<-lobraModelSelectionObject@quality[[q]]
    oldpar<-graphics::par()$mar;
    graphics::par(mar = c(5, 10, 4, 6)+0.1)
    if(q =="logLik"){
      vals<-sort(apply(data, 1, FUN=stats::median), decreasing = TRUE)
    }else{
      vals<-sort(apply(data, 1, FUN=stats::median))
    }
    
    graphics::boxplot(t(data)[,rev(names(vals))], horizontal=TRUE, col="lightblue", xlab=q, las=2,  main=mtitle)
    graphics::axis(4,at=1:length(rownames(data)),adj=1,labels=rev(round(vals)), las=2)
  }
  graphics::par(mar =oldpar);
}

#' @title Extract the optimal spline model parameters from the ModelSelection Object.
#' 
#' @description The method calculates which spline model and parameters worked best with respect to the median of the specified quality measure. The median is calculated among all component models.
#' @param lobraModelSelectionObject LDOmodelselection created by the 'lobraModelSelection' function. It stores all evaluated Spline models to chose from.
#' @param qualityMeasure Quality measure to be used to select the optimal spline.
#' @param summeryfun Define the Summery function to be used. Default value is set to stats::median. Other possible functions would be mean, for instance.
#' @return The function returns a 'lobraModelSelectionObject' that contains the optimal model according to the specified quality measure. @seealso plot.modelSelectionEvaluation
#' @importFrom methods new
#' @export 
#' @examples \dontrun{
#' } 
#'   data(LoBraExample)
#'   selectedLDO <- selectComponents(ldo, components)
#'   potentialBreaks=c(8, 12)
#'   nknots=c(1, 2)
#'   qualityMeasure=c("AIC", "BIC")
#'   ldoSelect<- lobraModelSelection(selectedLDO, potentialBreaks, nknots, qualityMeasure)
#'   
#'   optimalAIC<-getOptimalSpline(ldoSelect, qualityMeasure="AIC", summeryfun=stats::median)
#'   message(optimalAIC@breaks);
#'   
#'   optimalBIC<-getOptimalSpline(ldoSelect, qualityMeasure="BIC", summeryfun=base::mean)
#'   hist(unlist(optimalBIC@quality));
#' 
getOptimalSpline<-function(lobraModelSelectionObject, qualityMeasure="AIC", summeryfun=stats::median){
  meds<-apply(lobraModelSelectionObject@quality[[qualityMeasure]], 1, summeryfun)
  which(meds==min(meds))
  if(qualityMeasure =="logLik"){
    optimalID<-names(which(meds==max(meds)))
  }else{
    optimalID<-names(which(meds==min(meds)))
  }
  optimalModel<- lobraModelSelectionObject@modelList[[optimalID]]
  length(optimalModel)
  optimalSpline<- lobraModelSelectionObject@breaks[[optimalID]]
  
  
  optimalModelObject<- methods::new("LDOmodelselection", ldo=lobraModelSelectionObject@ldo,
                             potentialBreaks=optimalSpline,
                             splinetype=lobraModelSelectionObject@splinetype,
                             qualityMeasure=qualityMeasure,
                             modelList=optimalModel,
                             quality=list(lobraModelSelectionObject@quality[[qualityMeasure]][optimalID,]),
                             breaks=list(optimalID=optimalSpline))
  
  
  return(optimalModelObject)
}


#' @title Creating the power set of a set.
#' @param set Set of numbers of potential spline break points.
#' @return Returns power set of the given set.
powerSet <- function(set) { 
  n <- length(set);
  masks <- 2^(1:n-1);
  powerset<-lapply( 1:2^n-1, function(u) set[ bitwAnd(u, masks) != 0 ] );
  return(powerset);
}

