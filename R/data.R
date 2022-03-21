#' 'LoBrA'  Example Data Set
#'
#' 'LoBrA'  example data set created by the function 'createExampleData'. #' It consist of a single matrix for all experiments, time points, types (background, experiment), class and the intensity values of all components created.
#' The artificial data consist of 20 experiments and 100 components with 18 measurements (3 background, 15 sample). The 10 experiments are each associated to on of 2 groups (ONE and TWO). The components comprise 70 noise components and 30 components that randomly vary in their trajectories in one of three segments. Random noise is added to all intercepts, propagated and added to each time point for all samples and components separately.
#'
#' @name longDataExample
#' @docType data
#' @author Anne-Christin Hauschild \email{hauschild@@uni-marburg.de}
#' @keywords Breathomics Metabolomics Expression Longitudinal-Data Confounders
#' @format A matrix representing 20 experiments. It contains values for 100 variables at 18 time points for each experiment.
#' \describe{
#'   \item{id}{Experiment identifier}
#'   \item{time}{Time Point of Measurement}
#'   \item{type}{Type of Measurement (e.g. Background, or Sample measurement for each experiment)}
#'   \item{class}{Class or Group id of the sample/ experiment}
#'   \item{bgcomponent-x}{70 random variables that represent the background noise of the experiments}
#'   \item{components-x-x}{30 components that randomly vary in their trajectories in one of three time periods, (1:4-8, 2:9-13, 3:14-18). } 
#'   ...
#' }
"longDataExample"

#' 'LoBrA'  Data Object (LDO) for Example data set
#'
#' 'LoBrA'  example LDO created by the function 'createExampleData' and converted to an LDO by 'as.LOBdataset' function.
#' It consist of a single matrix for all experiments, time points, types (background, experiment), class and the intensity values of all components created.
#' The artificial data consist of 20 experiments and 100 components with 18 measurements (3 background, 15 sample). The 10 experiments are each associated to on of 2 groups (ONE and TWO). The components comprise 70 noise components and 30 components that randomly vary in their trajectories in one of three segments. Random noise is added to all intercepts, propagated and added to each time point for all samples and components separately.
#'
#' @name ldo
#' @docType data
#' @author Anne-Christin Hauschild \email{hauschild@@uni-marburg.de}
#' @keywords Breathomics Metabolomics Expression Longitudinal-Data Confounders
#' @format A matrix representing 20 experiments. It contains values for 100 variables at 18 time points for each experiment.
#' Object of class \code{LDO}.
"ldo"


#' 'LoBrA'  Data Object (LDO) for Example data set
#'
#' 'LoBrA'  example LDO created by the function 'createExampleData' and converted to an LDO by 'as.LOBdataset' function.
#' It consist of a single matrix for all experiments, time points, types (background, experiment), class and the intensity values of all components created.
#' The artificial data consist of 20 experiments and 100 components with 18 measurements (3 background, 15 sample). The 10 experiments are each associated to on of 2 groups (ONE and TWO). The components comprise 70 noise components and 30 components that randomly vary in their trajectories in one of three segments. Random noise is added to all intercepts, propagated and added to each time point for all samples and components separately.
#'
#' @name components
#' @docType data
#' @author Anne-Christin Hauschild \email{hauschild@@uni-marburg.de}
#' @keywords Breathomics Metabolomics Expression Longitudinal-Data Confounders
#' @format A vector of selected components from the longitudinal example data set.
"components"
