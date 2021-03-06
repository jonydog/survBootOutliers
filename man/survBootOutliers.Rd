% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survBootOutliers.R
\name{survBootOutliers}
\alias{survBootOutliers}
\title{Extract the most outlying observations following a criteria based on the bootstrapped concordance with parallel processing}
\usage{
survBootOutliers(surv.object, covariate.data, sod.method, B, B.N = NULL,
  max.outliers, parallel.param = NULL)
}
\arguments{
\item{surv.object}{An obect of type survival::Surv containing lifetimes and right-censoring status}

\item{covariate.data}{A data frame containing the data with covariate values for each individual}

\item{sod.method}{One of c("osd","bht","dbht","ld","martingale","deviance")}

\item{B}{The number of bootstrap samples generated only applicable for "bht" and "dbht" methods. 
Typically at least 10x the size of the dataset, ideally should be increased until convergence.}

\item{B.N}{the number of observations in each bootstrap sample}

\item{max.outliers}{This parameter is only used for the "osd" method}

\item{parallel.param}{(Optional) A BiocParallel object, examples: SerialParam(), MulticoreParam()}
}
\value{
For all methods except for "bht" and "dbht" the value returned is a data.frame containing the most outlying observations sorted by outlying score.
        For the "bht" method the value returned is a list of two members: 
            "outlier_set": the most outlygin observations sorted by p-values; 
            "histograms": histogram of concordance variation for each observation.
        For the "dbht" method the value returned is a list of two members: 
            "outlier_set": the most outlygin observations sorted by p-values;
           "histograms": histogrms of concordance for each observations for the two types of bootstap: "poison" and "antidote".
}
\description{
Extract the most outlying observations following a criteria based on the bootstrapped concordance with parallel processing
}
\examples{
## One Step Deletion "osd" method
\dontrun{
whas <- get.whas100.dataset()
print( getwd() )
outliers_osd <- survBootOutliers( 
   surv.object=Surv(time = whas$times,event = whas$status ), 
   covariate.data = whas[,2:5], 
   sod.method = "osd", 
   max.outliers = 5
 )
}

## Bootstrap Hypothesis Test "bht" with 1000 bootstrap samples, 
## each with 100 individuals, running without parallelism.
\dontrun{  
whas <- get.whas100.dataset()
outliers_bht <- survBootOutliers( 
     surv.object=Surv(time = whas$times,event = whas$status ), 
     covariate.data = whas[,2:5], 
     sod.method = "bht", 
     B = 1000, 
     B.N = 100,
     parallel.param = BiocParallel::MulticoreParam() 
)
}

## Dual Bootstrap Hypothesis Test "dbht" with 1000 bootstrap samples,
## each with 50 individuals and running on all available cores.
\dontrun{ whas <- get.whas100.dataset()
outliers_dbht <- survBootOutliers( 
   surv.object=Surv(time = whas$times,event = whas$status ), 
   covariate.data = whas[,2:5], 
   sod.method = "dbht",
   B = 1000, 
   B.N = 50,
   parallel.param = BiocParallel::MulticoreParam() 
)
}
## One Step Deletion "osd" with an amount of 10 for maximum outlier count
whas <- get.whas100.dataset()
outliers_osd <- survBootOutliers( 
   surv.object=Surv(time = whas$times,event = whas$status ), 
   covariate.data = whas[,2:5], 
   sod.method = "osd", 
   max.outliers = 10
)

## Likelihood displacement criterion for outlier ranking
whas <- get.whas100.dataset()
outliers_ld <- survBootOutliers( 
   surv.object=Surv(time = whas$times,event = whas$status ), 
   covariate.data = whas[,2:5], 
   sod.method = "ld"
)

## Cox regression deviance residuals criterion for outlier ranking
whas <- get.whas100.dataset()
outliers_deviance <- survBootOutliers( 
   surv.object=Surv(time = whas$times,event = whas$status ), 
   covariate.data = whas[,2:5], 
   sod.method = "deviance"
)

## Cox regression Martingale residuals criterion for outlier ranking
whas <- get.whas100.dataset()
outliers_martingale <- survBootOutliers( 
   surv.object=Surv(time = whas$times,event = whas$status ), 
   covariate.data = whas[,2:5], 
   sod.method = "martingale"
)

}
