\name{display.obs.histogram}
\alias{display.obs.histogram}
\title{Auxiliar function that displays the concordance histogram associated with the observation.}
\usage{
display.obs.histogram(histograms, type, obs.index)
}
\description{
  This function plots the histogram of the histogram or histograms that result in the p-values associated with each observation. For type "bht"
  the histogram contains the counts of the concordance displacement when fiting the model to each bootstrap resample without the target observation. For type "dbht" there are   two histograms present: poison and antidore histograms. The poison histograms is calculated with each boostrap resample containing the target observation, and the antidote the same histogram as in "bht".
}
\arguments{
  \item{histograms}{
    The histograms object returned by the survBootOutliers function when the method selected is "bht" or "dbht".
  }
  \item{type}{
    The type of histogram that is given as input, possible choices are again "bht" or "dbht".
  }
  \item{obs.index}{
    The original index of the observation of the concordance histograms to be displayed
  }
}
\value{
  There is no return value
}
\examples{
whas <- get.whas100.dataset()
## using "bht"" as outlier detection method
outliers_bht <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "bht" , B = 10 , B.N = 100  )
## now display the histogram of observation 67 (for example) for the 
display.obs.histogram(outliers_bht$histograms, "bht", 67)
## using "dbht"" as outlier detection method
whas <- get.whas100.dataset()
outliers_dbht <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "dbht" , B = 10 , B.N = 100  )
## now display the dual histograms again for observation 67
display.obs.histogram(outliers_dbht$histograms, "dbht", 67)
}
