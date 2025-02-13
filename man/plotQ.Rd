% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotQ.R
\name{plotQ}
\alias{plotQ}
\title{plotQ}
\usage{
plotQ(imp, debug.flag = 0, temp.dir = NULL, try.counter = 3, start_at = 37)
}
\arguments{
\item{imp}{"reduction in impurity" importances from a random forest model}

\item{debug.flag}{either 0 (no debugging information), 1 or 2}

\item{temp.dir}{if debug flag is >0 then information is written to temp.dir}

\item{try.counter}{where to explain this?}

\item{start_at}{--  x <- f_fit$midpoints  is of length 119 (quite arbitrary). We use the first start_at
values of x to fit the skew-normal distribution.}
}
\value{
\itemize{
\item df, contains x and y, midpoints and counts from a histogram of imp
\item final.estimates_C_0.95, the output from the fitting routine nlsLM in minpack.lm, where df has been truncated
at the value C_0.95 (the 0.95 quantile of the skew-Normal distribution fitted to the imp histogram
\item final.estimates_cc, as for final.estimates_C_0.95 but with cc determined by the procedure of Gauran et.al 2018
\item temp.dir,  the directory where debugging information may be written
\item C_0.95, the 0.95 quantile of the skew-Normal distribution fitted to the imp histogram
\item cc,  determined by the procedure of Gauran et.al 2018
\item fileConn, a file connectin for writing debugging information
\item f_fit, a spline fit to the histogram
\item ww the minimum value of the local fdr
}
}
\description{
produces a plot showing the q values
\itemize{
\item q_95, the 95th quantile of the data
\item q using the penalized selection method of Gauran et.al 2018
}
}
\details{
We estiamte a value "q" such that:
to the left of "q", the density is composed solely of NULL importance values
to the right of "q" we have  a density that is a mixture of null and non-null  importance values.
The method of Gauran et.al 2018 may not work in cases where the data distribution is not well modelled by a skew-normal.
The q_95 value can be uses as a workaround in these case.
In many cases they will be very similar
}
\examples{
data(imp20000)
imp <- log(imp20000$importances)
t2 <- imp20000$counts
plot(density((imp)))
hist(imp,col=6,lwd=2,breaks=100,main="histogram of importances")
res.temp <- determine_cutoff(imp, t2 ,cutoff=c(0,1,2,3),plot=c(0,1,2,3),Q=0.75,try.counter=1)
plot(c(0,1,2,3),res.temp[,3])
imp<-imp[t2 > 1]
qq <- plotQ(imp,debug.flag = 0)                                                          
ppp<-run.it.importances(qq,imp,debug=0)                                                       
aa<-significant.genes(ppp,imp,cutoff=0.2,debug.flag=0,do.plot=2, use_95_q=TRUE)                           
length(aa$probabilities) #11#                                                          
names(aa$probabilities)

\donttest{
library(RFlocalfdr.data)
data(ch22)                                                                                 
?ch22                                                                                     
#document how the data set is created                                                      
plot(density(log(ch22$imp)))                                                               
t2 <-ch22$C                                                                                
imp<-log(ch22$imp)                                                                         
#Detemine a cutoff to get a unimodal density.
# This was calculated previously. See determine_cutoff
imp<-imp[t2 > 30]
qq <- plotQ(imp,debug.flag = 0)

data(smoking)
?smoking 
y<-smoking$y
smoking_data<-smoking$rma
y.numeric <-ifelse((y=="never-smoked"),0,1)

library(ranger)
rf1 <-ranger::ranger(y=y.numeric ,x=smoking_data,importance="impurity",seed=123, num.trees = 10000,
             classification=TRUE)
t2 <-count_variables(rf1)
imp<-log(rf1$variable.importance)
plot(density(imp),xlab="log importances",main="")
cutoffs <- c(2,3,4,5)
res.con<- determine_cutoff(imp,t2,cutoff=cutoffs,plot=c(2,3,4,5))

plot(cutoffs,res.con[,3],pch=15,col="red",cex=1.5,ylab="max(abs(y - t1))")
cutoffs[which.min(res.con[,3])]

temp<-imp[t2 > 3]
temp <- temp - min(temp) + .Machine$double.eps
qq <- plotQ(temp)
ppp<-run.it.importances(qq,temp,debug.flag = 0)
aa<-significant.genes(ppp,temp,cutoff=0.05,debug.flag=0,do.plot=TRUE,use_95_q=TRUE)
length(aa$probabilities) # 17
}
}
\keyword{importance}
\keyword{variable}
