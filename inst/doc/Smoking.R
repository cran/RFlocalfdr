## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
knitr::opts_chunk$set(fig.width=6, fig.height=6, dpi=300,echo = FALSE)
knitr::opts_chunk$set(fig.pos = "H", out.extra = "")


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  library(RFlocalfdr.data)
#  data(smoking)
#  ?smoking
#  y<-smoking$y
#  smoking_data<-smoking$rma
#  y.numeric <-ifelse((y=="never-smoked"),0,1)
#  

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  install.packages("geneExpressionFromGEO")
#  library(geneExpressionFromGEO)
#  DF1 <- getGeneExpressionFromGEO("GSE994", FALSE, FALSE)  #retrieveGeneSymbols, verbose = FALSE)
#  
#  aa<-match(gsub("([A-Z0-9]*).*","\\1", rownames(smoking_data)), colnames(DF1))
#  all.equal(aa,1:57)
#  $[1] TRUE
#  
#  temp <- t(DF1[1:57])
#  #but class temp is a data.frame not a ‘AffyBatch’ object. So how do we normalize it?
#  
#  ## library(affy)
#  ## rma.data <- rma(data)
#  ## or
#  ## rma.data <- expresso(data,
#  ## bgcorrect.method = "rma",
#  ## normalize.method = "quantiles",
#  ## pmcorrect.method = "pmonly",
#  ## summary.method = "medianpolish")
#  

## ----  echo=TRUE,  eval=FALSE-------------------------------------------------
#  library(ranger)
#  rf1 <-ranger(y=y.numeric ,x=smoking_data,importance="impurity",seed=123, num.trees = 10000,
#               classification=TRUE)
#  t2 <-count_variables(rf1)
#  imp<-log(rf1$variable.importance)
#  #png("./supp_figures/smoking_log_importances.png")
#  plot(density(imp),xlab="log importances",main="")
#  #dev.off()
#  

## ----simulation2, echo=FALSE, fig.cap="A small simulated data set. Each band contains blocks of size {1, 2, 4, 8, 16, 32, 64}, and each block consists of correlated (identical variables).", fig.align="center", out.width = '50%'----
knitr::include_graphics("./supp_figures/smoking_log_importances.png")


## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  cutoffs <- c(2,3,4,5)
#  #png("./supp_figures/smoking_data_determine_cutoff.png")
#  res.con<- determine_cutoff(imp,t2,cutoff=cutoffs,plot=c(2,3,4,5))
#  #dev.off()
#  
#  #png("./supp_figures/smoking_data_determine_cutoffs_2.png")
#  plot(cutoffs,res.con[,3],pch=15,col="red",cex=1.5,ylab="max(abs(y - t1))")
#  #dev.off()
#  cutoffs[which.min(res.con[,3])]
#  

## ---- echo=TRUE,eval=FALSE----------------------------------------------------
#  temp<-imp[t2 > 3]
#  temp <- temp - min(temp) + .Machine$double.eps
#  qq <- plotQ(temp,debug.flag = 1)
#  ppp<-run.it.importances(qq,temp,debug.flag = 0)
#  
#  
#  png("./supp_figures/smoking_significant_genes.png")
#  #aa<-significant.genes(ppp,temp,cutoff=0.05,do.plot=1)
#  aa<-significant.genes(ppp,temp,cutoff=0.05,debug.flag=0,do.plot=TRUE,use_95_q=TRUE)
#  dev.off()
#  length(aa$probabilities) # 17
#  
#  aa<-significant.genes(ppp,temp,cutoff=0.05,debug.flag=0,do.plot=TRUE,use_95_q=FALSE)
#  length(aa$probabilities) # 19
#  
#  
#  

## ---- echo=TRUE,eval=FALSE----------------------------------------------------
#  sessionInfo()
#  

## ---- echo=TRUE,eval=FALSE----------------------------------------------------
#  devtools::install_github("parsifal9/RFlocalfdr", build_vignettes = TRUE, force = TRUE)
#  
#  library(RFlocalfdr)
#  data(smoking)
#  ?smoking
#  y<-smoking$y
#  smoking_data<-smoking$rma
#  y.numeric <-ifelse((y=="never-smoked"),0,1)
#  
#  library(ranger)
#  rf1 <-ranger(y=y.numeric ,x=smoking_data,importance="impurity",seed=123, num.trees = 10000,
#               classification=TRUE)
#  t2 <-count_variables(rf1)
#  imp<-log(rf1$variable.importance)
#  #png("./supp_figures/smoking_log_importances.png")
#  plot(density(imp),xlab="log importances",main="")
#  #dev.off()
#  
#  
#  
#  cutoffs <- c(2,3,4,5)
#  #png("./supp_figures/smoking_data_determine_cutoff.png")
#  res.con<- determine_cutoff(imp,t2,cutoff=cutoffs,plot=c(2,3,4,5))
#  #dev.off()
#  
#  #png("./supp_figures/smoking_data_determine_cutoffs_2.png")
#  plot(cutoffs,res.con[,3],pch=15,col="red",cex=1.5,ylab="max(abs(y - t1))")
#  #dev.off()
#  cutoffs[which.min(res.con[,3])]
#  
#  
#  temp<-imp[t2 > 3]
#  temp <- temp - min(temp) + .Machine$double.eps
#  qq <- plotQ(temp,debug.flag = 1)
#  ppp<-run.it.importances(qq,temp,debug.flag = 0)
#  
#  
#  #png("./supp_figures/smoking_significant_genes.png")
#  #aa<-significant.genes(ppp,temp,cutoff=0.05,do.plot=1)
#  aa<-significant.genes(ppp,temp,cutoff=0.05,debug.flag=0,do.plot=TRUE,use_95_q=TRUE)
#  #dev.off()
#  length(aa$probabilities) # 17  -- Roc gets 30
#  
#  aa<-significant.genes(ppp,temp,cutoff=0.05,debug.flag=0,do.plot=TRUE,use_95_q=FALSE)
#  length(aa$probabilities) # 19-- Roc gets 30
#  

