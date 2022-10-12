## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- PAFGRS------------------------------------------------------------------
library("PAFGRS")

## ---- data--------------------------------------------------------------------
library("FamAgg")
data("minnbreast")

head(minnbreast[minnbreast$famid==4,1:7])

## ---- kinship-----------------------------------------------------------------
ped_data <- minnbreast[minnbreast$famid==4,c(1:7,14)]
# generate pedigree
ped <- pedigree(id = ped_data$id,dadid = ped_data$fatherid,momid = ped_data$motherid,
                sex = ped_data$sex,affected = ped_data$cancer&ped_data$sex=="F", status=)
# compute kinsip matrix: 
k <- kinship2::kinship(ped)

plot(ped)

## ---- covmat------------------------------------------------------------------
r <- k*2
h2= 0.5 
covmat <- h2*r 
diag(covmat) <- 1

## ---- fgrs-full-obs-----------------------------------------------------------
covmat <- covmat[c(3,(1:nrow(covmat))[-3]),c(3,(1:nrow(covmat))[-3])]
covmat[1,1] <- h2 
ped_data[ped_data$sex=="M",]$cancer <- 0 
status <- ped_data[-3,]$cancer
prev <- .05
thr <- qnorm(1-prev)

pa_fgrs(rel_status = status, covmat = covmat ,thr = thr)

## ---- fgrs-sex----------------------------------------------------------------
ped_data$prev <- .1
ped_data$prev[ped_data$sex=="M"] <- 0
ped_data$thr <- qnorm(1-ped_data$prev)

pa_fgrs(rel_status = status, covmat = covmat,rel_thr = ped_data$thr[-3])

## ---- fgrs-age----------------------------------------------------------------
plot(sapply(1:100,function(x) {
  w=(x -20 )/80
  ifelse(w<0,0,w)*.1}),type = "l",ylab = "Cumulative Incidence in Females",xlab = "Age")
ped_data$w <-(ped_data$endage -20 )/80
ped_data$w[ped_data$w<0] <- 0
ped_data$w[ped_data$cancer==1] <- 1

pa_fgrs(rel_status = status, covmat = covmat,rel_thr = ped_data$thr[-3],rel_w = ped_data$w[-3])

