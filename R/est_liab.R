
# Here are functions that summarize genetic liabilities from pedigree records: 

#' PA-FGRS[mix], PA-FGRS[adt] and PA-FGRS[cont]  
#' 
#' Compute posterior mean and variance of liability of index person given disease status, thresholds and proportion of risk observed in relatives.
#' @param rel_status vector with status of relatives 1=affected, 0=not affected, NA= missing
#' @param thr threshold value used for all individuals (relatives and index person)
#' @param rel_thr optional vector of specific lower truncation for each relative
#' @param rel_w optional vector of proportion of risk experienced by each relative (used only in PA-FGRS[mix])
#' @param covmat covariance matrix of liabilities. First row and column correspoding to index person,subsequent rows and columns corresponding to relatives.
#' @param i_status optional status of index person, if to be considered 1=affected, 0=not affected, NA= missing 
#' @param i_thr optional specific threshold for index person 
#' @param i_w  optional proportion of risk experienced by index person (used only in PA-FGRS[mix])
#' @param conditional.mix logical indicator of wether the mixture of proportions should be conditional (default=TRUE).
#' @name pa_fgrs
#' @examples 
#' pa_fgrs(c(0,1),qnorm(.9),covmat = matrix(c(.5,.25,.25,.25,1,.25,.25,.25,1),3),rel_w=c(.5,1))
#' pa_fgrs_adt(c(0,1),qnorm(.9),covmat = matrix(c(.5,.25,.25,.25,1,.25,.25,.25,1),3))

#' @rdname pa_fgrs
#' @export
pa_fgrs_adt <- function(rel_status,thr=NA,rel_thr=rep(thr,length(rel_status)),
                        covmat,i_status=NA,i_thr=thr) {
  pa_fgrs2thr(rel_t1 = ifelse(rel_status,rel_thr,-Inf),rel_t2 = rel_thr,
              covmat = covmat,i_t1 = ifelse(i_status,i_thr,-Inf),i_t2 = i_thr)}


#' @rdname pa_fgrs
#' @export
pa_fgrs <- function(rel_status,thr=NA,rel_thr=rep(thr,length(rel_status)),rel_w=rep(1,length(rel_status)),
                        covmat,i_status=NA,i_thr=thr,i_w=1) {
  pa_fgrs2thr(rel_t1 = ifelse(rel_status,rel_thr,-Inf),rel_t2 = ifelse(rel_status,Inf,rel_thr),rel_w = rel_w,
              covmat = covmat,i_t1 = ifelse(i_status,i_thr,-Inf),i_t2 = ifelse(i_status,Inf,i_thr),i_w=i_w)}

#' @rdname pa_fgrs
#' @export
pa_fgrs_cont <- function(rel_value,covmat,i_value=NA) {
  pa_fgrs2thr(rel_t1 = rel_value,rel_t2 = rel_value,
              covmat = covmat,i_t1 = i_value,i_t2 = i_value)}

# this reduces to normal regression estimates:
#covmat = matrix(c(.7,.35,.35,.35,1,.25,.35,.35,1),3)
#pa_fgrs_cont(rel_value = c(1,3),covmat = covmat)
#covmat[-1,1]%*%solve(covmat[-1,-1])%*%c(1,3)
#covmat[1,1]-covmat[-1,1]%*%solve(covmat[-1,-1])%*%covmat[-1,1]

##### Helpers
trunc_norm_below <- function (mean, var, trunc) 
{
  mu = mean
  std.dev = sqrt(var)
  beta = (trunc - mu)/std.dev
  c(mean = mu - std.dev * dnorm(beta)/(pnorm(beta)), 
    var = std.dev^2 * (1 - beta * dnorm(beta)/pnorm(beta) - 
                               (dnorm(beta)/pnorm(beta))^2))
}
trunc_norm_above <- function (mean, var, trunc) 
{
  mu = mean
  std.dev = sqrt(var)
  alpha = (trunc - mu)/std.dev
  c(mean = mu + std.dev * dnorm(alpha)/(1 - pnorm(alpha)), 
    var = std.dev^2 * (1 + alpha * dnorm(alpha)/(1 - 
                                                         pnorm(alpha)) - (dnorm(alpha)/(1 - pnorm(alpha)))^2))
}

trunc_norm <- function(mean,var,trunc_l,trunc_u){
  if(trunc_l==trunc_u) out <- c(mean=ifelse(is.infinite(trunc_l),1e10,trunc_l),var=0) else
  if(is.infinite(trunc_l)) out <- trunc_norm_below(mean,var,trunc_u) else
  if(is.infinite(trunc_u)) out <- trunc_norm_above(mean,var,trunc_l) else{
    mu=mean
    std.dev = sqrt(var)
    alpha= (trunc_l-mu)/std.dev
    beta= (trunc_u-mu)/std.dev
  
    out <- c(
      mean=mu-std.dev*(dnorm(beta)-dnorm(alpha))/(pnorm(beta)-pnorm(alpha)),
      var=std.dev^2*(1-(beta*dnorm(beta)-alpha*dnorm(alpha))/
                     (pnorm(beta)-pnorm(alpha))-
                     ((dnorm(beta)-dnorm(alpha))/(pnorm(beta)-pnorm(alpha)))^2))
  }
  out
}
trunc_norm_mixture_conditional2  <- function(mean,var,trunc_l,trunc_u,Kp){
  mu=mean
  if(Kp==0) w_below= pnorm(trunc_u,mean = mean,sd = sqrt(var)) else
  if(trunc_u==Inf) w_below=1 else
  w_below= pnorm(trunc_u,mean = mean,sd = sqrt(var))/
    (1-pnorm(trunc_u,mean = mean,sd = sqrt(var),lower.tail = F)*
       Kp/pnorm(trunc_u,lower.tail = F))
  w_above=1-w_below
  m0=trunc_norm(mean=mu,var=var,trunc_l= trunc_l,trunc_u = trunc_u)
  m1=trunc_norm(mean=mu,var=var,trunc_l= trunc_u,trunc_u= Inf)
  new_mean=w_below*m0[1]+w_above*m1[1]
  c(mean=new_mean,
    var=w_below*(m0[1]^2+m0[2])+w_above*(m1[1]^2+m1[2])-new_mean^2)}

# This is an internal function that can both be used 
# for the standard PA-FGRS, for the PA-FGRS_adt model 
# and for a quantitative trait 
# specific functions for those three use cases are provided below
pa_fgrs2thr <- function(rel_t1=NA,rel_t2=NA,rel_w=rep(1,length(rel_t1)),covmat,i_t1=NA,i_t2=NA,i_w=1) {
  # remove NA's
  nonnas <- !is.na(rel_t1)&!is.na(rel_t2)&!is.na(rel_w)&!rel_w==0
  if(sum(nonnas)==0&is.na(i_t1+i_w))  c(postM=0,postVar=covmat[1,1]) else {
    if(sum(nonnas)==0) {
      rel_t1 <- i_t1
      rel_w <- i_w
      covmat <- matrix(c(rep(covmat[1,1],3),1),nrow=2)} else {
        rel_t1<-rel_t1[nonnas]
        rel_t2<-rel_t2[nonnas]
        rel_w <- rel_w[nonnas]
        covmat <- covmat[c(T,nonnas),c(T,nonnas),drop=F]
        
        # order 
        ord <- order(-rel_w,-covmat[-1,1],-rowSums(covmat[-1,,drop=F]))
        rel_t1 <- rel_t1[ord]
        rel_t2 <- rel_t2[ord]
        rel_w <- rel_w[ord]
        covmat <- covmat[c(1,ord+1),c(1,ord+1)]

        if(!is.na(i_t1+i_t2+i_w)){ # Add t1 and t2 of individual i if that's to be conditioned on
          rel_t1 <- c(i_t1,rel_t1)
          rel_t2 <- c(i_t2,rel_t2)
          rel_w <- c(i_w,rel_w)
          covmat <- rbind(c(covmat[1,1],covmat[1,1:ncol(covmat)]),
                          cbind(covmat[1:nrow(covmat),1],covmat))
          covmat[2,2] <- 1
        }}
    
    new_m <-t(t(rep(0,nrow(covmat))))   # E(l_j)
    new_cov <- covmat
    
    
    # Selection
    for(i in 1:length(rel_t1)){
      j <- nrow(new_cov)
      
      update= trunc_norm_mixture_conditional2(mean = new_m[j,],var = new_cov[j,j],
                                              trunc_l= rel_t1[j-1],
                                              trunc_u= rel_t2[j-1],
                                              Kp=rel_w[j-1]*(1-pnorm(rel_t2[j-1])))
      updated_m <-update[1]
      updated_var <-update[2]
      new_m  = new_m[-j,] + new_cov[-j,j]%*%solve(new_cov[j,j])%*%(updated_m-new_m[j,1])
      new_cov = new_cov[-j,-j,drop=F] -
        new_cov[-j,j]%*%(solve(new_cov[j,j])-solve(new_cov[j,j])%*%updated_var%*%solve(new_cov[j,j]))%*%new_cov[j,-j]
    }
    return(c(postM=new_m,postVar=new_cov))}
}

#' FGRS 
#' 
#' Our implementations of Kendler et al.

#' @param df data.frame with one row per individual with columns containing the phenotype information on the relatives
#' @param stat_cols vector of names of the columns containing the disease status of the relatives
#' @param age_cols vector of names of the columns containing the age of the relatives
#' @param rel_mat relatedness matrix with the first row and column correponding to the index person a subsequent rows to the relatives.
#' @param prev numeric indicating disease prevalence
#' @param env_cor optional numeric. Factor by which the first degree relatives are downweighted.
#' @param aoo a function that transforms the age values into a proportion of risk experienced by the individual
#' @export
FGRS <- function(df, stat_cols=which(grepl("stat",colnames(df))) ,age_cols=which(grepl("age",colnames(df))) , rel_mat, prev, env_cor=1, aoo){ 
  N=nrow(df)
  age= as.matrix(df[,age_cols])
  status= as.matrix(df[,stat_cols])
  
  w <- apply(age,2,aoo)
  w[!is.na(status)&status==1] <- 1 
  
  zhat_above = integrate(qnorm,1-prev,1)[[1]]/prev
  zhat_below = integrate(qnorm,0,1-prev)[[1]]/(1-prev)
  
  zhat  <- apply(status, 2, function(x) ifelse(x==1,zhat_above,zhat_below)) 
  
  df <- cbind(df, age= age, w= w , zhat=zhat)
  
  # number of relatives >0
  df$n_rel <- apply(df[,2:nrow(rel_mat)],1,function(x) sum(!is.na(x)))
  #df <- df[!df$n_rel==0,]
  
  stat_cols <- which(grepl("stat",colnames(df)))  
  zhat_cols <- which(grepl("zhat",colnames(df)))  
  w_cols <- which(grepl("w",colnames(df)))  
  pr_cols <- which(grepl("prop_risk",colnames(df)))  
  max_rel <- length(zhat_cols)
  
  # Summaries of family risk:
  df$sum_rzw <- rowSums(df[,zhat_cols]*df[,w_cols]*t(matrix(rel_mat[1,-1]*ifelse(rel_mat[1,-1]==0.5,env_cor,1),max_rel,N)),na.rm = T)
  
  df$sum_r <- rowSums(as.numeric(!is.na(df[,zhat_cols]))*t(matrix(rel_mat[1,-1],max_rel,N)),na.rm = T)
  
  df$mean_rzw <- df$sum_rzw/df$n_rel
  var_rel<-   var(unlist(df[,zhat_cols]),na.rm=T)  
  df$mean_rzw_shrunk <- df$mean_rzw*(var(df$mean_rzw,na.rm = T)/(var(df$mean_rzw,na.rm = T)+var_rel/df$sum_r))
  df$mean_rzw_shrunk  
}