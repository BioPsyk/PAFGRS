# Here are functions that summarize genetic liabilities from pedigree records: 

#' PA-FGRS 
#' 
#' Compute posterior mean and variance of liability of index person given disease status, thresholds and proportion of risk observed in relatives.
#' @param rel_status vector with status of relatives 1=affected, 0=not affected, NA= missing
#' @param thr threshold value used for all individuals (relatives and index person)
#' @param rel_thr optional vector of specific threshold for each relative
#' @param rel_w optional vector of proportion of risk experienced by each relative
#' @param covmat covariance matrix of liabilities. First row and column correspoding to index person,subsequent rows and columns corresponding to relatives.
#' @param i_status optional status of index person, if to be considered 1=affected, 0=not affected, NA= missing 
#' @param i_thr optional specific threshold for index person 
#' @param i_w  optional proportion of risk experienced by index person
#' @examples 
#' pa_fgrs(c(0,1),qnorm(.9),covmat = matrix(c(.5,.25,.25,.25,1,.25,.25,.25,1),3))
#' 
#' @importFrom stats dnorm 
#' @importFrom stats integrate
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats var
#' @export
pa_fgrs <- function(rel_status,thr=NA,rel_thr=rep(thr,length(rel_status)),rel_w=rep(1,length(rel_status)),covmat,i_status=NA,i_thr=thr,i_w=1,conditional.mix=T) {
  # remove NA's
  nonnas <- !is.na(rel_status)&!is.infinite(rel_thr)&!is.na(rel_w)&!rel_w==0
  if(sum(nonnas)==0&is.na(i_status+i_thr+i_w))  c(postM=0,postVar=covmat[1,1]) else {
    if(sum(nonnas)==0) {
      rel_status <- i_status
      rel_thr <- i_thr
      rel_w <- i_w
      covmat <- matrix(c(rep(covmat[1,1],3),1),nrow=2)} else {
        rel_status<-rel_status[nonnas]
        rel_thr<-rel_thr[nonnas]
        rel_w <- rel_w[nonnas]
        covmat <- covmat[c(T,nonnas),c(T,nonnas),drop=F]
        # order 
        ord <- order(-rel_w,-covmat[-1,1],-rowSums(covmat[-1,,drop=F]))
        rel_status <- rel_status[ord]
        rel_thr <- rel_thr[ord]
        rel_w <- rel_w[ord]
        covmat <- covmat[c(1,ord+1),c(1,ord+1)]
        
        if(!is.na(i_status+i_thr+i_w)){ # Add status of individual i if that's to be conditioned on
          rel_status <- c(i_status,rel_status)
          rel_thr <- c(i_thr,rel_thr)
          rel_w <- c(i_w,rel_w)
          covmat <- rbind(c(covmat[1,1],covmat[1,1:ncol(covmat)]),
                          cbind(covmat[1:nrow(covmat),1],covmat))
          covmat[2,2] <- 1
        }}
    rel_w[rel_status==1] <- 1

    new_m <-t(t(rep(0,nrow(covmat))))   # E(l_j)
    new_cov <- covmat
    
if(conditional.mix)    trunc_norm_mixture_fun <-trunc_norm_mixture_conditional else
  trunc_norm_mixture_fun <-trunc_norm_mixture
    
    # Selection
    for(i in 1:length(rel_status)){
      j <- nrow(new_cov)
            updated_m <-ifelse(rel_status[j-1]==1,
                               trunc_norm_above(mean = new_m[j,],var = new_cov[j,j],
                                                trunc = rel_thr[j-1] )[1]
                         ,trunc_norm_mixture_fun(mean = new_m[j,],var = new_cov[j,j],
                                             trunc = rel_thr[j-1],Kp =rel_w[j-1]*
                                               (1-pnorm(rel_thr[j-1]))  )[1])
      updated_var <-ifelse(rel_status[j-1]==1,
                           trunc_norm_above(mean = new_m[j,],var = new_cov[j,j],
                                            trunc = rel_thr[j-1] )[2]
                           ,trunc_norm_mixture_fun(mean = new_m[j,],var = new_cov[j,j],
                                               trunc = rel_thr[j-1],Kp =rel_w[j-1]*
                                                 (1-pnorm(rel_thr[j-1]))  )[2])
      new_m  = new_m[-j,] + new_cov[-j,j]%*%solve(new_cov[j,j])%*%(updated_m-new_m[j,1])
      new_cov = new_cov[-j,-j,drop=F] -
        new_cov[-j,j]%*%(solve(new_cov[j,j])-solve(new_cov[j,j])%*%updated_var%*%solve(new_cov[j,j]))%*%new_cov[j,-j]
    }
    return(c(postM=new_m,postVar=new_cov))}
}

##### Helpers

# function that estimates the expected mean and variance of a truncated normal distribution with left truncation
trunc_norm_above <- function(mean,var,trunc){
  mu=mean
  std.dev = sqrt(var)
  alpha= (trunc-mu)/std.dev
  c(
  mean_above=mu+std.dev*dnorm(alpha)/(1-pnorm(alpha)),
  var_above=std.dev^2*(1+alpha*dnorm(alpha)/(1-pnorm(alpha))-(dnorm(alpha)/(1-pnorm(alpha)))^2)
  )
}

# function that estimates the expected mean and variance of a truncated normal distribution with right truncation
trunc_norm_below <- function(mean,var,trunc){
  mu=mean
  std.dev = sqrt(var)
  beta= (trunc-mu)/std.dev
  c(
    mean_below=mu-std.dev*dnorm(beta)/(pnorm(beta)),
    var_below=std.dev^2*(1-beta*dnorm(beta)/pnorm(beta)-(dnorm(beta)/pnorm(beta))^2)
  )
}


#https://en.wikipedia.org/wiki/Mixture_distribution#Moments 

trunc_norm_mixture <- function(mean,var,trunc,Kp){
  mu=mean
  w_below= (pnorm(trunc))/(1-Kp)
  w_above=1-w_below
  
    new_mean=w_below*trunc_norm_below(mean=mu,var=var,trunc = trunc)[1]+
      w_above*trunc_norm_above(mean=mu,var=var,trunc = trunc)[1]
    c(mean=new_mean,
      var=w_below*(trunc_norm_below(mean=mu,var=var,trunc = trunc)[1]^2+
      trunc_norm_below(mean=mu,var=var,trunc = trunc)[2])+
      w_above*(trunc_norm_above(mean=mu,var=var,trunc = trunc)[1]^2+
      trunc_norm_above(mean=mu,var=var,trunc = trunc)[2])-new_mean^2
  )
}

trunc_norm_mixture_conditional  <- function(mean,var,trunc,Kp){
  mu=mean
  w_below= pnorm(trunc,mean = mean,sd = sqrt(var))/
    (1-pnorm(trunc,mean = mean,sd = sqrt(var),lower.tail = F)*
       Kp/pnorm(trunc,lower.tail = F))
  w_above=1-w_below
  new_mean=w_below*trunc_norm_below(mean=mu,var=var,trunc = trunc)[1]+
      w_above*trunc_norm_above(mean=mu,var=var,trunc = trunc)[1]
      c(mean=new_mean,
        var=w_below*(trunc_norm_below(mean=mu,var=var,trunc = trunc)[1]^2+
                       trunc_norm_below(mean=mu,var=var,trunc = trunc)[2])+
          w_above*(trunc_norm_above(mean=mu,var=var,trunc = trunc)[1]^2+
                     trunc_norm_above(mean=mu,var=var,trunc = trunc)[2])-new_mean^2
  )
}


#' FGRS 
#' 
#' Our implementations of Kendler et al.

#' @param df data.frame TODO
#' @param stat_cols TODO
#' @param age_cold TODO
#' @param rel_mat TODO
#' @param prev TODO
#' @param env_cor TODO
#' @param aoo TODO
#' @param ped TODO
#' @export
FGRS <- function(df, stat_cols=which(grepl("stat",colnames(df))) ,age_cols=which(grepl("age",colnames(df))) , rel_mat, prev, env_cor, aoo,ped){ 
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




pa_thr_hier_old <- function(rel_status,rel_thr,covmat,i_status=NA,i_thr=NA) {
  # remove NA's
  nonnas <- !is.na(rel_status)&!is.infinite(rel_thr)
  if(sum(nonnas)==0)  c(postM=0,postVar=covmat[1,1]) else{
    rel_status<-rel_status[nonnas]
    rel_thr<-rel_thr[nonnas]
    covmat <- covmat[c(T,nonnas),c(T,nonnas),drop=F]
    
    # order 
    ord <- order(rel_thr,-covmat[-1,1],-rowSums(covmat[-1,,drop=F]))
    rel_status <- rel_status[ord]
    rel_thr <- rel_thr[ord]
    covmat <- covmat[c(1,ord+1),c(1,ord+1)]
    #  
    updated_m <- t(t(c(0,sapply(1:length(rel_status), function(x){ # E(l_j|Y_j,T_j)
      ifelse(rel_status[x]==1, 
             lambda_fun_above(rel_thr[x]),
             lambda_fun_below(rel_thr[x])
      )}) )))
    new_m <-t(t(rep(0,nrow(covmat))))   # E(l_j)
    orig_cov <- new_cov <- covmat
    
    delta = c(0,sapply(1:length(rel_status), function(x) # 
      ifelse(rel_status[x]==1, 
             delta_fun( lambda_fun_above(rel_thr[x]),rel_thr[x]),
             delta_fun( lambda_fun_below(rel_thr[x]),rel_thr[x]))) )
    
    # delta = c(0,sapply(1:length(rel_status), function(x) # 
    #   ifelse(rel_status[x]==1, 
    #          delta_fun( lambda_fun_above(rel_thr[x]),rel_thr[x]),
    #          ifelse(is.infinite(rel_thr[x]),0,delta_fun( lambda_fun_below(rel_thr[x]),rel_thr[x]))) ))
    # 
    
    updated_var <- diag(covmat)-delta
    
    # Selection
    for(i in 1:length(rel_status)){
      j <- nrow(new_cov)
      updated_cov = new_cov[-j,-j,drop=F] -
        new_cov[-j,j]%*%(solve(new_cov[j,j])-solve(new_cov[j,j])%*%updated_var[j]%*%solve(new_cov[j,j]))%*%new_cov[j,-j]
      new_m  = new_m[-j,] + new_cov[-j,j]%*%solve(new_cov[j,j])%*%(updated_m[j,1]-new_m[j,1])
      new_cov=updated_cov
      if(j>2){
        alpha = rel_thr[j-1] - new_m[j-1] / sqrt(covmat[j-1,j-1])
        lambda_alpha = ifelse(rel_status[j-2] == 1, lambda_fun_above(alpha),lambda_fun_below(alpha))
        #delta_alpha =  ifelse(is.infinite(alpha),0,delta_fun(lambda_alpha,alpha))
        delta_alpha =  delta_fun(lambda_alpha,alpha)
        updated_m[j-1,] = new_m[j-1,1] + sqrt(updated_cov[j-1,j-1])*lambda_alpha
        updated_var[j-1] = updated_cov[j-1,j-1]*(1-delta_alpha)}
      # update mean
      #   new_mean = mean_matrix[2,1]+ cov_matrix[2,1]%*%solve(cov_matrix[1,1])%*%(mean_after_selec-mean_matrix[1,1])
      # update variance
      
    }
    
    
    c(postM=new_m,postVar=new_cov)}
}



pa_thr_hier <- function(rel_status,thr,rel_thr=rep(thr,length(rel_status)),rel_w=rep(1,length(rel_status)),covmat,i_status=NA,i_thr=NA,i_w=NA) {
  # remove NA's
  nonnas <- !is.na(rel_status)&!is.infinite(rel_thr)&!is.na(rel_w)&!rel_w==0
  if(sum(nonnas)==0&is.na(i_status+i_thr+i_w))  c(postM=0,postVar=covmat[1,1]) else {
    if(sum(nonnas)==0) {
      rel_status <- i_status
      rel_thr <- i_thr
      rel_w <- i_w
      covmat <- matrix(c(rep(covmat[1,1],3),1),nrow=2)} else {
        rel_status<-rel_status[nonnas]
        rel_thr<-rel_thr[nonnas]
        rel_w <- rel_w[nonnas]
        covmat <- covmat[c(T,nonnas),c(T,nonnas),drop=F]
        # order 
        ord <- order(-rel_w,-covmat[-1,1],-rowSums(covmat[-1,,drop=F]))
        rel_status <- rel_status[ord]
        rel_thr <- rel_thr[ord]
        rel_w <- rel_w[ord]
        covmat <- covmat[c(1,ord+1),c(1,ord+1)]
        
        if(!is.na(i_status+i_thr+i_w)){ # Add status of individual i if that's to be conditioned on
          rel_status <- c(i_status,rel_status)
          rel_thr <- c(i_thr,rel_thr)
          rel_w <- c(i_w,rel_w)
          covmat <- rbind(c(covmat[1,1],covmat[1,1:ncol(covmat)]),
                          cbind(covmat[1:nrow(covmat),1],covmat))
          covmat[2,2] <- 1
        }}
    rel_w[rel_status==1] <- 1
    #  
    updated_m <- t(t(c(0,sapply(1:length(rel_status), function(x){ # E(l_j|Y_j,T_j)
      ifelse(rel_status[x]==1, 
             lambda_fun_above(rel_thr[x]),
             -lambda_fun_above(rel_thr[x])/(1/(pnorm(-rel_thr[x])*rel_w[x])-1)
      )}) )))
    new_m <-t(t(rep(0,nrow(covmat))))   # E(l_j)
    new_cov <- covmat
    
    updated_var = c(0,sapply(1:length(rel_status), function(x) # 
      1+rel_thr[x]*updated_m[x+1]-updated_m[x+1]^2
    )) #https://en.wikipedia.org/wiki/Mixture_distribution#Moments 
    
    
    # Selection
    for(i in 1:length(rel_status)){
      j <- nrow(new_cov)
      new_m  = new_m[-j,] + new_cov[-j,j]%*%solve(new_cov[j,j])%*%(updated_m[j,1]-new_m[j,1])
      new_cov = new_cov[-j,-j,drop=F] -
        new_cov[-j,j]%*%(solve(new_cov[j,j])-solve(new_cov[j,j])%*%updated_var[j]%*%solve(new_cov[j,j]))%*%new_cov[j,-j]
    }
    c(postM=new_m,postVar=new_cov)}
}

pa_thr_hier_2022 <- function(rel_status,thr,rel_thr=rep(thr,length(rel_status)),rel_w=rep(1,length(rel_status)),covmat,i_status=NA,i_thr=NA,i_w=NA) {
  # remove NA's
  nonnas <- !is.na(rel_status)&!is.infinite(rel_thr)&!is.na(rel_w)&!rel_w==0
  if(sum(nonnas)==0&is.na(i_status+i_thr+i_w))  c(postM=0,postVar=covmat[1,1]) else {
    if(sum(nonnas)==0) {
      rel_status <- i_status
      rel_thr <- i_thr
      rel_w <- i_w
      covmat <- matrix(c(rep(covmat[1,1],3),1),nrow=2)} else {
        rel_status<-rel_status[nonnas]
        rel_thr<-rel_thr[nonnas]
        rel_w <- rel_w[nonnas]
        covmat <- covmat[c(T,nonnas),c(T,nonnas),drop=F]
        # order 
        ord <- order(-rel_w,-covmat[-1,1],-rowSums(covmat[-1,,drop=F]))
        rel_status <- rel_status[ord]
        rel_thr <- rel_thr[ord]
        rel_w <- rel_w[ord]
        covmat <- covmat[c(1,ord+1),c(1,ord+1)]
        
        if(!is.na(i_status+i_thr+i_w)){ # Add status of individual i if that's to be conditioned on
          rel_status <- c(i_status,rel_status)
          rel_thr <- c(i_thr,rel_thr)
          rel_w <- c(i_w,rel_w)
          covmat <- rbind(c(covmat[1,1],covmat[1,1:ncol(covmat)]),
                          cbind(covmat[1:nrow(covmat),1],covmat))
          covmat[2,2] <- 1
        }}
    rel_w[rel_status==1] <- 1
    #  
    updated_m <- t(t(c(0,sapply(1:length(rel_status), function(x){ # E(l_j|Y_j,T_j)
      ifelse(rel_status[x]==1, 
             lambda_fun_above(rel_thr[x]),
             -lambda_fun_above(rel_thr[x])/(1/(pnorm(-rel_thr[x])*rel_w[x])-1)
      )}) )))
    new_m <-t(t(rep(0,nrow(covmat))))   # E(l_j)
    new_cov <- covmat
    
    updated_var = c(0,sapply(1:length(rel_status), function(x) # 
      1+rel_thr[x]*updated_m[x+1]-updated_m[x+1]^2
    )) #https://en.wikipedia.org/wiki/Mixture_distribution#Moments 
    
    
    # Selection
    for(i in 1:length(rel_status)){
      j <- nrow(new_cov)
      new_m  = new_m[-j,] + new_cov[-j,j]%*%solve(new_cov[j,j])%*%(updated_m[j,1]*sqrt(new_cov[j,j]))
      new_cov = new_cov[-j,-j,drop=F] -
        new_cov[-j,j]%*%(solve(new_cov[j,j])-updated_var[j]%*%solve(new_cov[j,j]))%*%new_cov[j,-j]
    }
    c(postM=new_m,postVar=new_cov)}
}


pa_thr_hier_2022_2 <- function(rel_status,thr,rel_thr=rep(thr,length(rel_status)),rel_w=rep(1,length(rel_status)),covmat,i_status=NA,i_thr=NA,i_w=NA) {
  # remove NA's
  nonnas <- !is.na(rel_status)&!is.infinite(rel_thr)&!is.na(rel_w)&!rel_w==0
  if(sum(nonnas)==0&is.na(i_status+i_thr+i_w))  c(postM=0,postVar=covmat[1,1]) else {
    if(sum(nonnas)==0) {
      rel_status <- i_status
      rel_thr <- i_thr
      rel_w <- i_w
      covmat <- matrix(c(rep(covmat[1,1],3),1),nrow=2)} else {
        rel_status<-rel_status[nonnas]
        rel_thr<-rel_thr[nonnas]
        rel_w <- rel_w[nonnas]
        covmat <- covmat[c(T,nonnas),c(T,nonnas),drop=F]
        # order 
        ord <- order(-rel_w,-covmat[-1,1],-rowSums(covmat[-1,,drop=F]))
        rel_status <- rel_status[ord]
        rel_thr <- rel_thr[ord]
        rel_w <- rel_w[ord]
        covmat <- covmat[c(1,ord+1),c(1,ord+1)]
        
        if(!is.na(i_status+i_thr+i_w)){ # Add status of individual i if that's to be conditioned on
          rel_status <- c(i_status,rel_status)
          rel_thr <- c(i_thr,rel_thr)
          rel_w <- c(i_w,rel_w)
          covmat <- rbind(c(covmat[1,1],covmat[1,1:ncol(covmat)]),
                          cbind(covmat[1:nrow(covmat),1],covmat))
          covmat[2,2] <- 1
        }}
    rel_w[rel_status==1] <- 1
    #  
    updated_m <- t(t(c(0,sapply(1:length(rel_status), function(x){ # E(l_j|Y_j,T_j)
      ifelse(rel_status[x]==1, 
             lambda_fun_above(rel_thr[x]),
             -lambda_fun_above(rel_thr[x])/(1/(pnorm(-rel_thr[x])*rel_w[x])-1)
      )}) )))
    new_m <-t(t(rep(0,nrow(covmat))))   # E(l_j)
    new_cov <- covmat
    
    updated_var = c(0,sapply(1:length(rel_status), function(x) # 
      1+rel_thr[x]*updated_m[x+1]-updated_m[x+1]^2
    )) #https://en.wikipedia.org/wiki/Mixture_distribution#Moments 
    
    
    # Selection
    for(i in 1:length(rel_status)){
      j <- nrow(new_cov)
      lambda_alpha <- ifelse(rel_status[j-1]==1,
                             lambda_fun_above((rel_thr[j-1]-new_m[j,])*sqrt(new_cov[j,j])),
                             -lambda_fun_above((rel_thr[j-1]-new_m[j,])*sqrt(new_cov[j,j]))/
                               (1/(pnorm(-rel_thr[j-1])*rel_w[j-1])-1))
      lambda_alpha <- ifelse(rel_status[j-1]==1,
                             lambda_fun_above((rel_thr[j-1]-new_m[j,])*sqrt(new_cov[j,j])),
                             -lambda_fun_above((rel_thr[j-1]-new_m[j,])*sqrt(new_cov[j,j]))/
                               (1/(pnorm(-(rel_thr[j-1]-new_m[j,])*sqrt(new_cov[j,j]))*rel_w[j-1])-1))
      #lambda_alpha <- ifelse(rel_status[j-1]==1,
      #                       lambda_fun_above((rel_thr[j-1]-new_m[j,])*sqrt(new_cov[j,j])),
      #                       lambda_fun_below((rel_thr[j-1]-new_m[j,])*sqrt(new_cov[j,j])))
      updated_m <-lambda_alpha *sqrt(new_cov[j,j])+new_m[j,]
      delta_alpha =  delta_fun(lambda_alpha,(rel_thr[j-1]-new_m[j,])*sqrt(new_cov[j,j]))
      updated_var <-new_cov[j,j]*(1-delta_alpha)
      
      new_m  = new_m[-j,] + new_cov[-j,j]%*%solve(new_cov[j,j])%*%(updated_m-new_m[j,1])
      new_cov = new_cov[-j,-j,drop=F] -
        new_cov[-j,j]%*%(solve(new_cov[j,j])-solve(new_cov[j,j])%*%updated_var%*%solve(new_cov[j,j]))%*%new_cov[j,-j]
    }
    c(postM=new_m,postVar=new_cov)}
}




PA_thr <- function(df, stat_cols=which(grepl("stat",colnames(df))) ,age_cols=which(grepl("age",colnames(df))) , rel_mat,h2, prev, env_cor, aoo,ped,norm_weights=T){ 
  data.frame(Reduce('rbind',lapply(1:nrow(df), function(y){
    age= df[y,age_cols]
    status= df[y,stat_cols]
    
    w <- sapply(age,aoo)
    
    #if(norm_weights) w <- w/mean(w)
    #w[!is.na(status)&status==1] <- 1 
    
    #thr <- sapply(w, function(x) qnorm(1-prev*x))
    thr <- rep(qnorm(1-prev),length(status))
    
    
    
    covmatrix= rel_mat*h2
    diag(covmatrix) <- 1
    covmatrix[1,1] = h2
    
    
    pa_thr_hier(status,rel_thr = thr,rel_w = w,covmat = covmatrix)})))
}


PA_thr_2022 <- function(df, stat_cols=which(grepl("stat",colnames(df))) ,age_cols=which(grepl("age",colnames(df))) , rel_mat,h2, prev, env_cor, aoo,ped,norm_weights=T){ 
  data.frame(Reduce('rbind',lapply(1:nrow(df), function(y){
    age= df[y,age_cols]
    status= df[y,stat_cols]
    
    w <- sapply(age,aoo)
    
    #if(norm_weights) w <- w/mean(w)
    #w[!is.na(status)&status==1] <- 1 
    
    #thr <- sapply(w, function(x) qnorm(1-prev*x))
    thr <- rep(qnorm(1-prev),length(status))
    
    
    
    covmatrix= rel_mat*h2
    diag(covmatrix) <- 1
    covmatrix[1,1] = h2
    
    
    pa_thr_hier_2022(status,rel_thr = thr,rel_w = w,covmat = covmatrix)})))
}

PA_thr_2022_2 <- function(df, stat_cols=which(grepl("stat",colnames(df))) ,age_cols=which(grepl("age",colnames(df))) , rel_mat,h2, prev, env_cor, aoo,ped,norm_weights=T){ 
  data.frame(Reduce('rbind',lapply(1:nrow(df), function(y){
    age= df[y,age_cols]
    status= df[y,stat_cols]
    
    w <- sapply(age,aoo)
    
    #if(norm_weights) w <- w/mean(w)
    #w[!is.na(status)&status==1] <- 1 
    
    #thr <- sapply(w, function(x) qnorm(1-prev*x))
    thr <- rep(qnorm(1-prev),length(status))
    
    
    
    covmatrix= rel_mat*h2
    diag(covmatrix) <- 1
    covmatrix[1,1] = h2
    
    
    pa_thr_hier_2022_2(status,rel_thr = thr,rel_w = w,covmat = covmatrix)})))
}

# Our implementation of So et al AJHG2011
PA <- function(df, stat_cols=which(grepl("stat",colnames(df))) ,age_cols=which(grepl("age",colnames(df))) , rel_mat,h2, prev, env_cor, aoo,ped){ 
  data.frame(Reduce('rbind',lapply(1:nrow(df), function(y){
    #  age= df[y,age_cols]
    status= df[y,stat_cols]
    
    #    w <- sapply(age,aoo)
    #   w[!is.na(status)&status==1] <- 1 
    
    thr <- rep(qnorm(1-prev),length(status))
    
    covmatrix= rel_mat*h2
    diag(covmatrix) <- 1
    covmatrix[1,1] = h2
    
    
    pa_thr_hier(status,rel_thr = thr,rel_w=rep(1,length(status)),covmat = covmatrix)})))
}

# devtools::install_github("DesmondCampbell/diseaseRiskPredictor")
library(diseaseRiskPredictor)


gibbs <-  function(df, stat_cols=which(grepl("stat",colnames(df))) ,age_cols=which(grepl("age",colnames(df))) , rel_mat,h2, prev, env_cor, aoo,ped,iter=30000){ 
  data.frame(Reduce('rbind',lapply(1:nrow(df), function(y){
    age= df[y,age_cols]
    status= unlist(df[y,stat_cols])
    
    dfPed=ped  
    dfPed$age=unlist(c(0,age))
    thr=qnorm(1-prev)
    dfPed$thr=thr
    dfPed$expressedProportionOfLifetimeRisk = sapply(dfPed$age, aoo) 
    dfPed$affected[2:nrow(dfPed)] = status
    
    nonna <- unlist(c(T, !is.na(dfPed$affected[2:nrow(dfPed)])))
    
    dfPed<- dfPed[nonna,]
    
    
    covmatrix= rel_mat*h2
    diag(covmatrix) <- 1
    covmatrix[1,1] = h2
    
    
    mPriorCov = covmatrix[nonna,nonna]
    status <- status[nonna[-1]]
    rownames(mPriorCov) <- colnames(mPriorCov) <- paste0("liability_Tot_",dfPed$id)
    # mPriorMean =   updated_m <- t(t(c(0,sapply(1:length(status), function(x){ # E(l_j|Y_j,T_j)
    #   ifelse(status[x]==1, 
    #          lambda_fun_above(thr),
    #          lambda_fun_below(thr)
    #   )}) )))
    mPriorMean =   updated_m <- t(t(c(0,rep(0,length(status) ))))
    mPriorMean[is.na(mPriorMean),] <-0
    rownames(mPriorMean) <- paste0("liability_Tot_",dfPed$id)
    lConditioningOnRest <- diseaseRiskPredictor:::fnCalcConditioningOnRestOfVariables( mPriorCov )
    lmRegCoeff <- lConditioningOnRest$lmRegCoeff
    vPostStdDevB <- lConditioningOnRest$vPostStdDevB
    mean(diseaseRiskPredictor:::fnGibbsSampler( dfPed, mPriorMean, lmRegCoeff, vPostStdDevB, 100, nofIter = iter, F )[,1])
    
  })))}
