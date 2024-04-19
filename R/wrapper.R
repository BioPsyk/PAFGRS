#' FGRS wrapper 
#' 
#' Compute (PA-)FGRS from for a large a dataset
#' @param proband_ids vector of sample ids indicating which individuals FGRS should be estimated for. 
#' @param K kinship matrix provided either as a \code{matrix}, a \code{dsCMatrix} or a \code{data.frame} with column names \code{i}, \code{j} and \code{x}.
#' @param pheno \code{data.frame} containing the phenotype information on the relatives. Must contain columns called \code{id} and \code{aff}. 
#' @param method character indicating which FGRS method to use. Default is \code{method="PAFGRS"}. 
#' @param thr numeric vector wiht threshold value used for each relative in \code{pheno}.
#' @param w numeric vector of proportion of risk experienced by each relative in \code{pheno}.
#' @param h2 numeric heritability estimate of the phenotype (liability scale). 
#' @param env_cor_s numeric indicating the down-weighting of siblings (only used for Ohlsson-Kendler FGRS). Default =1 (no correction).
#' @param env_cor_f numeric indicating the down.weighting of fathers (only used for Ohlsson-Kendler FGRS). Default =1 (no correction).
#' @param env_cor_m numeric indicating the down-weighting of mothers (only used for Ohlsson-Kendler FGRS). Default =1 (no correction).
#' @param sib_mat optional \code{data.frame}  with colunm names \code{i} and \code{j} giving the \code{id}s of individuals that are siblings. 
#' @param father_mat optional \code{data.frame}  with colunm names \code{i} and \code{j} giving the \code{id}s of individuals that are father-child relations.  
#' @param mother_mat optional \code{data.frame}  with colunm names \code{i} and \code{j} giving the \code{id}s of individuals that are mother-child relations. 
#' @return postM posterior mean liability 
#' @return postVar posterior variance
#' @examples 
#' pa_fgrs(c(0,1),qnorm(.9),covmat = matrix(c(.5,.25,.25,.25,1,.25,.25,.25,1),3))
#' @export
FGRS_wrapper <- function(proband_ids,K,pheno,method="PAFGRS",thr=NULL,w=NULL,h2=NULL,env_cor_s=1,env_cor_f=1,env_cor_m=1,sib_mat=NULL,father_mat=NULL,mother_mat=NULL){
  if(is.numeric(proband_ids)) proband_ids <- as.integer(proband_ids)
  if(class(K)[1]=="matrix") K <- as(K, "sparseMatrix")
  if(class(K)[1]=="dsCMatrix"){
    k_dt <- data.table(summary(K))
    k_dt[,i:=rownames(K)[i]]
    k_dt[,j:=rownames(K)[j]] 
    K <- k_dt 
    rm(k_dt)}
  if(class(K)[1]=="data.frame") K <- data.table(K) 
  if(is.integer(proband_ids)){
    K[,i:=as.integer(i)]
    K[,j:=as.integer(j)]}else{
      K <- rbind(K[!j<i],K[j<i,.(i=j,j=i,x)])[!duplicated(paste(i,j))]
    }
  if(!all(c("i","j","x") %in% colnames(K))) stop("K should contain columns 'i','j' and 'x'")
  if(K[,any(i==j)]) if(K[i==j,.N,x][order(-N)][1,x==1]) stop("most i==j entries in K are 1, did you provide relatedness matrix instead of kinship matrix?")
  if(any(!proband_ids %in% K[,c(unique(i),unique(j))])) warning("some proband_ids are not in K")
  K <-K[!i==j]
  if(K[,any(j<i)]){
    if(!K[,any(i<j)]) K <- K[,.(i=j,j=i,x)] # just flip i and j
    ij_repeats=K[K[j<i,.(i=j,j=i,x)],on=c("i","j"),nomatch=0]
    if(nrow(ij_repeats)>0){
      if(ij_repeats[,all(x==i.x)]) K <- rbind(K[i<j],K[i>j,.(i=j,j=i,x)])[,.(x=mean(x)),.(i,j)] else
        stop("conflicting entries in K")} else 
          K <- rbind(K[i<j],K[i>j,.(i=j,j=i,x)])}
  if(K[,any(x==1)]) stop("some i!=j entries in K have x==1, did you provide relatedness matrix instead of kinship matrix?")
  if(K[x==.5,.N]>K[x==.25,.N]) stop("more i!=j entries in K have x==0.5 than x==0.25, did you provide relatedness matrix instead of kinship matrix?")
  setkey(K,i,j)
  
  # sib_mat 
  if(!is.null(sib_mat)){
  if(class(sib_mat)[1]=="data.frame") sib_mat <- data.table(sib_mat) else
    if(!class(sib_mat)[1]=="data.table") stop("sib_mat should be either a data.frame or a data.table")
  setkey(sib_mat,i,j)}
  
  if(!is.null(thr)) pheno$thr = thr 
  if(!"thr" %in% colnames(pheno)) 
    stop("'thr' should be provided either as a column in 'pheno' or by the 'thr' argument")
  if(!is.null(w)) pheno$w = w
  if(!"w" %in% colnames(pheno)) 
    stop("'w' should be provided either as column in 'pheno' or by the 'w' argument")

  pheno <- data.table(pheno)
  if(!all(c("id","aff") %in% colnames(pheno))) stop("'pheno' should contain columns 'id' and 'aff'")
  setkey(pheno,id)
  if(pheno[,any(aff==1&w!=1)]) { pheno[aff==1,w:=1]
    message("setting w=1, for all relatives with aff==1")}
  #if(any(!proband_ids %in% pheno$id)) stop("some proband ids are not in pheno$id")
  if(any(!proband_ids %in% K[,c(unique(i),unique(j))])) warning("some proband_ids do not seem to have any relatives, returning postM=0")
  if(any( colnames(pheno) %in% c("fatherid","motherid","id_f","id_m","momid","dadid"))){
    par_id <- min(which( colnames(pheno) %in% c("fatherid","motherid","id_f","id_m","momid","dadid")))
    par_kin25 =K[setNames(pheno[,c(which(colnames(pheno)=="id"),par_id),with=F],c("i","j"))[i<j],on=.(i,j),mean(x==.25,na.rm=T)]
    if(abs(par_kin25-1)>0.05) 
      stop(paste(round(par_kin25*100,1), "% of", colnames(pheno)[par_id], "have kinship==0.25 with 'id'. Is something wrong?"))
  }
  if(method!="OK2") { 
    rels <- K[.(proband_ids),.(list(c(j))),i,nomatch=0]
    more_rels <- K[j %in% proband_ids,.(list(c(i))),j]
    rels=merge(rels,more_rels,by.x = "i",by.y = "j",all = T)
    rels=merge(rels,data.table(i=proband_ids,key="i"),all.y = T)
    rels_and_self=rels[,.(list(c(i,unlist(V1.y),unlist(V1.x)))),i] 
    }

  if(method=="OK") { 
    if(!is.null(h2)) stop("Ohlsson-Kendler method does not require h2 paramter")
    if(env_cor_s<1 & is.null(sib_mat)) stop("when 'env_cor_s<1' you must provide 'sib_mat'")
    if(env_cor_f<1 & is.null(father_mat)) stop("when 'env_cor_f<1' you must provide 'father_mat'")
    if(env_cor_m<1 & is.null(mother_mat)) stop("when 'env_cor_m<1' you must provide 'mother_mat'")
    pheno[aff==0,z:=trunc_norm_below(0,1,thr)[1]]
    pheno[aff==1,z:=trunc_norm_above(0,1,thr)[1]]
    }

  ghat <- sapply(rels_and_self[,V1], function(k){
    
    
    k_proband <- K[CJ(i=k,j=k,unique = T)[i<j],on=.(i,j),nomatch=0]
     
    pheno_k <- pheno[.(k)]
    k_proband[,i_ind:=match(i,sort(k))]
    k_proband[,j_ind:=match(j,sort(k))]
    if(method=="OK") 
      h2 <- 1
    if(length(k)>1) 
      covmat <- sparseMatrix(k_proband$i_ind,
                             k_proband$j_ind,
                             x=k_proband$x*2*h2,
                             symmetric = T,
                             dimnames=list(1:k_proband[,max(c(i_ind,j_ind))],1:k_proband[,max(c(i_ind,j_ind))]))[match(k,sort(k)),match(k,sort(k))] else
                               covmat <- matrix(h2)
  
    diag(covmat) <- 1  
    covmat[1,1] <- h2
      # 
      # print(k)
      # print(pheno_k$id)
      # print(match(k,sort(k)))
      # print(pheno_k$id[match(k,sort(k))])
      # print(dimnames(covmat))
      # 
    if(method=="OK") { 
      if(length(k)>1) {
        pheno_k$r = covmat[,1]/h2
        pheno_k$c=1
        
        if(!is.null(sib_mat)){ pheno_k$c = pheno_k$c*ifelse(pheno_k[,id] %in% sib_mat[i==pheno_k[1,id],j],env_cor_s,1)}
        if(!is.null(father_mat)){ 
          pheno_k$c = 
            pheno_k$c*ifelse(pheno_k[,id] %in% father_mat[i==pheno_k[1,id],j],
                             env_cor_f,1)}
        if(!is.null(mother_mat)){ 
          pheno_k$c = 
            pheno_k$c*ifelse(pheno_k[,id] %in% mother_mat[i==pheno_k[1,id],j],
                             env_cor_m,1)}
        out1 <- data.frame(pheno_k[-1,.(mean(z*w*r*c,na.rm=T),sum(r[!is.na(z*w*r*c)]))])} else
          out1 <- data.frame(0,0)
      
      } else  if(method=="PAFGRS") out1 <- pa_fgrs(rel_status = pheno_k$aff[-1],rel_thr =pheno_k$thr[-1],rel_w = pheno_k$w[-1],covmat = as.matrix(covmat))
    return(out1)  
  
  })
  if(method=="OK") {
    out <- data.table(id=rels_and_self[,i], unlist(t(ghat)[,1]),unlist(t(ghat)[,2]),n_rels=rels_and_self[,sapply(V1,length)-1])
    vz = var(pheno$z,na.rm = T)
    vs = var(out$V2,na.rm = T)
    out[,FGRS := out[,V2*vs/(vs+vz/V3)]]
    out = out[,.(id,FGRS,n_rels,s=V2,vs=vs,vz=vz,sum_r=V3)]

      } else if(method=="PAFGRS") {out=data.frame(id=rels_and_self[,i], t(ghat),n_rels=rels_and_self[,sapply(V1,length)-1])}
  if(method=="OK2"){
  K <- data.table(K)
  K[,i_ind:=match(i,pheno$id)]
  K[,j_ind:=match(j,pheno$id)]
  k_sparse_t <- sparseMatrix(i = K$i_ind,j = K$j_ind,x = K$x,dims = rep(length(pheno$id),2),dimnames =list(pheno$id,pheno$id) ,symmetric = T)
  
  out <-FGRS_kendler(pheno_t = pheno,k_sparse_t = k_sparse_t,sib_mat_t = sib_mat,father_mat_t = father_child_mat,mother_mat_t = mother_child_mat,env_cor_f = 0.5,env_cor_m = 0.5,env_cor_sib =0.5 )
  }
  return(data.frame(out))
}

FGRS_kendler = function(pheno_t,k_sparse_t,sib_mat_t,father_mat_t,mother_mat_t,env_cor_f,env_cor_m,env_cor_sib){
  pheno_t <- data.table(pheno_t)
  pheno_t[aff==0,z:=trunc_norm_below(0,1,thr)[1]]
  pheno_t[aff==1,z:=trunc_norm_above(0,1,thr)[1]]
  
  k_temp = k_sparse_t*2
  #k_temp =k_temp[order(as.numeric(row.names(k_temp))),order(as.numeric(row.names(k_temp)))] 
  diag(k_temp) = 0
  sum_r <- rowSums(k_temp)
  n_rel <- rowSums(k_temp>0)
  
  s_temp= sparseMatrix(i = sib_mat_t$i,j = sib_mat_t$j,x = 1-0.5,dims = dim(k_temp))
  f_temp= sparseMatrix(i = father_mat_t$i,j = father_mat_t$j,x = 1-0.5,dims = dim(k_temp))
  m_temp= sparseMatrix(i = mother_mat_t$i,j = mother_mat_t$j,x = 1-0.5,dims = dim(k_temp))
  
  k_temp = k_temp - (k_temp* s_temp)  - (k_temp* f_temp)  - (k_temp* m_temp)
  S <- k_temp %*%  (pheno_t$z * pheno_t$w)
  S <- S/n_rel
  vs = var(S[,1],na.rm = T)
  vz = var(pheno_t$z)
  FGRS = S*vs/(vs+vz/sum_r)
  FGRS_out <- data.table(id= as.numeric(rownames(FGRS) ),FGRS=FGRS[,1])
  FGRS_out
}
