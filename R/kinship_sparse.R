#' kinship sparse 
#' 
#' Compute a sparse kinship matrix from pedigree data, a pedigree object or a pedigreeList object
#' @param id either a pedigree object, pedigreeList object, or a vector of subject identifiers.  Subject identifiers may be numeric or character. 
#' @param dadid for each subject, the identifier of the biological father. This is only used if \code{id} is a vector.
#' @param momid for each subject, the identifier of the biological mother. This is only used if \code{id} is a vector.
#' @param sex vector of sex values coded as 1=male, 2=female
#' @param chrtype chromosome type.  The currently supported types are "autosome" and "X" or "x". 
#' @param ... extra arguments to pass forward to internal functions.
#' @details This is very minor modification of the \code{kinship2::kinship()} function. The only difference is that the kinship matrix returned by this function is always a sparse matrix. This can be usefull when working with large interconnected pedigrees, which cannot be cut into smaller pedigrees by the \code{kinship2::makefamid} function.  
#' 
#' @references Sinnwell J, P, Therneau T, M, Schaid D, J: The kinship2 R Package for Pedigree Data. Hum Hered 2014;78:91-93. doi: 10.1159/000363105
#' 
#' @examples 
#' library(kinship2)
#' data(minnbreast)
#' mb <- minnbreast
#' bpeds <- pedigree(mb$id, mb$fatherid, mb$motherid, 
#'                   mb$sex, affected=mb$proband, 
#' 		     famid=mb$famid)
#'
#' identical(as.matrix(kinship_sparse(bpeds[1])),kinship(bpeds[1]))  
#'
#' mb2 <- minnbreast[minnbreast$famid==4,]
#' 
#' identical(
#'    kinship(mb2$id,mb2$fatherid,mb2$motherid,mb2$sex),
#'    as.matrix(kinship_sparse(mb2$id, mb2$fatherid, mb2$motherid, mb2$sex))
#' )
#' 
#' # Note that these are not very efficient:
#' mb3 <- minnbreast[minnbreast$famid %in% 4:50,]
#' system.time(kinship_sparse(mb3$id, mb3$fatherid, mb3$motherid, mb3$sex))
#' system.time(kinship_sparse(pedigree(mb3$id, mb3$fatherid, mb3$motherid, mb3$sex)))
#' 
#' # Splitting into families and computing per per family is much faster: 
#' system.time(makefamid(mb3$id, mb3$fatherid, mb3$motherid))
#' system.time(kinship(pedigree(mb3$id, mb3$fatherid, mb3$motherid, mb3$sex, famid=mb3$famid)))
#' @importFrom methods as 
#' @importFrom methods is 
#' @export
kinship_sparse <- function(id, ...) {
  UseMethod('kinship_sparse')
}

#' @rdname kinship_sparse
#' @importFrom Matrix Diagonal 
#' @export
kinship_sparse.default <- function(id, dadid, momid, sex, chrtype="autosome", ...) {
  chrtype <- match.arg(casefold(chrtype), c("autosome", "x"))
  if (any(duplicated(id))) stop("All id values must be unique")
  n <- length(id)
  pdepth <- kindepth(id, dadid, momid)
  if (chrtype == "autosome") {
    if (n==1) 
      return(matrix(.5,1,1, dimnames=list(id, id)))
    
    kmat <- Diagonal(n+1,c(rep(.5, n), 0))  #founders
    
    mrow <- match(momid, id, nomatch=n+1) #row number of the mother
    drow <- match(dadid, id, nomatch=n+1) #row number of the dad 
    ## When all unrelateds, pdepth all=0. 
    ## Put c(1,) to make guard from iter 1:0
    for (depth in 1:max(c(1,pdepth))) {
      for (j in  (1:n)[pdepth==depth]) {
        kmat[,j] <-kmat[j,] <- (kmat[mrow[j],]  + kmat[drow[j],]) /2
        kmat[j,j] <- (1 + kmat[mrow[j], drow[j]]) /2
      }
    }
  }
  else if (chrtype == "x") {
    if (missing(sex) || length(sex) !=n) 
      stop("invalid sex vector")
    #1 = female, 2=male
    if (n==1) 
      return(matrix(ifelse(sex>2,sex/2,NA), 1,1, dimnames=list(id, id)))
    
    # kmat <- diag(c((3-sex)/2, 0)) #founders
    kmat <- Diagonal(length(sex)+1,c((3-sex)/2, 0))  #1 for males, 1/2 for females
    mrow <- match(momid, id, nomatch=n+1) #row number of the mother
    drow <- match(dadid, id, nomatch=n+1) #row number of the dad 
    
    for (depth in 1:max(c(1,pdepth))) {
      for (j in (1:n)[pdepth==depth]) {
        if (sex[j] ==1) {
          kmat[,j] <- kmat[j,] <- kmat[mrow[j],]
          kmat[j,j]<- 1
        } 
        else if(sex[j] == 2) {
          kmat[,j] <-kmat[j,] <- (kmat[mrow[j],]  + kmat[drow[j],]) /2
          kmat[j,j] <- (1 + kmat[mrow[j], drow[j]]) /2
        } 
        else {
          kmat[,j] <-kmat[j,] <- NA
          kmat[j,j] <- NA 
        }
      }
    }
  }
  kmat <- kmat[1:n,1:n]
  dimnames(kmat) <- list(id, id)
  kmat
}

#' @rdname kinship_sparse
#' @export
kinship_sparse.pedigree <- function(id, chrtype="autosome", ...) {
  chrtype <- match.arg(casefold(chrtype), c("autosome", "x"))
  if (any(duplicated(id$id))) stop("All id values must be unique")
  n <- length(id$id)
  pdepth <- kindepth(id)
  
  # Are there any MZ twins to worry about?
  havemz <- FALSE
  if (!is.null(id$relation) && any(id$relation$code=="MZ twin")) {
    havemz <- TRUE
    ## Doc: MakeMZIndex
    temp <- which(id$relation$code=="MZ twin")
    ## drop=FALSE added in case only one MZ twin set
    mzmat <- as.matrix(id$relation[,c("indx1", "indx2")])[temp,,drop=FALSE]
    mzgrp <- 1:max(mzmat) #everyone starts in their own group
    ## The loop below will take k-1 iterations for a set labeled as
    ##   (k-1):k, ..., 4:3, 3:2, 2:1;  this is the worst case.
    while(1) {
      if (all(mzgrp[mzmat[,1]] == mzgrp[mzmat[,2]])) break
      for (i in 1:nrow(mzmat)) 
        mzgrp[mzmat[i,1]] <- mzgrp[mzmat[i,2]] <- min(mzgrp[mzmat[i,]])
    }
    mzindex <- cbind(unlist(tapply(mzmat, mzgrp[mzmat], function(x) {
      z <- unique(x)
      rep(z, length(z))})),
      unlist(tapply(mzmat, mzgrp[mzmat], function(x) {
        z <- unique(x)
        rep(z, each=length(z))})))
    mzindex <- mzindex[mzindex[,1] != mzindex[,2],]
  }
  
  if (chrtype == "autosome") {
    if (n==1) 
      return(matrix(.5,1,1, dimnames=list(id$id, id$id)))
    
    kmat <- Diagonal(n+1,x=c(rep(.5,n),0)) #founders
    mrow <- ifelse(id$mindex ==0, n+1, id$mindex)
    drow <- ifelse(id$findex ==0, n+1, id$findex)
    
    for (depth in 1:max(pdepth)) {
      indx <- which(pdepth == depth)
      kmat[indx,] <- (kmat[mrow[indx],] + kmat[drow[indx],]) /2
      kmat[,indx] <- (kmat[,mrow[indx]] + kmat[,drow[indx]]) /2
      for (j in indx) kmat[j,j] <- (1 + kmat[mrow[j], drow[j]])/2
      if (havemz) kmat[mzindex] <- (diag(kmat))[mzindex[,1]]
    }
  }
  else if (chrtype == "x") {
    sex <- as.numeric(id$sex) # 1 = female, 2=male
    if (n==1) 
      return(matrix(sex/2, 1,1, dimnames=list(id$id, id$id)))
    
    kmat <- Diagonal(length(sex)+1,c((3-sex)/2, 0))  #1 for males, 1/2 for females
    mrow <- ifelse(id$mindex ==0, n+1, id$mindex)
    drow <- ifelse(id$findex ==0, n+1, id$findex)
    
    for (depth in 1:max(pdepth)) {
      for (j in (1:n)[pdepth==depth]) {
        if (sex[j] ==1) {
          kmat[,j] <- kmat[j,] <- kmat[mrow[j],]
          kmat[j,j]<- 1
        }
        else if(sex[j]==2) {
          kmat[,j] <-kmat[j,] <- (kmat[mrow[j],]  + kmat[drow[j],]) /2
          kmat[j,j] <- (1 + kmat[drow[j],mrow[j]]) /2
        } else {
          kmat[,j] <-kmat[j,] <- NA
          kmat[j,j] <- NA
        }
        if (havemz) kmat[mzindex] <- (diag(kmat))[mzindex[,1]]
      }
    }
  }
  kmat <- kmat[1:n,1:n]
  dimnames(kmat) <- list(id$id, id$id)
  kmat
}

#' @rdname kinship_sparse
#' @importFrom Matrix forceSymmetric 
#' @importFrom Matrix bdiag 
#' @export
kinship_sparse.pedigreeList <- function(id, chrtype="autosome", ...) {
  famlist <- unique(id$famid)
  nfam <- length(famlist)
  matlist <- vector("list", nfam)
  idlist  <- vector("list", nfam) #the possibly reorderd list of id values
  
  for (i in 1:length(famlist)) {
    tped <- id[i]  #pedigree for this family
    temp <- try(kinship(tped, chrtype=chrtype, ...), silent=TRUE)
    if (class(temp)=="try-error") 
      stop(paste("In family", famlist[i], ":", temp))
    else matlist[[i]] <- as(forceSymmetric(temp), "dsCMatrix")
    idlist[[i]] <- tped$id
  }
  
  result <- bdiag(matlist)
  if (any(duplicated(id$id))) 
    temp <-paste(rep(famlist, sapply(idlist, length)),
                 unlist(idlist), sep='/') 
  else temp <- unlist(idlist)
  
  dimnames(result) <- list(temp, temp)
  result
}

#' kinship sparse path  
#' 
#' Compute a sparse kinship matrix from a pedigree object using path counting
#' @param ped a pedigree object
#' @details This is an alternative way of computing the kinhips matrix which is usually slower, but less memory intense.  
#'
#' @examples 
#' library(kinship2)
#' data(minnbreast)
#' mb <- minnbreast[minnbreast$famid %in% 4:5,]
#' system.time(k_path <- kinship_path(
#'      pedigree(mb$id, mb$fatherid, mb$motherid, mb$sex)))
#' system.time(k_sparse <- kinship_sparse(
#'      pedigree(mb$id, mb$fatherid, mb$motherid, mb$sex)))
#' system.time(k_per_fam <- kinship(
#'      pedigree(mb$id, mb$fatherid, mb$motherid, mb$sex, famid=mb$famid)))
#' identical(k_path,k_sparse) 
#' val1=!k_path[order(rownames(k_path)),order(rownames(k_path))]
#' val2=k_per_fam[order(rownames(k_per_fam)),order(rownames(k_per_fam))]
#' sum(val1==val2)
#' @rdname kinship_sparse_path
#' @importFrom igraph all_shortest_paths
#' @importFrom igraph vertex_attr
#' @importFrom utils tail
#' @importFrom Matrix sparseMatrix
#' @import kinship2
#' @export
kinship_path <- function(ped){
  ids <- ped$id
  # First convert the pedigree to a graph: 
  ped_graph <- ped2graph(ped)
  # We create a list of ancendents by finding all incomming paths: 
  anc <- Reduce('rbind.data.frame',lapply(vertex_attr(ped_graph)[[1]],function(i){
    p <- all_shortest_paths(ped_graph,from = i, mode="in")$res  
    if(length(p)!=0)
      cbind.data.frame(id1=as.numeric(i),id2=as.numeric(names(sapply(p,tail,n=1))),gen=sapply(p,length)-1,path=I(lapply(p,function(x) as.numeric(names(x)))))
  } 
  ))
  anc <- anc[!anc$id1==anc$id2,]
  # decendents are just the opposite: 
  dec <- setNames(anc,c("id2","id1","gen","path"))
  # now we merge these to find the descendents of the ascencendents 
  other <- merge(data.table(anc),data.table(dec),by.x = "id2",by.y = "id1",allow.cartesian = T)
  anc$gen.x =0
  dec$gen.y =0
  names(anc)[3] <- "gen.y"
  names(dec)[3] <- "gen.x"
  other$id2 <- NULL
  names(other)[4] <- "id2"
  
  # remove self: 
  other <- other[!id1==id2,]
  # remove dec: 
  other <- other[!paste(id1,id2,sep = "_") %in% paste(dec$id1,dec$id2,sep = "_"),]
  other <- other[!paste(id1,id2,sep = "_") %in% paste(dec$id2,dec$id1,sep = "_"),]
  # set data.table variables to make R checker happy (recomended solution by data.table team)
  path.y <- shortest <- contained <- NULL
  # concatenate path  
  other[,path:=Map(c,other$path.x,lapply(path.y,function(x) rev(x)[-1]))]
  # find the shortest paths (there are often 2) between each pair  
  setkey(other,id2)
  setkey(other,id1)
  # set data.table function to make R checker happy (recomended solution by data.table team)
  . <- NULL
  other[,shortest:=!other[,(gen.x+gen.y)>min(gen.x+gen.y),.(id1,id2)]$V1]
  other_shortest <- other[shortest==T]
  other_shortest <- other_shortest[,.(short_path=(list(path))),.(id1,id2)]
  # merge all "other"-paths with their corresponding shortest path: 
  other <- other_shortest[other,on=c("id1","id2")] 
  
  # find those that are not the shortest, but are not contained in the shortest:
  ## this is necessary to catch e.g. if someone is both siblings and cousins
  other_not_shortest <- other[shortest==F]
  other_not_shortest[,contained:=sapply(1:nrow(other_not_shortest),function(z) lapply(other_not_shortest[z]$short_path,function(x) any(unlist(lapply(x, function(y) all(y %in% other_not_shortest[z]$path[[1]]))))))]
  other_not_shortest[,contained:=unlist(contained)]
  
  other <- other[shortest==T]
  
  # Paths not contained by other paths
  non_contained <- other_not_shortest[contained==F]
  if(nrow(non_contained)>0){
    setkey(non_contained,id2)
    setkey(non_contained,id1)
    non_contained <- non_contained[!non_contained[,(gen.x+gen.y)>min(gen.x+gen.y),.(id1,id2)]$V1]
    other <- rbind.data.frame(other,non_contained[,-c("contained")])
  }
  # set data.table variables to make R checker happy (recomended solution by data.table team)
  id1 <- id2 <- gen.x <- gen.y <- path <- r <- NULL
  # combine anc, dec and other 
  rel <- data.table(rbind.data.frame(data.table(anc)[,.(id1,id2,gen.x,gen.y,path)],
                                     data.table(dec)[,.(id1,id2,gen.x,gen.y,path)],
                                     other[,.(id1,id2,gen.x,gen.y,path)]))
  
  # computed relatedness per path 
  rel[,r:= 0.5^(gen.x)*0.5^(gen.y)] 
  
  # sum all paths per relative pair: 
  k <- rel[,.(k=sum(r)/2,gen.x=mean(gen.x),gen.y=mean(gen.y)),.(id1,id2)]
  k <- rbind.data.frame(k,data.frame(id1=ids,id2=ids,k=.5,gen.x=0,gen.y=0))
  k$ind1=match(k$id1,ids)
  k$ind2=match(k$id2,ids)
  k_mat <- sparseMatrix(k$ind1,k$ind2,x = k$k,dimnames = list(ids,ids))
  k_mat}


