#' @title
#' Estimates Regulatory Activity
#'
#' @aliases sciraEstRegAct
#'
#' @description
#' This function takes as input a set of regulons, and a matrix of RNA-Seq profiles (single cell or bulk) and
#' outputs a regulatory activity matrix over the transcription factors and samples
#' 
#' @param data
#' Either a vector with names labeling genes (Entrez gene IDs) or a matrix with rows labeling genes (Entrez gene IDs) and columns labeling samples (cells).
#'
#' @param regnet.m
#' The regulon matrix with rows labeling genes (Entrez IDs) and columns labeling the transcription factors (also Entrez gene IDs) and entries either 1 (positive), -1 (inhibitory) or 0 (no association). Typically, the output of either 'sciraSelReg' or 'sciraInfReg'
#' 
#' @param norm
#' Type of normalization to use if data is provided as a matrix. "z" stands for z-score transformation, "c" just row centres the data matrix.
#'
#' @param ncores
#' Number of parallel cores to use
#' 
#' 
#' @return A matrix of regulatory activity scores, with rows labeling the transcription factors and columns labeling samples (cells).
#'
#' 
#' @examples
#' pheno.v <- c(rep(1,50),rep(2,50));
#' data.m <- matrix(rnorm(10000,0,1),nrow=1000,ncol=100);
#' data.m[1:100,1:50] <- matrix(rnorm(5000,2,1),nrow=100,ncol=50);
#' out.o <- LimmaFn(pheno.v,data.m);
#'
#' @importFrom parallel mclapply
#'
#' @export
#' 

sciraEstRegAct <- function(data,regnet.m,norm=c("z","c"),ncores=4){

 if(is.vector(data)){
   common.v <- intersect(names(data),rownames(regnet.m));
   match(common.v,names(data)) -> map1.idx;
   match(common.v,rownames(regnet.m)) -> map2.idx;   
   actTF <- InferTFact(data[map1.idx],regnet.m[map2.idx,]);
   names(actTF) <- colnames(regnet.m);
 }
 else if (is.matrix(data)){
   common.v <- intersect(rownames(data),rownames(regnet.m));
   match(common.v,rownames(data)) -> map1.idx;
   match(common.v,rownames(regnet.m)) -> map2.idx;        
   ndata <- data[map1.idx,] - rowMeans(data[map1.idx,]);

   if(norm=="z"){
    sd.v <- apply(data[map1.idx,],1,sd);
    nz.idx <- which(sd.v>0);
    z.idx <- which(sd.v==0);
    ndata <- data[map1.idx,];       
    ndata[nz.idx,] <- (data[map1.idx[nz.idx],] - rowMeans(data[map1.idx[nz.idx],]))/sd.v[nz.idx];
    ndata[z.idx,] <- 0;
   }
   idx.l <- as.list(1:ncol(data));
   prl.o <- mclapply(idx.l,InferTFactPRL,ndata,regnet.m[map2.idx,],mc.cores=ncores);
   actTF <- matrix(unlist(prl.o),nrow=ncol(regnet.m),ncol=length(prl.o),byrow=FALSE)
   rownames(actTF) <- colnames(regnet.m);
   colnames(actTF) <- colnames(data);
 }

 return(actTF);
}

### Auxilliary Functions

InferTFact <- function(exp.v,regnet.m){
  act.v <- apply(regnet.m,2,function(tmp.v){lm.o <- lm(exp.v ~ tmp.v); act <- summary(lm.o)$coeff[2,3];return(act);})
  return(act.v);
} ### end of InferTFact function

InferTFactPRL <- function(idx,tmp.m,regnet.m){
  exp.v <- tmp.m[,idx];
  act.v <- apply(regnet.m,2,function(tmp.v){lm.o <- lm(exp.v ~ tmp.v); act <- summary(lm.o)$coeff[2,3];return(act);})
  return(act.v);
} ### end of InferTFactPRL function
