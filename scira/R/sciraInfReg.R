#' @title
#' Infers Regulatory Network
#'
#' @aliases sciraInfReg
#'
#' @description
#' This function infers a general regulatory network given a bulk RNA-Seq dataset (e.g. GTEX)
#'
#' @param data.m
#' The normalized gene expression data with rownames annotated to Entrez gene IDs (ensure these are unique) and columns labeling the samples
#'
#' @param regEID.v
#' A vector listing the regulators/transcription factors, and annotated to Entrez gene IDs
#' 
#' @param sdth
#' A threshold on the standard deviation of gene expression in data.m, to select genes with standard deviation larger than sdth.
#'
#' @param sigth
#' The significance P-value threshold for the marginal correlations. If not specified it uses the Bonferroni threshold.
#'
#' @param spTH
#' This parameter controls the sparsity of the regulatory network. It is a threshold on the number of allowed marginal associations (expressed as a percentage). The default value is 0.01, i.e. if at the given value of sigth, the fraction of marginal correlations is larger than 0.01, we cap the fraction at 0.01 and only retain the top 1% of associations
#'
#' @param pcorth
#' A significance threshold on the partial correlation values. By default this is set to 0.2. User needs to select or estimate this for each dataset and choose an appropriate threshold designed to find the appropriate number of target genes per regulator.
#'
#' @param minNtgts
#' The minimum number of gene targets required for each transcription factor
#'
#' @param ncores
#' Number of parallel cores to use.
#' 
#' @return A list of two elements:
#'
#' @return regnet
#' The regulon matrix with columns labeling the transcription factors and rows labeling the target geens, with entries either 1 (positive interaction), -1 (inhibitory interaction) or 0 (no association).
#'
#' @return exp
#' The bulk tissue expression matrix.
#' 
#' @examples
#' pheno.v <- c(rep(1,50),rep(2,50));
#' data.m <- matrix(rnorm(10000,0,1),nrow=1000,ncol=100);
#' data.m[1:100,1:50] <- matrix(rnorm(5000,2,1),nrow=100,ncol=50);
#' out.o <- LimmaFn(pheno.v,data.m);
#'
#' @importFrom corpcor cor2pcor
#' @importFrom parallel mclapply
#' 
#' @export
#' 

sciraInfReg <- function(data.m,regEID.v,sdth=0.25,sigth=NULL,spTH=0.01,pcorth=0.2,minNtgts=10,ncores=4){

   ### remove genes with no or little variance
   sd.v <- apply(data.m,1,sd);
   selG.idx <- which(sd.v>sdth);
   exp.m <- data.m[selG.idx,];
   ### find representation of regulators in data, and define regulatees/targets
   tfEID.v <- regEID.v;
   repTF.v <- intersect(tfEID.v,rownames(exp.m));
   tgtsEID.v <- setdiff(rownames(exp.m),tfEID.v);
   match(repTF.v,rownames(exp.m)) -> mapTF.idx;
   match(tgtsEID.v,rownames(exp.m)) -> mapTGTS.idx;
   ### compute correlations and estimate P-values
   corNET.m <- cor(t(exp.m[mapTGTS.idx,]),t(exp.m[mapTF.idx,]));
   zNET.m <- 0.5*log( (1+corNET.m)/(1-corNET.m) );
   stdev <- 1/sqrt(ncol(exp.m)-3); 
   pvNET.m <- 2*pnorm(abs(zNET.m),0,stdev,lower.tail=FALSE)
   ### for each gene, now identify the TFs which are correlated univariately- for these then run multivariate regression
   if(is.null(sigth)){
    sigth <- 0.05/prod(dim(pvNET.m)); ### Bonferroni
   }
   binNET.m <- pvNET.m;
   binNET.m[pvNET.m < sigth] <- 1;
   binNET.m[pvNET.m >= sigth] <- 0;
   if(sum(binNET.m) > spTH*prod(dim(pvNET.m))){### e.g. if more than 1% cap at 1%
       topE <- floor(spTH*prod(dim(pvNET.m)));
       binNET.m <- genNetTopE(zNET.m,topE);
   }
    
   ### number of targets per tf
   ntgTF.v <- apply(binNET.m,2,sum)
   ### number of regulators per gene
   nregG.v <- apply(binNET.m,1,sum)

   ### select TFs with at least minNtgts targets
   selTF.idx <- which(ntgTF.v>=minNtgts);
   selbinNET.m <- binNET.m[,selTF.idx];
   selpvNET.m <- pvNET.m[,selTF.idx];
   selzNET.m <- zNET.m[,selTF.idx];
   selcorNET.m <- corNET.m[,selTF.idx];

   mapTG.idx <- match(rownames(selbinNET.m),rownames(exp.m));
   mapTF.idx <- match(colnames(selbinNET.m),rownames(exp.m));

   idx.l <- as.list(1:nrow(selbinNET.m));
   print("Computing Partial Correlations");
   pcor.l <- mclapply(idx.l,ComputePCOR,mapTG.idx,mapTF.idx,selbinNET.m,exp.m,mc.cores=ncores);

   pcorNET.m <- matrix(0,nrow=nrow(selbinNET.m),ncol=ncol(selbinNET.m));
   rownames(pcorNET.m) <- rownames(selbinNET.m);
   colnames(pcorNET.m) <- colnames(selbinNET.m);
   for(g in 1:length(pcor.l)){
    reg.idx <- which(selbinNET.m[g,]==1);
    if(length(reg.idx)>=2){
        pcorNET.m[g,reg.idx] <- pcor.l[[g]][1,-1];
    }
    else if (length(reg.idx)==1){
        pcorNET.m[g,reg.idx] <- selcorNET.m[g,reg.idx];
    }
    print(g);
   }

   gtexNET.m <- sign(pcorNET.m);
   gtexNET.m[abs(pcorNET.m) < pcorth] <- 0;

   sumnetTF.m <- matrix(nrow=4,ncol=ncol(gtexNET.m));
   rownames(sumnetTF.m) <- c("nTG","nUP","nDN","P");
   colnames(sumnetTF.m) <- colnames(gtexNET.m);
   sumnetTF.m[1,] <- apply(abs(gtexNET.m),2,sum);

   for(tf in 1:ncol(gtexNET.m)){
    sumnetTF.m[2,tf] <- length(which(gtexNET.m[,tf]==1));
    sumnetTF.m[3,tf] <- length(which(gtexNET.m[,tf]==-1));    
   }
   pv.v <- pbinom(apply(sumnetTF.m[2:3,],2,max),size=sumnetTF.m[1,],prob=0.5,lower.tail=FALSE);
   sumnetTF.m[4,] <- pv.v;

   gtexNETf.m <- gtexNET.m[,which(sumnetTF.m[1,]>= minNtgts)];
   return(list(regnet=gtexNETf.m,exp=exp.m));
} ## EOF

#### Auxillary internal functions
ComputePCOR <- function(idx,mapTG.idx,mapTF.idx,selbinNET.m,exp.m){
    g <- idx;
    reg.idx <- which(selbinNET.m[g,]==1);
    if(length(reg.idx)>=2){
        tmp.idx <- c(mapTG.idx[g],mapTF.idx[reg.idx]);
        cor.m <- cor(t(exp.m[tmp.idx,]));
        pcor.m <- cor2pcor(cor.m);
    }
    else{
        pcor.m <- NULL;
    }
    return(pcor.m);
} ## end of ComputePCOR function

genNetTopE <- function(z.m,topE){

    z.v <- as.vector(abs(z.m));
    tmp.s <- sort(z.v,decreasing=TRUE,index.return=TRUE);
    out.v <- rep(0,length(z.v));
    out.v[tmp.s$ix[1:topE]] <- 1;
    out.m <- matrix(out.v,nrow=nrow(z.m));
    rownames(out.m) <- rownames(z.m);
    colnames(out.m) <- colnames(z.m);
    return(out.m);

}
