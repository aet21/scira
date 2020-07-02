#' @title
#' Builds tissue-specific regulatory network
#'
#' @aliases sciraSelReg
#'
#' @description
#' This function derives the tissue-specific regulatory network
#'
#' @param regnet.o
#' The output from 'sciraInfReg'
#'
#' @param tissue.v
#' A vector listing the tissue-type of each sample. Note, this vector must match the columns of the data matrix used as input to the 'sciraInfReg' function.
#' 
#' @param toi
#' The tissue-type of interest.
#'
#' @param cft.v
#' A vector specifying the tissue types that should be treated as confounders. Often this will include blood since immune cells infiltrate most tissues
#'
#' @param degth.v
#' A vector specifying the adjusted P-value thresholds to use. The length of this vector equals 1+number of entries in 'cft.v'. The first P-value threshold refers to the comparison of the tissue type of interest to all other tissues together, whereas subsequent P-values are for the comparison of the tissue type of interest to each of the confounding tissue types.
#' 
#' @param lfcth.v 
#' A corresponding threshold on the log2 fold changes for each of the comparison. The length of this vector equals 1+number of entries in 'cft.v'.
#'
#' 
#' @return A list of 3 elements
#'
#' @return netTOI
#' The tissue-specific regulon matrix with columns labeling the transcription factors and rows labeling the target genes, with entries either 1 (positive interaction), -1 (inhibitory interaction) or 0 (no association).
#'
#' @return sumnet
#' A summary of the tissue-specific regulatory network listing number of targets and numbers of positive and inhibitory interactions
#'
#' @return top
#' A list of matrices, each matrix listing the top-ranked DEGs for each comparison of interest. First entry corresponds to the comparison of the tissue type of interest against all other tissue-types, the remaining entries correspond to the comparison of the tissue type of interest to the confounding tissue types.
#'
#' 
#' @examples
#' pheno.v <- c(rep(1,50),rep(2,50));
#' data.m <- matrix(rnorm(10000,0,1),nrow=1000,ncol=100);
#' data.m[1:100,1:50] <- matrix(rnorm(5000,2,1),nrow=100,ncol=50);
#' out.o <- LimmaFn(pheno.v,data.m);
#'
#' 
#' @export
#' 

sciraSelReg <- function(regnet.o,tissue.v,toi,cft.v=NULL,degth.v=rep(0.05,length(cft.v)),lfcth.v=c(1,rep(log2(1.5),length(cft.v)))){

   tt.v <- levels(factor(tissue.v));
   if( length(intersect(toi,tt.v))==0){
       print("Your tissue of interest is not in tissue.v, or you have mispelled toi");
       stop;
   }
   exp.m <- regnet.o$exp;
   gtexNETf.m <- regnet.o$regnet;
   ### now which TFs are overexpressed in toi?
   topTOI.lm <- list();
   ### compare toi to all other tissues
   pheno.v <- rep(0,ncol(exp.m));
   pheno.v[which(tissue.v==toi)] <- 1;
   lim.o <- LimmaFn(pheno.v,exp.m);
   topTOI.lm[[1]] <- lim.o$top[[1]];

   ### now compare toi to other tissues (in order to avoid confounding by immune or stromal cell infiltrates)
   if(!is.null(cft.v)){
    ti <- 2;
    for(t in cft.v){
     sel.idx <- which(tissue.v %in% c(toi,t));
     tmp.v <- tissue.v[sel.idx];
     tmpPH.v <- rep(0,length(sel.idx));
     tmpPH.v[which(tmp.v==toi)] <- 1;
     lim.o <- LimmaFn(tmpPH.v,exp.m[,sel.idx]);
     topTOI.lm[[ti]] <- lim.o$top[[1]];
     ti <- ti+1;
    }
   }

   ### now find tissue-specific TFs
   toiTF.lv <- list();
   for(i in 1:length(topTOI.lm)){
    statTF.m <- topTOI.lm[[i]][match(colnames(gtexNETf.m),rownames(topTOI.lm[[i]])),c("logFC","t","P.Value","adj.P.Val")];
    toiTF.idx <- intersect(which(statTF.m[,4] < degth.v[i]),which(statTF.m[,1]>lfcth.v[i]));
    toiTF.lv[[i]] <- rownames(statTF.m[toiTF.idx,]);
   }

   toiTF.v <- toiTF.lv[[1]];
   if(length(toiTF.lv)>1){
     for(i in 2:length(toiTF.lv)){
       toiTF.v <- intersect(toiTF.v,toiTF.lv[[i]]);
     }
   }
   map.idx <- match(toiTF.v,colnames(gtexNETf.m));
   #### tissue-specific regulatory network is:
   netTOI.m <- gtexNETf.m[,map.idx];

   distNet.m <- matrix(nrow=ncol(netTOI.m),ncol=3);
   rownames(distNet.m) <- colnames(netTOI.m);
   colnames(distNet.m) <- c("nTGTS","Act","Rep");
   distNet.m[,1] <- apply(abs(netTOI.m),2,sum);
   for(c in 1:ncol(netTOI.m)){
    act.idx <- which(netTOI.m[,c]==1)
    inact.idx <- which(netTOI.m[,c]==-1)
    distNet.m[c,2:3] <- c(length(act.idx),length(inact.idx));
   }

   return(list(netTOI=netTOI.m,sumnet=distNet.m,top=topTOI.lm));

} ### EOF


### Auxilliary function
#' @import limma
### LimmaFn
LimmaFn <- function(pheno.v,data.m){

### construct model matrix
sampletype.f <- as.factor(pheno.v);
design.sample <- model.matrix(~0 + sampletype.f);
colnames(design.sample) <- levels(sampletype.f);
sampletypes.v <- levels(sampletype.f);

### do linear model fit
lmf.o <- lmFit(data.m,design.sample);

### construct contrast matrix
ntypes <- length(levels(sampletype.f));
ncomp <- 0.5*ntypes*(ntypes-1);
cont.m <- matrix(0,nrow=ncol(design.sample),ncol=ncomp);
tmp.v <- vector();
c <- 1;
for(i1 in 1:(ntypes-1)){
 for(i2 in (i1+1):ntypes){
   cont.m[i1,c] <- -1;
   cont.m[i2,c] <- 1;
   tmp.v[c] <- paste(sampletypes.v[i2],"--",sampletypes.v[i1],sep="");
   c <- c+1;
 }
}
rownames(cont.m) <- sampletypes.v; # sampletype.v determined separately
colnames(cont.m) <- tmp.v;

### do linear model to contrasts
lmf2.o <- contrasts.fit(lmf.o,cont.m);

### empirical Bayesian estimation of differentially expressed genes (DEGs)
bay.o <- eBayes(lmf2.o);

### build ranked list of DEGs for each comparison
top.lm <- list();
for(c in 1:ncol(cont.m)){
top.lm[[c]] <- topTable(bay.o,coef=c,adjust="fdr",number=nrow(data.m));
}

return(list(top=top.lm,cont=cont.m));

} ## end of limma function

