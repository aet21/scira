#' @title
#' Limma differential gene expression analysis for bulk tissue data
#'
#' @aliases LimmaFn
#'
#' @description
#' This function runs limma analysis and outputs ranked list of DEGs for each contrast of interest
#'
#' @param pheno.v
#' A phenotype vector. In our case, this will label the tissue-type in the bulk GTEX dataset
#'
#' @param data.m
#' The normalized gene expression data with rownames annotated to Entrez gene IDs (ensure these are unique) and columns labeling the bulk tissue samples.
#'
#' @return A list with two elements.
#'
#' @return top
#' A list itself with as many entries as their contrasts, and each entry being a matrix of top-ranked DEGs.
#'
#' @return cont
#' The matrix of contrasts.
#'
#' @examples
#' pheno.v <- c(rep(1,50),rep(2,50));
#' data.m <- matrix(rnorm(10000,0,1),nrow=1000,ncol=100);
#' data.m[1:100,1:50] <- matrix(rnorm(5000,2,1),nrow=100,ncol=50);
#' out.o <- LimmaFn(pheno.v,data.m);
#'
#' @import limma 
#'
#' @export
#' 
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
