## ----chunk1, eval=T, echo=T---------------------------------------------------
library(parallel);
library(corpcor);
library(limma);
library(scira);
data(tfEID); ### loads in the list of transcriptions factors (annotated to Entrez gene IDs)
data(tissue); ### loads in tissue-type information of GTEX dataset (needed later)
#regnet.o <- sciraInfReg(gtex.m,tfEID.v,sdth=0.25,sigth=1e-6,spTH=0.01,pcorth=0.2,minNtgts=10,ncores=4); ### this infers the non-tissue specific regulatory network but we comment this out as the GTEX data matrix, labeled here as gtex.m, is large and the inference takes a few minutes. 

## ----chunk2, eval=T, echo=T---------------------------------------------------
#selreg.o <- sciraSelReg(regnet.o,tissue.v,toi="Liver",cft.v=c("Blood"),degth.v=rep(0.05,2),lfcth.v=c(log2(1.5),log2(1))); ### generates tissue-specific network
data(selreg);
print(names(selreg.o));

## ----chunk2b, eval=T, echo=T--------------------------------------------------
dim(selreg.o$netTOI);
head(selreg.o$sumnet);

## ----chunk2c, eval=T, echo=T, fig.cap = "The liver-specific regulatory network"----
netplot <- sciraPlotNet(selreg.o$netTOI);
netplot;

## ----chunk3, message = FALSE,warning=F,fig.dim = c(8,6),fig.align="center",fig.cap = "Validation of liver-specific regulons in ProteinAtlas dataset"----
###load in independent data
data(PrtAtlasExp)
tfaPrtAtl.m <- sciraEstRegAct(data = PrtAtlasExp.m, regnet = selreg.o$netTOI, norm = "z", ncores = 4)
avAct.v <- colMeans(tfaPrtAtl.m);
medact<-tapply(colMeans(tfaPrtAtl.m), colnames(tfaPrtAtl.m), median);
smedact<- sort(medact,decreasing = T)
pheno.v <- vector(length=ncol(tfaPrtAtl.m));
for(i in 1:length(smedact)){
  pheno.v[which(colnames(tfaPrtAtl.m)==names(smedact)[i])] <- i;
}
par(mar=c(11,4,1,1));
boxplot(avAct.v ~ pheno.v,names=names(smedact),las=2,ylab="AvAct(LiverTFs)",outline=FALSE,col=c(rep("grey",2),"red",rep("grey",length(smedact)-3)),cex.axis=0.8,xlab="");


## ----chunk4, message = FALSE,warning=F,fig.dim = c(6,5),fig.align="center",fig.cap = "Boxplot of the average activity level between liver and all other tissues"----

### now plot liver against all other tissues and compute P-value
pheno2.v <- pheno.v;
pheno2.v[pheno.v!=3] <- 4;
pheno2.v <- pheno2.v-2;
par(mar=c(4,4,2,1));
boxplot(avAct.v ~ pheno2.v,names=c("Liver","AllOther"),ylab="AvAct(LiverTFs)",outline=FALSE,col=c("red","grey"),cex.axis=1.25,xlab="");
pv <- t.test(avAct.v ~ pheno2.v,alt="gr")$p.value;
text(x=1.5,y=3,paste("P=",signif(pv,3),sep=""),font=2,cex=1.25);
nspg.v <- summary(factor(pheno2.v));
mtext(side=1,at=1:2,line=1.85,paste("n=",nspg.v,sep=""));


## ----chunk5, eval=T, echo=T---------------------------------------------------
data(scdataYML)
dim(avlscHEIDyml.m)

## ----chunk6, eval=T, echo=T---------------------------------------------------
act <- sciraEstRegAct(avlscHEIDyml.m,regnet = selreg.o$netTOI, norm = "z", ncores = 4);

## ----chunk7, eval=T, echo=T---------------------------------------------------
statTFAz.m <- act[,1:2];
colnames(statTFAz.m) <- c("t","P");
for(r in 1:nrow(act)){
    lm.o <- lm(act[r,] ~ phenoYML.v);
    statTFAz.m[r,] <- summary(lm.o)$coeff[2,3:4];
}
sig.idx <- which(statTFAz.m[,2] < 0.05/nrow(statTFAz.m));
print(summary(factor(sign(statTFAz.m[sig.idx,1]))));

## ----chunk8, message = FALSE,warning=F,fig.dim = c(8, 8),fig.align="center",fig.cap = "Heatmaps of TFA and gene expression in scRNA-Seq liver set"----
library(marray);
library(plotrix);
library(annotationTools);
library(Hmisc);
nTF <- nrow(act);

par(mfrow=c(1,2));

### cells already ordered according to developmental timepoint
actcolor.v <- maPalette(low="cyan",high="magenta",mid="grey",k=11);
breaks.v <- c(-10^6,-2,-1.5,-1,-0.5,-0.25,0.25,0.5,1,1.5,2,5,10^6);
dist.o <- as.dist(1-cor(t(act)));
hcl.o <- hclust(dist.o,method="complete");

par(mar=c(2,7,4,1));
image(y=1:nrow(act),x=1:ncol(act),z=t(act[hcl.o$order,]),col=actcolor.v,breaks=breaks.v,axes=FALSE,ylab="",xlab="");
#mtext("A)",side=3,at=-200,line=1,cex=1.25,font=2);
abline(v=0.5,col="black");
abline(v=length(phenoYML.v)+0.5,col="black");
for(h in seq(0.5,nTF+.5,1)){
abline(h=h,col="black",lwd=0.5);
}
for(i in 1:length(timept.v)){
   abline(v=length(which(phenoYML.v <= i))+0.5,col="black");
}
library(org.Hs.eg.db);
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
 
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
 
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}
geneS.v <- convertIDs(rownames(act[hcl.o$order,]), "ENTREZID","SYMBOL", org.Hs.eg.db,ifMultiple="useFirst");
axis(2,at=1:nrow(act),labels=geneS.v,font=3,las=2,cex=0.5);
mtext(side=2,"Liver-specific TFs",line=5,font=4);
mtext(side=1,"447 single cells",line=1,font=4,cex=0.75);
par(xpd=TRUE);
text(x=c(54/2,54+70/2,54+70+41/2,54+70+41+65/2,165+65+70/2,230+70+77/2,300+77+70/2),y=nTF+1,srt=45,labels=timept.v,cex=0.9)
par(xpd=FALSE);

expcolor.v <- maPalette(low="grey",mid="orange",high="red",k=8);
expbreaks.v <- c(0,0.25,0.5,1,2,4,6,8,100);

par(mar=c(2,5,4,3));
tmp.m <-  avlscHEIDyml.m[match(rownames(act[hcl.o$order,]),rownames(avlscHEIDyml.m)),];
image(y=1:nrow(tmp.m),x=1:ncol(tmp.m),z=t(tmp.m),col=expcolor.v,breaks=expbreaks.v,axes=FALSE,ylab="",xlab="");
#mtext("B)",side=3,at=-200,line=1,cex=1.25,font=2);
abline(v=0.5,col="black");
abline(v=length(phenoYML.v)+0.5,col="black");
for(h in seq(0.5,nTF+0.5,1)){
abline(h=h,col="black",lwd=0.5);
}
for(i in 1:length(timept.v)){
   abline(v=length(which(phenoYML.v <= i))+0.5,col="black");
}
axis(2,at=1:nrow(tmp.m),labels=geneS.v,font=3,las=2,cex=0.5);
mtext(side=1,"447 single cells",line=1,font=4,cex=0.75);
par(xpd=TRUE);
text(x=c(54/2,54+70/2,54+70+41/2,54+70+41+65/2,165+65+70/2,230+70+77/2,300+77+70/2),y=nTF+1,srt=45,labels=timept.v,cex=0.9)
par(xpd=FALSE);


par(xpd=TRUE);
text(x=525,y=18.5,labels=c("TFA"),cex=0.75,font=2);
text(x=525,y=10.5,labels="log(E)",cex=0.75,font=2);
text(x=c(54/2,54+70/2,54+70+41/2,54+70+41+65/2,165+65+70/2,230+70+77/2,300+77+70/2),y=nTF+1,srt=45,labels=timept.v,cex=0.9)
breaksPRX.v <- breaks.v;
breaksPRX.v[1] <- "< -10";
breaksPRX.v[length(breaks.v)] <- ">10";
color.legend(xl=550,yb=12,xr=580,yt=18,legend=breaksPRX.v,rect.col=actcolor.v,cex=0.5,align="rt",gradient="y")

breaksPRX.v <- expbreaks.v;
breaksPRX.v[length(expbreaks.v)] <- ">8";
color.legend(xl=550,yb=4,xr=580,yt=10,legend=breaksPRX.v,rect.col=expcolor.v,cex=0.5,align="rt",gradient="y");
par(xpd=FALSE);

## ----chunk9, message = FALSE,warning=F,fig.dim = c(8,4),fig.align="center",fig.cap = "Dynamics of HNF1A activity and comparison to DE"----

### HNF1A plot
par(mar=c(4,4,4,1));
act.v <- act[match(6927,rownames(act)),];
boxplot(act.v ~ phenoYML.v,names=timept.v,col=maPalette(low="grey",mid="skyblue",high="darkblue",k=7),pch=23,outline=FALSE,ylim=range(act.v),ylab="TFact",xlab="Stage",cex.lab=1.25,cex.axis=1);
mtext("HNF1A",side=3,line=0.5,cex=1,font=4);
ncpg.v <- summary(factor(phenoYML.v));
for(tp in 1:length(timept.v)){
    points(x=jitter(rep(tp,ncpg.v[tp]),amount=0.25),y=act.v[which(phenoYML.v==tp)],pch=23,cex=0.25,col="red");
}
text(x=2,y=4,paste("P=",signif(statTFAz.m[match(6927,rownames(statTFAz.m)),2],2),sep=""),font=2,cex=1.25);
mtext(side=1,at=1:7,line=2,paste("n=",ncpg.v,sep=""),cex=0.5);

## ----chunk10, message = FALSE,warning=F,fig.dim = c(8, 8),fig.align="center",fig.cap = "Scatterplot of PC1 vs PC2 from PCA applied to TF-activity matrix"----

act.pca <- prcomp(t(act)) 
color.v <-  rainbow(length(timept.v));
plot(act.pca$x[,1:2],col=color.v[phenoYML.v],pch=18)
legend("topleft",legend=timept.v,col = color.v,pch=18)


## ----chunk11, message = FALSE,warning=F,fig.dim = c(8, 4),fig.align="center",fig.cap = "The activity of known hepatocyte/cholangiocyte TFs"----

par(mfrow=c(1,2));
par(mar=c(4,4,2,1));
actcolor.v <- maPalette(low="cyan",high="magenta",mid="grey",k=11);
breaks.v <- c(-10^6,-2,-1.5,-1,-0.5,-0.25,0.25,0.5,1,1.5,2,5,10^6);

tmpEID.v <- convertIDs(c("HNF4A","IRF6"),"SYMBOL","ENTREZID",org.Hs.eg.db,ifMultiple="useFirst");
map.idx <- match(tmpEID.v,rownames(act));
i <- 1;
for(tf in map.idx){
    tfa.v <- act[tf,];
    colorTFA.v <- vector();
    for(br in 2:length(breaks.v)){
        colorTFA.v[which(tfa.v > breaks.v[br-1])] <- actcolor.v[br-1];
    }
    plot(act.pca$x[,1],act.pca$x[,2],col=colorTFA.v,pch=23,bg=colorTFA.v,xlab="PC1(12%)",ylab="PC2(7%)",cex=0.5);
    mtext(side=3,line=0.1,c("HNF4A","IRF6")[i],font=4);
    text(x=0,y=10,"Cholangiocyte-branch",font=4,cex=1);
    text(x=5,y=-5,"Hepatocyte-branch",font=4,cex=1);
    i <- i+1;
}


## ----chunk12, eval=T, echo=T--------------------------------------------------
# load in data
data(colonDATA);
# we now estimate regulatory activity
actTFz.m <- sciraEstRegAct(data = avlfpkmEID.m, regnet = netCOL.m, norm = "z", ncores = 4); ###

## ----chunk13, eval=T, echo=T--------------------------------------------------
statTFAz.m <- actTFz.m[,1:2];
colnames(statTFAz.m) <- c("t","P");
for(r in 1:nrow(actTFz.m)){
    lm.o <- lm(actTFz.m[r,] ~ statusNC.v);
    statTFAz.m[r,] <- summary(lm.o)$coeff[2,3:4];
}
    
nTF <- ncol(netCOL.m);
sig.idx <- which(statTFAz.m[,2] < 0.05/nTF);
print(summary(factor(sign(statTFAz.m[sig.idx,1]))));

## ----chunk14, eval=T, echo=T--------------------------------------------------
DoDEGwt <- function(exp.v,pheno.v){
       nspg.v <- summary(factor(pheno.v));
       wt.o <- wilcox.test(exp.v ~ pheno.v);
       auc <- 1-wt.o$stat/prod(nspg.v);
       out.v <- c(auc,wt.o$p.value);       
       return(out.v);
}
map.idx <- match(colnames(netCOL.m),rownames(avlfpkmEID.m));
statTF.m <- matrix(unlist(apply(avlfpkmEID.m[map.idx,],1,DoDEGwt,statusNC.v)),ncol=2,nrow=nrow(statTFAz.m),byrow=TRUE);
rownames(statTF.m) <- colnames(netCOL.m);
colnames(statTF.m) <- c("AUC","P");
    
nTF <- ncol(netCOL.m);
sig.idx <- which(statTF.m[,2] < 0.05/nTF);
print(summary(factor(sign(statTF.m[sig.idx,1]-0.5))));

## ----chunk15, message = FALSE,warning=F,fig.dim = c(8, 4),fig.align="center",fig.cap = "Differential activity vs differential expression patterns in colon cancer"----


layout(matrix(c(1,1,1,2,3,4,5,6),nrow=2,ncol=4,byrow=TRUE));

sigDNtfa.idx <- intersect(which(statTFAz.m[,2] < 0.05/nTF),which(statTFAz.m[,1] < 0));
sigUPtfa.idx <- intersect(which(statTFAz.m[,2] < 0.05/nTF),which(statTFAz.m[,1] > 0));
sigDN.idx <- intersect(which(statTF.m[,2] < 0.05/nTF),which(statTF.m[,1] < 0.5));
sigUP.idx <- intersect(which(statTF.m[,2] < 0.05/nTF),which(statTF.m[,1] > 0.5));

rest.idx <- setdiff(1:nrow(statTFAz.m),c(sigDNtfa.idx,sigUPtfa.idx));
ordTF.idx <- c(sigDNtfa.idx,sigUPtfa.idx,rest.idx);

geneS.v <- convertIDs(rownames(actTFz.m),from="ENTREZID",to="SYMBOL",org.Hs.eg.db,ifMultiple="useFirst");

bin.m <- matrix(0,ncol=2,nrow=nTF);
bin.m[sigDNtfa.idx,1] <- -1;
bin.m[sigUPtfa.idx,1] <- 1;
bin.m[sigDN.idx,2] <- -1;
bin.m[sigUP.idx,2] <- 1;

par(mar=c(3,4,5,1));
image(x=1:nrow(bin.m),y=1:ncol(bin.m),z=bin.m[ordTF.idx,2:1],col=c("blue","grey","brown"),breaks=c(-1.5,-0.5,0.5,1.5),axes=FALSE,xlab="",ylab="");
axis(2,at=1:2,labels=c("SCIRA","DE")[2:1],las=2);
axis(3,at=1:nrow(bin.m),labels=geneS.v[ordTF.idx],las=2,cex.axis=0.75);

bar.m <- matrix(c(20,3,8,5),nrow=2,ncol=2,byrow=TRUE);
barplot(t(bar.m),beside=TRUE,horiz=FALSE,ylab="#TF",col=c("blue","brown"),legend=TRUE,names=c("SCIRA","DE"));
pv <- pbinom(20,23,0.5,lower.tail=FALSE);
par(xpd=TRUE);
text(x=2,y=22,paste("P=",signif(pv,2),sep=""),font=2);
par(xpd=FALSE);
pv <- pbinom(8,13,0.5,lower.tail=FALSE);
text(x=5,y=10,paste("P=",signif(pv,2),sep=""),font=2);


tmpTFA.m <- actTFz.m[match(c("KLF5","ATOH1"),geneS.v),];
tmpEXP.m <- avlfpkmEID.m[match(rownames(actTFz.m),rownames(avlfpkmEID.m)),];
tmpTFE.m <- tmpEXP.m[match(c("KLF5","ATOH1"),geneS.v),];

tmpA.m <- statTFAz.m[match(c("KLF5","ATOH1"),geneS.v),]
tmpE.m <- statTF.m[match(c("KLF5","ATOH1"),geneS.v),];

par(mar=c(4,4,3,1));
for(i in 1:2){
nspg.v <- summary(factor(statusNC.v));
ylim.v <- range(tmpTFA.m[i,]);    
boxplot(tmpTFA.m[i,] ~ statusNC.v,names=c("N","C"),col=c("darkgreen","darkred"),pch=23,cex=0.25,outline=FALSE,ylab="TFA",xlab="",ylim=ylim.v);
text(x=1.5,y=c(3.5,4)[i],paste("P=",signif(tmpA.m[i,2],2),sep=""),font=2);
mtext(at=1:2,side=1,paste("n=",nspg.v,sep=""),line=1.75,cex=0.5);
mtext(side=3,line=0.1,c("KLF5","ATOH1")[i],font=4);
points(x=jitter(rep(1,nspg.v[1]),4),y=tmpTFA.m[i,statusNC.v==0],pch=23,cex=0.25)
points(x=jitter(rep(2,nspg.v[2]),3),y=tmpTFA.m[i,statusNC.v==1],pch=23,cex=0.25)

ylim.v <- range(tmpTFE.m[i,]);    
boxplot(tmpTFE.m[i,] ~ statusNC.v,names=c("N","C"),col=c("darkgreen","darkred"),pch=23,cex=0.25,outline=FALSE,ylab="log(FPKM+1)",xlab="",ylim=ylim.v);
text(x=1.5,y=c(12,2)[i],paste("P=",signif(tmpE.m[i,2],2),sep=""),font=2);
mtext(side=3,line=0.1,c("KLF5","ATOH1")[i],font=4);
mtext(at=1:2,side=1,paste("n=",nspg.v,sep=""),line=1.75,cex=0.5);
points(x=jitter(rep(1,nspg.v[1]),4),y=tmpTFE.m[i,statusNC.v==0],pch=23,cex=0.25)
points(x=jitter(rep(2,nspg.v[2]),3),y=tmpTFE.m[i,statusNC.v==1],pch=23,cex=0.25)

}


## ----sessionInfo, eval=T, echo=T----------------------------------------------
sessionInfo()

