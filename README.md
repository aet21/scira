> **Abstract: scira is an R-package aimed at estimating regulatory
> activity of transcription factors in scRNA-Seq data. It leverages the
> high-power of a large scale multi-tissue bulk gene expression dataset
> (GTEX) to build tissue-specific regulons, which can subsequently be
> applied to scRNA-Seq data to infer regulatory activity of
> tissue-specific transcription factors at single-cell resolution. It is
> particular aimed at scRNA-Seq studies profiling cancer or
> preneoplastic cells, or normal cells exposed to some disease risk
> factor, and may help to identify the tissue-specific regulatory
> networks that are disrupted early on in diseases like cancer.**

Motivation and Background
=========================

<p>
An early event in carcinogenesis are blocks in differentiation, driven
by disruption or inactivation of tissue-specific regulatory networks. It
is therefore of interest to map out the patterns of regulatory activity
for relevant transcripton factors (TFs) in cancer, preneoplastic lesions
and normal cells exposed to disease risk factors (Chen Y, Widschwendter
M, and Teschendorff AE 2017). Ideally, the inference of regulatory
activity should be performed at single-cell resolution, because
otherwise regulatory activity may be confounded with cell-type
heterogeneity. However, a key challenge of inferring regulatory activity
at single-cell resolution is the relatively high dropout rate, which can
in particular affect transcription factors themselves. We have developed
'scira' to help obtain regulatory activity estimates in scRNA-Seq data.
</p>
<p>
**scira** estimates regulatory activity in single cells by using
TF-regulons (Teschendorff AE and Wang 2020). A regulon is a set of
genes, enriched for direct binding targets of the TF, which can be used
to estimate the TF-binding activity of the TF. The members of the
regulon are faithful markers of upstream regulatory activity and the
mode of regulation can be positive or inhibitory. scira builds the
regulons using a high quality large-scale multi-tissue bulk gene
expression dataset, whilst also adjusting for stromal heterogeneity,
notably immune-cell contamination, which can otherwise confound
analysis. Specifically, we use the large multi-tissue GTEX dataset
encompassing 8555 samples from 30 different tissue types, to build the
tissue-specific regulons (Chen Y, Widschwendter M, and Teschendorff AE
2017). Because of the very large size, there is sufficient power for
most tissue-types to detect TFs and regulons that may however only be
expressed in a relatively low fraction of the resident cells within a
tissue.
</p>
<p>
In effect, **scira** consists of 2 steps:
<ol type="1">
<li>
Construction of a tissue-specific regulatory network consisting of
tissue-specific TFs and their regulons
</li>
<li>
Estimation of regulatory activity of the tissue-specific TFs in a
relevant scRNA-Seq dataset.
</li>
</ol>
</p>
<p>
The purpose of this vignette is to show how to implement **scira**, and
to demonstrate some of its potential applications. We shall first focus
on <span style="color:red">liver tissue</span>, inferring a
liver-specific regulatory network. We then validate the regulons using
independent bulk RNA-Seq data. We then showcase an application to a
timecourse differentiation scRNA-Seq experiment where hepatoblasts are
differentiated into hepatocytes and cholangiocytes, which serves to
demonstrate how **scira** can reveal this bifurcation, despite
cholangiocytes only making up 5-10% of liver cells in bulk liver tissue.
Finally, we consider another application to identify tumor suppressor
events at single-cell resolution in colon cancer.
</p>

Construction and validation of a tissue-specific regulatory network
===================================================================

In principle, the starting point of **scira** is the normalized
multi-tissue bulk RNA-Seq dataset from GTEX and a list of regulators
(transcription factors). These data are necessary to build the tissue
specific regulatory network. Because the GTEX dataset is a very large
matrix, encompassing 8555 samples and over 20,000 genes, and the
inference of the regulatory network takes a few minutes, for the
purposes of this vignette, we comment out the operations where the
network is built, and instead read in the inferred network.

Load in necessary libraries and build a liver-specific regulatory network
-------------------------------------------------------------------------

    library(parallel);
    library(corpcor);
    library(limma);
    library(scira);

    ## 

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

    data(tfEID); ### loads in the list of transcriptions factors (annotated to Entrez gene IDs)
    data(tissue); ### loads in tissue-type information of GTEX dataset (needed later)
    #regnet.o <- sciraInfReg(gtex.m,tfEID.v,sdth=0.25,sigth=1e-6,spTH=0.01,pcorth=0.2,minNtgts=10,ncores=4); ### this infers the non-tissue specific regulatory network but we comment this out as the GTEX data matrix, labeled here as gtex.m, is large and the inference takes a few minutes. 

As one can see, we have commented out the operation of the function
*sciraInfReg* which generates an initial non tissue specific network.
Next, one builds the tissue-specific network. In order to build it, we
need to specify a tissue of interest, which in this case is going to be
liver. The function to build the tissue-specific regulons is
*sciraSelReg* which uses the output from *sciraInfReg* as input. In
principle, we would thus run the following code, but because this step
still requires the large GTEX dataset, we instead read in the output
which is stored in the data object *selreg*:

    #selreg.o <- sciraSelReg(regnet.o,tissue.v,toi="Liver",cft.v=c("Blood"),degth.v=rep(0.05,2),lfcth.v=c(log2(1.5),log2(1))); ### generates tissue-specific network
    data(selreg);
    print(names(selreg.o));

    ## [1] "netTOI" "sumnet" "top"

The liver-specific regulons are specified in the *netTOI* element,
whereas a summary of the regulons is available in *sumnet*.

    dim(selreg.o$netTOI);

    ## [1] 18165    22

    head(selreg.o$sumnet);

    ##       nTGTS Act Rep
    ## 7494     18  18   0
    ## 51599    72  72   0
    ## 3172     31  28   3
    ## 633      93  93   0
    ## 3169     10  10   0
    ## 3175     15  15   0

So, we can see from this that the liver-specific regulatory network that
we have derived at the particular choice of parameters consists of 22
liver-specific TFs. Although there are 18165 targets, most of these have
zero entries and are not part of any liver-specific regulon. We can also
see that e.g. the regulon for transcription factor XBP1 (with Entrez
gene ID 7494) consists of 18 target genes, all of which are under
positive regulation upon activation of XBP1. In order to visualize the
liver specific regulatory network we can use the *sciraPlotNet*
function:

    netplot <- sciraPlotNet(selreg.o$netTOI);

    ## 'select()' returned 1:1 mapping between keys and columns
    ## 'select()' returned 1:1 mapping between keys and columns

    netplot;

    ## Warning: Removed 782 rows containing missing values (geom_text).

![The liver-specific regulatory
network](https://github.com/aet21/scira/tree/master/scira/vignettes/scira_files/figure-markdown_strict/chunk2c-1.png)

Estimation of regulatory activity and validation of liver-specific regulons
---------------------------------------------------------------------------

It is important to first validate the regulons in an independent bulk
tissue expression dataset, before estimating regulatory activity in
single cells. To this purpose, we use an independent RNA-seq dataset
from the ProteinAtlas project. Below, we load this data in, subsequently
estimate regulatory activity in each tissue sample from the
ProteinAtlas, and display the output:

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

<img src="scira_files/figure-markdown_strict/chunk3-1.png" alt="Validation of liver-specific regulons in ProteinAtlas dataset"  />
<p class="caption">
Validation of liver-specific regulons in ProteinAtlas dataset
</p>

As we can see, despite the small number of samples in each tissue, liver
is ranked very highly and is only superseeded by developmentally very
similar tissues like duodenum and small intestine. We can now compare
liver to all other tissue-types:

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

<img src="scira_files/figure-markdown_strict/chunk4-1.png" alt="Boxplot of the average activity level between liver and all other tissues"  />
<p class="caption">
Boxplot of the average activity level between liver and all other
tissues
</p>

We can see that the regulatory activity scores of the liver-specific TFs
are indeed significantly higher in liver-tissue compared to all other
tissues, as required.

Estimating regulatory activity in scRNA-Seq data: a liver differentiation timecourse
====================================================================================

Having validated the liver specific regulons, let us now estimate
regulatory activity of the liver-specific TFs in a scRNA-Seq dataset. We
choose a Fluidigm C1 timecourse differentiation experiment from Yang et
al (Yang L et al. 2017), that differentiated hepatoblasts at embyronic
day-10 (E10) into hepatocytes and cholangiocytes. There are seven
timepoints in total: E10, E11, E12, E13, E14, E15, E17. Although the
experiment was performed in mouse, for convenience we load in the data
with genes mapped to human homologs (Entrez gene IDs). The data has also
already been appropriately normalized.

    data(scdataYML)
    dim(avlscHEIDyml.m)

    ## [1] 16119   447

We can see that the scRNA-Seq data matrix is defined over 16119 genes
and 447 single-cells. The timepoint information is provided in the
*timept.v* and *phenoYML.v* vector objects. Let us now estimate
regulatory activity in these single cells:

    act <- sciraEstRegAct(avlscHEIDyml.m,regnet = selreg.o$netTOI, norm = "z", ncores = 4);

We now determine how many of the 22 liver specific TFs exhibit increased
regulatory activity with differentiation timepoint:

    statTFAz.m <- act[,1:2];
    colnames(statTFAz.m) <- c("t","P");
    for(r in 1:nrow(act)){
        lm.o <- lm(act[r,] ~ phenoYML.v);
        statTFAz.m[r,] <- summary(lm.o)$coeff[2,3:4];
    }
    sig.idx <- which(statTFAz.m[,2] < 0.05/nrow(statTFAz.m));
    print(summary(factor(sign(statTFAz.m[sig.idx,1]))));

    ## -1  1 
    ##  2 16

<p>
We note that estimated activity levels, unlike expression levels, are
not affected by dropouts and that therefore a linear model is sensible.
We can see that 16 of the 22 TFs exhibit significant increased
regulatory activity in the timecourse differentiation from hepatablasts
into hepatocytes and cholangiocytes. This makes a lot of biological
sense, because hepatocytes and cholangiocytes are two of the main liver
epithelial subtypes present in adult liver tissues, such as those from
GTEX where the original regulons were derived from. The fact that such a
high proportion of TFs are predicted to increase in activity further
validates the SCIRA algorithm.
</p>
Let us now display a heatmap of regulatory activity versus timepoint
alongside a corresponding heatmap of gene expression (log2(TMP+1.1)).

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

<img src="scira_files/figure-markdown_strict/chunk8-1.png" alt="Heatmaps of TFA and gene expression in scRNA-Seq liver set"  />
<p class="caption">
Heatmaps of TFA and gene expression in scRNA-Seq liver set
</p>

    par(xpd=FALSE);

As we can see, whilst the majority of liver-specific TFs do exhibit
increased regulatory activity with differentiation timepoint, this is
generally not the case when assessed using TF-expression levels (see
(Teschendorff AE and Wang 2020)). The next plot displays the TF-activity
score for one of the key TFs, HNF1A.

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

<img src="scira_files/figure-markdown_strict/chunk9-1.png" alt="Dynamics of HNF1A activity and comparison to DE"  />
<p class="caption">
Dynamics of HNF1A activity and comparison to DE
</p>

Now let us see if SCIRA can recapitulate the known bifurcation into
hepatocytes and cholangiocytes. We perform a PCA on the SCIRA-derived
transcription factor activity matrix defined over the 22 liver-specific
TFs and 447 single cells, which reveals a branching pattern emerging at
E14/E15.

    act.pca <- prcomp(t(act)) 
    color.v <-  rainbow(length(timept.v));
    plot(act.pca$x[,1:2],col=color.v[phenoYML.v],pch=18)
    legend("topleft",legend=timept.v,col = color.v,pch=18)

<img src="scira_files/figure-markdown_strict/chunk10-1.png" alt="Scatterplot of PC1 vs PC2 from PCA applied to TF-activity matrix"  />
<p class="caption">
Scatterplot of PC1 vs PC2 from PCA applied to TF-activity matrix
</p>

In order to check that the two branches emerging during E14-E15 indeed
relate to cholangiocytes and hepatocytes, we can superimpose the
activity of known hepatocyte/cholangiocyte factors. Here we use HNF4A
for hepatocytes, and IRF6 for cholangiocytes:

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

<img src="scira_files/figure-markdown_strict/chunk11-1.png" alt="The activity of known hepatocyte/cholangiocyte TFs"  />
<p class="caption">
The activity of known hepatocyte/cholangiocyte TFs
</p>

Thus, this confirms that PCA on the TF-activity matrix can recapitulate
the known bifurcation into hepatocytes and cholangiocytes. This is
extremely important observation, because we note that the original
regulons were derived from bulk liver tissue samples from GTEX, where
cholangiocytes only make up 5-10% of the cells in liver tissue. It is
because GTEX is highly powered (over 100 liver samples and over 8000
non-liver samples) that we can detect cholangiocyte-specific factors and
their regulons.

Application to colon cancer
===========================

Let us now consider an application to colon cancer. As before, we would
run the *sciraInfReg* and *sciraSelReg* functions to generate a
corresponding colon specific regulatory network, but instead load in the
colon-specific regulons. We then aim to use these regulons to infer
regulatory activity for the colon-specific TFs in a scRNA-Seq dataset
encompassing both normal colon and colon cancer cells (Li H et al.
2017). We have collated all the relevant data in the data object
*colonDATA*, which contains the normalized log-transformed scRNA-Seq
data set *avlfpkmEID.m*, the colon-specific regulons *netCOL.m* (a total
of 56 colon-TFs), as well as the normal-cancer status *statusNC.v* of
432 cells.

    # load in data
    data(colonDATA);
    # we now estimate regulatory activity
    actTFz.m <- sciraEstRegAct(data = avlfpkmEID.m, regnet = netCOL.m, norm = "z", ncores = 4); ###

As before, let us now determine how many of the 56 colon-specific TFs
exhibit differential activity between normal and cancer, and in which
direction they change:

    statTFAz.m <- actTFz.m[,1:2];
    colnames(statTFAz.m) <- c("t","P");
    for(r in 1:nrow(actTFz.m)){
        lm.o <- lm(actTFz.m[r,] ~ statusNC.v);
        statTFAz.m[r,] <- summary(lm.o)$coeff[2,3:4];
    }
        
    nTF <- ncol(netCOL.m);
    sig.idx <- which(statTFAz.m[,2] < 0.05/nTF);
    print(summary(factor(sign(statTFAz.m[sig.idx,1]))));

    ## -1  1 
    ## 20  3

Thus, of the 56 colon-specific TFs, 23 exhibit differential activity,
with the overwhelming majority of these (i.e 20 out of 23) exhibiting
inactivity in cancer compared to normal cells. Let us now contrast this
with differential expression:

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

    ## -1  1 
    ##  8  5

In the above, we have a Wilcoxon rank sum test, which is non-parametric,
to call differential expression (DE). Note that the statistic of a
Wilcoxon test can be related to the AUC, which is why the value 0.5
corresponds to the null-case of no difference. Using Wilcoxon rank sum
tests for DE in scRNA-Seq data is a well-accepted and popular procedure.
As we can see, according to DE only 13 are differentially expressed,
with only a moderate skew towards inactivation/downregulation.

It turns out that many of the colon-specific TFs are known or suspected
tumor suppressors in colon cancer, and thus the results obtained with
SCIRA are more credible and also consistent with DE changes in the bulk
tissue TCGA dataset (Teschendorff AE and Wang 2020). Let us display the
patterns of differential activity and differential expression, focusing
on two TFs (KLF5, ATOH1) which are known tumor suppressors in colon
cancer, but which DE fails to predict downregulation:

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

<img src="scira_files/figure-markdown_strict/chunk15-1.png" alt="Differential activity vs differential expression patterns in colon cancer"  />
<p class="caption">
Differential activity vs differential expression patterns in colon
cancer
</p>

As we can see, the two TFs KLF5 and ATOH1 do not exhibit convinding
downregulation patterns, yet SCIRA correctly predicts their inactivation
in cancer cells.

Session Info
============

    sessionInfo()

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] org.Hs.eg.db_3.8.2     AnnotationDbi_1.46.1   IRanges_2.18.3        
    ##  [4] S4Vectors_0.22.1       Biobase_2.44.0         BiocGenerics_0.30.0   
    ##  [7] Hmisc_4.4-0            ggplot2_3.3.1          Formula_1.2-3         
    ## [10] survival_3.1-12        lattice_0.20-41        annotationTools_1.58.0
    ## [13] plotrix_3.7-8          marray_1.62.0          scira_0.9.1           
    ## [16] limma_3.40.6           corpcor_1.6.9          BiocStyle_2.12.0      
    ## [19] markdown_1.1           knitr_1.28             rmarkdown_2.2         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit64_0.9-7          splines_3.6.3        network_1.16.0      
    ##  [4] BiocManager_1.30.10  highr_0.8            latticeExtra_0.6-29 
    ##  [7] blob_1.2.1           yaml_2.2.1           backports_1.1.7     
    ## [10] pillar_1.4.4         RSQLite_2.2.0        glue_1.4.1          
    ## [13] digest_0.6.25        checkmate_2.0.0      RColorBrewer_1.1-2  
    ## [16] colorspace_1.4-1     htmltools_0.4.0      Matrix_1.2-18       
    ## [19] plyr_1.8.6           pkgconfig_2.0.3      purrr_0.3.4         
    ## [22] scales_1.1.1         sna_2.5              jpeg_0.1-8.1        
    ## [25] htmlTable_1.13.3     tibble_3.0.1         generics_0.0.2      
    ## [28] farver_2.0.3         ellipsis_0.3.1       withr_2.2.0         
    ## [31] nnet_7.3-14          magrittr_1.5         crayon_1.3.4        
    ## [34] statnet.common_4.3.0 memoise_1.1.0        evaluate_0.14       
    ## [37] GGally_2.0.0         foreign_0.8-72       data.table_1.12.8   
    ## [40] tools_3.6.3          lifecycle_0.2.0      stringr_1.4.0       
    ## [43] munsell_0.5.0        cluster_2.1.0        compiler_3.6.3      
    ## [46] rlang_0.4.6          grid_3.6.3           rstudioapi_0.11     
    ## [49] htmlwidgets_1.5.1    base64enc_0.1-3      gtable_0.3.0        
    ## [52] DBI_1.1.0            reshape_0.8.8        R6_2.4.1            
    ## [55] gridExtra_2.3        dplyr_1.0.0          bit_1.1-15.2        
    ## [58] stringi_1.4.6        Rcpp_1.0.4.6         rpart_4.1-15        
    ## [61] vctrs_0.3.0          acepack_1.4.1        png_0.1-7           
    ## [64] tidyselect_1.1.0     xfun_0.14            coda_0.19-3

References
==========

Chen Y, Widschwendter M, and Teschendorff AE. 2017. “Systems-Epigenomics
Inference of Transcription Factor Activity Implicates
Aryl-Hydrocarbon-Receptor Inactivation as a Key Event in Lung Cancer
Development.” *Genome Biol* 18: 236.

Li H, Elise T Courtois, Debarka Sengupta, Yuliana Tan, Kok Hao Chen,
Jolene Jie Lin Goh, Say Li Kong, et al. 2017. “Reference Component
Analysis of Single-Cell Transcriptomes Elucidates Cellular Heterogeneity
in Human Colorectal Tumors.” *Nat Genet* 49 (5): 708–18.

Teschendorff AE, and Ning Wang. 2020. “Improved Detection of Tumor
Suppressor Events in Single-Cell Rna-Seq Data.” *Submitted*.

Yang L, Wei-Hua Wang, Wei-Lin Qiu, Zhen Guo, Erfei Bi, and Cheng-Ran Xu.
2017. “A Single-Cell Transcriptomic Analysis Reveals Precise Pathways
and Regulatory Mechanisms Underlying Hepatoblast Differentiation.”
*Hepatology* 66: 1387–1401.
