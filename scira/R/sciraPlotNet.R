#' @title
#' Visualization of a regulatory network
#'
#' @aliases sciraPlotNet
#'
#' @description
#' This function plots the regulatory network out
#'
#' @param net
#' the regulatory network in matrix format with columns labeling TFs and rows labeling gene targets and entries either 1, -1 or 0, depending on whether interactions are positive, negative or zero.
#'
#' @return A ggplot2 object
#' 
#' 
#' @examples
#' data(colonDATA)
#' sciraPlotNet(netCOL.m)
#'
#' @import ggplot2
#' @import GGally
#' @import network
#' @import org.Hs.eg.db
#' @import sna
#' @import scales
#' 
#' @export


sciraPlotNet <- function(net){
  

TF <- convertIDs(colnames(net), "ENTREZID","SYMBOL", org.Hs.eg.db,ifMultiple="useFirst")
TG <- convertIDs(rownames(net), "ENTREZID","SYMBOL", org.Hs.eg.db,ifMultiple="useFirst")

net <- net[-which(is.na(TG)==TRUE),]
TG <- TG[-which(is.na(TG)==TRUE)]


colnames(net) <- TF
rownames(net) <- TG

#########remove all zero TG
net <- net[-which(apply(net, 1, function(x) length(which(x==0)))== ncol(net)),]


#####square matrix

regulon.l<-apply(net,2,regulon.fun)

all.v <- c(names(regulon.l),unique(unname(unlist(regulon.l))))
net.m <- matrix(0,nrow =length(all.v),ncol = length(all.v) )
colnames(net.m) <- all.v
rownames(net.m) <- all.v

for (n in 1:length(regulon.l)) {
  net.m[n,which(rownames(net.m) %in% regulon.l[[n]])] <- 1
}

######network object
net.o <- network(net.m, matrix.type = "adjacency")
## set vertice type
type.v <- ifelse(get.vertex.attribute(net.o, "vertex.names") %in% colnames(net), "TF", "TG")
set.vertex.attribute(net.o,"mode",type.v)


Net.p <- ggnet2(net.o,node.size = 1.2,color= "mode",palette=c("TF" = "black", "TG" ="firebrick3"),edge.size = 0.3, edge.color = "grey",label= F )
Net.p$data$name <- Net.p$data$label
Net.p$data$name[!(Net.p$data$name %in% colnames(net))] <- NA

NetAdd.p <- Net.p + geom_text(aes(x = x, y = y, label = name),hjust = 0, vjust = 0, 
                              family="Times",  size = 3,fontface="bold")

return(NetAdd.p)

}

### Auxiliary Functions
#' @import org.Hs.eg.db

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

regulon.fun<-function(x){
  tf<-colnames(x)
  tg<-names(which(x==1))
}
