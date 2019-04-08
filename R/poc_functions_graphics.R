#' add.alpha
#'
#' Add an alpha value to a colour
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' add.alpha(col, alpha)
add.alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha))  
}

#' Get.graph
#'
#' Add an alpha value to a colour
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' Get.graph(CLIQUES, draw=TRUE)
Get.graph<-function(CLIQUES, draw=FALSE, prefix='')
{
    edges=unique(matrix(unlist(lapply(CLIQUES, function(CLIQUE) combn(CLIQUE,2))),ncol=2, byrow=T))
    nodes=sort(unique(unlist(CLIQUES)))
    G = graph.data.frame(edges, nodes, directed=F)

    if (draw)
    {
        png(paste('Figures/',prefix,'_G.png', sep="", collapse=""), height=1024, width=1024)
        plot(G,
             mark.groups=CLIQUES,
             mark.col=add.alpha(brewer.pal(length(CLIQUES),"Set3"),0.5),
             mark.border=NA,
             vertex.color='white',
             vertex.size=15,
             vertex.frame.color='black',
             edge.color='black',
             edge.width=2)
        dev.off()

    }
    return(G)
}

#v' Get.weighted.clique.graph
#'
#' Add an alpha value to a colour
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' Get.weighted.clique.graph(CLIQUES, draw=TRUE)
Get.weighted.clique.graph<-function(CLIQUES, draw=FALSE, prefix='')
{
    W=Wmat(CLIQUES)
    WCG=graph_from_adjacency_matrix(W, mode =  "undirected", weighted = TRUE)
    if (draw)
    {
        png(paste('Figures/',prefix,'_WCG.png', sep="", collapse=""), height=1024, width=1024)
        plot(WCG,
             vertex.color=add.alpha(brewer.pal(length(CLIQUES),"Set3"),0.75),
             vertex.size=prop_to_size(unlist(lapply(CLIQUES,length))*2,10,20,1),
             vertex.frame.color='black',
             vertex.label.cex=2,
             vertex.label.color='black',
             vertex.label=paste('C',V(WCG), sep=""),
             edge.color=add.alpha('black',0.5),
             edge.width=prop_to_size(E(WCG)$weight*2,10,20,1),
             edge.label=E(WCG)$weight,
             edge.label.cex=3,
             edge.label.color='black')
        dev.off()

    }
    
    return(WCG)
}

#' Get.ALLPOCs.graph
#'
#' Add an alpha value to a colour
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' Get.ALLPOCs.graph(CLIQUES, draw=TRUE)
Get.ALLPOCs.graph<-function(POClist, prefix='')
{
    dim.scaling=max(512,256*sqrt(dim(POClist)[1]))
    L=t(apply(POClist,1,function(POC) sapply(length(POC):1, function(i) paste(c(POC[1:i]), collapse=""))))

    edge.list=lapply(1:dim(L)[1],function(i) t(matrix(c(L[i,-1],L[i,-dim(L)[2]]), ncol=2)))
    edges=count(matrix(unlist(edge.list),ncol=2, byrow=TRUE ))
    edges=cbind(edges,sapply(edges[,1], function(x) substring(x,1,1)))
    colnames(edges)=c("V1","V2","freq","group.label")

    nodes=unique(c(unlist(edge.list)))
    nodes=cbind(nodes,sapply(nodes, function(x) substring(x,1,1)),sapply(nodes, function(x) nchar(x)))
    nodes=rbind(nodes, c(' ',"8","0.5"))
    colnames(nodes)=c('nodes','group.label','level')

    G = graph.data.frame(edges, nodes, directed=F)
    G=add_edges(G, unlist(lapply(1:dim(POClist)[2], function(i) c(" ",paste(i)))), freq=sum(edges[,'freq']) )
    png(paste('Figures/',prefix,'_Bouquet.png',sep=""), height=dim.scaling, width=dim.scaling)
    colors=brewer.pal(dim(POClist)[2]+1,"Set2")
    plot(G,layout=layout_as_tree(G , circular=TRUE, root=' '),
         vertex.color=colors[as.numeric(nodes[,'group.label'])],
         edge.width=prop_to_size(edges[,'freq'],10,30, 1),
         vertex.size=prop_to_size(sqrt(1/as.numeric(nodes[,'level'])),1,4,0.5)*dim(POClist)[2]*2/log(dim(POClist)[1]),
         vertex.label.cex=sqrt(1/as.numeric(nodes[,'level']))*dim(POClist)[2],

         vertex.label.dist=0)
    dev.off() 
}
#' prop_to_size
#'
#' Resize properties of a graph
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' prop_to_size()
prop_to_size<-function(prop, mi=0, ma=1, power=2)
{
    if (sd(prop)<0.001)
    {
        size=mi
    }
    else
    {
        size= mi + (ma - mi)*((prop-min(prop))/(max(prop)-min(prop)))^power;
    }
    return(size)
}
