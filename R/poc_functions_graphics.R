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
Get.graph<-function(CLIQUES, draw=FALSE)
{
    edges=unique(matrix(unlist(lapply(CLIQUES, function(CLIQUE) combn(CLIQUE,2))),ncol=2, byrow=T))
    nodes=sort(unique(unlist(CLIQUES)))
    G = graph.data.frame(edges, nodes, directed=F)

    if (draw)
    {
        png(paste('Figures/',length(CLIQUES),'_Clique_Graph.png', sep="", collapse=""), height=1024, width=1024)
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
Get.weighted.clique.graph<-function(CLIQUES, draw=FALSE)
{
    W=Wmat(CLIQUES)
    WCG=graph_from_adjacency_matrix(W$A, mode =  "undirected", weighted = TRUE)
    if (draw)
    {
        png(paste('Figures/',length(CLIQUES),'_Weighted_Clique_Graph.png', sep="", collapse=""), height=1024, width=1024)
        plot(WCG,
             vertex.color=add.alpha(brewer.pal(length(CLIQUES),"Set3"),0.75),
             vertex.size=unlist(lapply(CLIQUES,length)),
             vertex.frame.color='black',
             vertex.label.cex=2,
             vertex.label.color='black',
             vertex.label=paste('Clique #',V(WCG), sep=""),
             edge.color=add.alpha('black',0.5),
             edge.width=E(WCG)$weight*2,
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
Get.ALLPOCs.graph<-function(POCname)
{
    filename=POCname
    POClist=read.csv(filename, header=FALSE)
    L=t(apply(POClist,1,function(POC) sapply(length(POC):1, function(i) paste(c(POC[1:i]), collapse=""))))

    edge.list=lapply(1:dim(L)[1],function(i) t(matrix(c(L[i,-1],L[i,-dim(L)[2]]), ncol=2)))
    edges=count(matrix(unlist(edge.list),ncol=2, byrow=TRUE ))
    ## edges=matrix(unlist(edge.list),ncol=2, byrow=TRUE )
    edges=cbind(edges,sapply(edges[,1], function(x) substring(x,1,1)))
    colnames(edges)=c("V1","V2","freq","group.label")

    nodes=unique(c(unlist(edge.list)))
    nodes=cbind(nodes,sapply(nodes, function(x) substring(x,1,1)),sapply(nodes, function(x) nchar(x)))
    nodes=rbind(nodes, c('Root',"0","1"))
    colnames(nodes)=c('nodes','group.label','level')

    G = graph.data.frame(edges, nodes, directed=F)
    G=add_edges(G, unlist(lapply(1:dim(POClist)[2], function(i) c("Root",paste(i)))), freq=max(edges[,'freq']) )
    png(paste('Figures/',POCname,'_Summary.png',sep=""), height=256*sqrt(dim(POClist)[1]), width=256*sqrt(dim(POClist)[1]))
    colors=brewer.pal(dim(POClist)[2]+1,"Set2")
    plot(G,layout=layout_as_tree(G , circular=TRUE, root='Root'),
         vertex.color=colors[as.numeric(nodes[,'group.label'])],
         edge.width=edges[,'freq'],
         vertex.size=sqrt(1/as.numeric(nodes[,'level']))*dim(POClist)[2],
         vertex.label.cex=sqrt(1/as.numeric(nodes[,'level']))*dim(POClist)[2],

         vertex.label.dist=0)
    dev.off()
}
