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
Get.graph<-function(CLIQUES, draw=FALSE, prefix='', proba = '')
{
    CLIQUES_with_edges = CLIQUES[lengths(CLIQUES) > 1]
    edges=unique(matrix(unlist(lapply(CLIQUES_with_edges, function(CLIQUE) combn(sort(CLIQUE),2))),ncol=2, byrow=T))
    nodes=sort(unique(unlist(CLIQUES)))
    G = graph.data.frame(edges, nodes, directed=FALSE)
    nCLIQUES=length(CLIQUES)
    CLIQUE.strings=lapply(
        1:nCLIQUES,
        function(i)
            paste('C',i,' = ( ', paste(CLIQUES[[i]], collapse=' '), ' )', sep='',collapse=''))
 
    if (draw)
    {
        sub = paste(c("Transition Probability: ", proba), sep="", collapse = " ")
        main = paste(CLIQUE.strings, sep="", collapse = " ")
        print(main)
        print(sub)
        png(paste('Figures/',prefix,'_G.png', sep="", collapse=""), height=1024, width=1024)
        plot(G,
             mark.groups=CLIQUES,
             mark.col=add.alpha(brewer.pal(nCLIQUES,"Set3"),0.5),
             mark.border=NA,
             vertex.color='white',
             vertex.size=20,
             vertex.label.cex=4,
             vertex.frame.color='black',
             edge.color='black',
             edge.width=3,
             #main=main,
             #sub=sub,
             )
        title(sub,cex.main=3)
        legend('topleft', legend=CLIQUE.strings, cex=3, fill=brewer.pal(nCLIQUES,"Set3"))
        dev.off()

    }
    return(G)
}

#' Get.weighted.clique.graph
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
             vertex.color=add.alpha(brewer.pal(length(CLIQUES),"Set3"),1),
             vertex.size=prop_to_size(unlist(lapply(CLIQUES,length))*2,25,50,1),
             vertex.frame.color='black',
             vertex.label.cex=4,
             vertex.label.color='black',
             vertex.label=paste('C',V(WCG), sep=""),
             edge.color=add.alpha('black',0.5),
             edge.width=prop_to_size(E(WCG)$weight*2,60,80,1),
             edge.label=E(WCG)$weight,
             edge.label.cex=5,
             edge.label.color='black')
        dev.off()

    }
    
    return(WCG)
}

#' Get.reduced.weighted.clique.graph
#'
#' Add an alpha value to a colour
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' Get.reduced.weighted.clique.graph(CLIQUES, draw=TRUE)
Get.reduced.weighted.clique.graph<-function(CLIQUES, draw=FALSE, prefix='')
{
    W=Wmat(CLIQUES)
    WCG=graph_from_adjacency_matrix(W, mode =  "undirected", weighted = TRUE)
    RW=matmod(W,c())
    RWCG=graph_from_adjacency_matrix(RW, mode =  "undirected", weighted = TRUE)

    removable_edges=which(W-RW>0, arr.ind=TRUE)
    E(WCG)$color=add.alpha('black',0.5)
    WCG_edges=get.edgelist(WCG)
    removable_edges_id=c()
    if (dim(removable_edges)[1]>0)
    {
        removable_edges_id=which(apply(WCG_edges,1,function(WCG_edge) apply(removable_edges, 1, function(removable_edge) all(WCG_edge==removable_edge))), arr.ind=TRUE)[,2]
    }

    E(WCG)[removable_edges_id]$color=add.alpha('red',0.5)
    if (draw)
    {
        png(paste('Figures/',prefix,'_RWCG.png', sep="", collapse=""), height=1024, width=1024)
        plot(WCG,
             vertex.color=add.alpha(brewer.pal(length(CLIQUES),"Set3"),1),
             vertex.size=prop_to_size(unlist(lapply(CLIQUES,length))*2,25,50,1),
             vertex.frame.color='black',
             vertex.label.cex=4,
             vertex.label.color='black',
             vertex.label=paste('C',V(WCG), sep=""),
             edge.color=E(WCG)$color, #add.alpha('black',0.5),
             edge.width=prop_to_size(E(WCG)$weight*2,40,80,1),
             edge.label=E(WCG)$weight,
             edge.label.cex=5,
             edge.label.color='black')
        dev.off()

    }
    
    return(RWCG)
}

#' Get.ALLPOCs.graph
#'
#' Add an alpha value to a colour
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' Get.ALLPOCs.graph(CLIQUES, draw=TRUE)
Get.ALLPOCs.graph<-function(POClist, prefix='', background=FALSE)
{
    dim.scaling=512*dim(POClist)[2]
    L=t(apply(POClist,1,function(POC) sapply(length(POC):1, function(i) paste(c(POC[1:i]), collapse=""))))

    edge.list=lapply(1:dim(L)[1],function(i) t(matrix(c(L[i,-1],L[i,-dim(L)[2]]), ncol=2)))
    edges=count(matrix(unlist(edge.list),ncol=2, byrow=TRUE ))
    edges=cbind(edges,sapply(edges[,1], function(x) substring(x,1,1)))
    colnames(edges)=c("V1","V2","freq","group.label")

    nodes=unique(c(unlist(edge.list)))
    nodes=cbind(nodes,sapply(nodes, function(x) substring(x,1,1)),sapply(nodes, function(x) nchar(x)))
    nodes=rbind(nodes, c(' ',"8","0.5"))
    colnames(nodes)=c('nodes','group.label','level')

    if (background)
    {
        all.permutations=matrix(unlist(permn(1:dim(POClist)[2])),ncol=dim(POClist)[2], byrow=T)
        add_sequence=apply(all.permutations, 1, function(permutation)any(apply(POClist,1, function(POC) permutation==POC)))
    }

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
    legend('topleft', legend=paste('#POCs = ',dim(POClist)[1],sep=''), cex=dim(POClist)[2]+3)
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
