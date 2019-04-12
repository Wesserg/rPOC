#' lp1
#'
#' Function lp1 draws (by rejection) randomly any non-empty subset of the set 1:n with the same probability. 
#' @param A Description of the parameters
#' @keywords A
#' @export 
#' @examples
#' lp1()
lp1=function(A) {
	n=length(A)
	w=c()
	while (sum(w)==0) { # forbidding emptyset
		v=c()
		for (i in 1:n) {v=c(v,sample(0:1,1))} # generation of the 0-1 sequence of length n 
		w <- v
	}
	return(A[w==1]) # choose from A only elements at positions in which 1 was generated
}

#' lp2
#'
#' Function lp2 draws (by rejection) randomly any non-empty proper subset of the set 1:n with the same probability.
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' lp2(A)
lp2=function(A) {
	n=length(A)
	w=c()
	while (sum(w)==0 | sum(w)==n) { # forbidding emptyset and the improper subset
		v=c()
		for (i in 1:n) {v=c(v,sample(0:1,1))} # generation of the 0-1 sequence of length n 
		w <- v
	}
	return(A[w==1]) # choose from A only elements at positions in which 1 was generated
}


#' lp3
#'
#' Function lp3 draws (by rejection) randomly any subset of size >1 of the set 1:n with the same probability.
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' lp3(A)
lp3=function(A) {
	n=length(A)
	w=c()
	while (sum(w)==0 | sum(w)==1) {	 # forbidding emptyset and singletons
		v=c()
		for (i in 1:n) {v=c(v,sample(0:1,1))} # generation of the 0-1 sequence of length n 
		w <- v
	}
	return(A[w==1]) # choose from A only elements at positions in which 1 was generated
}


#' poc1
#'
#' FUNCTION poc1 DRAWS (RANDOMLY, Non-uniformly) A PERFECT ORDERING OF CLIQUES OF n NODES	 (AND THUS A DECOMPOSABLE GRAPH) ITS OUTPUT CLI IS THE INPUT TO THE MAIN PROCEDURE
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' poc1(10)
poc1=function(n) {
	B=1:n
	I<-n
	C=lp3(B) # draw the first clique
	CLI=list(C) # put it on the list of cliques
	B<-setdiff(B,C) # what remains
	I<-length(B) # size what remains
	if (I==0) {
		SEP=list() # in case the first clique is the whole set the set of separators is empty
	} else {
		S=lp2(C) # otherwise draw the separator S2 from the first clique
		SEP=list(S) # and put it on the list of separators
	}	 
	while (I>0) { #	 proceed further only if what remains is non-empty
		R=lp1(B) # draw subsequent residual from what remained
		C=sort(union(S,R)) # create the subsequent clique 
		CLI=c(CLI,list(C)) # put in on the list of cliques
		B=setdiff(B,C) # what remains
		I=length(B) # size of what remains
		if (I==0) break # stop if empty set remains
		k=length(CLI) # number of cliques already on the list 
		L=sample(1:k,1) # pick up at random one of the cliques which are on the list at present
		S=lp2(CLI[[L]]) # draw the subsequent separator from the chosen clique
		SEP=c(SEP,list(S)) # put it on the list of separators
	}
	#names(CLI)=c(paste0("C",1:length(CLI))) # labelling cliques
	#for (i in 1:length(SEP)) names(SEP[[i]])=paste0("S",i+1) # labelling separators
	
	return(list(CLI=CLI,SEP=SEP))
}

#' CHECK1
#'
#' function CHECK checks if CLIB is a perfect ordering of cliques
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' CHECK1(CLIQUES)
CHECK1=function(CLIQUES)
{
		n=length(CLIQUES)
    UCLI=lapply(1:n, function(k) unique(unlist(CLIQUES[1:k])))
    UCLIS=lapply(2:n, function(k) intersect(CLIQUES[[k]], UCLI[[k-1]]))

		for (j in 1:(n-1))
    {
				for (i in 1:j)
        {
            ij_condition=prod(UCLIS[[j]] %in% CLIQUES[[i]])==1
            if (ij_condition) break;
				}
        if (!(ij_condition)) return(FALSE);
		}
		return(TRUE)
}


#' CHECK1_OLD
#'
#' function CHECK checks if CLIB is a perfect ordering of cliques
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' CHECK1_OLD(CLIQUES)
CHECK1_OLD=function(CLI)
{
		UCLI=list(CLI[[1]])
		kk=length(CLI)
		for (k in 2:kk) UCLI=c(UCLI,list(union(UCLI[[k-1]],CLI[[k]])))
		UCLI
		UCLIS=list(c(1))
		for (k in 2:kk) UCLIS=c(UCLIS,list(intersect(CLI[[k]],UCLI[[k-1]])))
		UCLIS

		bb=0
		for (k in 2:kk) { 
				b=1
				for (j in 1:(k-1)) {
						b=min(b,length(setdiff(UCLIS[[k]],CLI[[j]])))
						if (b==0) {
								break
						}
				}
				bb=b+bb
		}
		bb
		if (bb==0) {bbb=1}
		if (bb>0) {bbb=0}
		return(bbb)
}

#' Wmat
#'
#' FUNCTIONS Wmat AND matmod FORM THE TECHNICAL CORE OF THE ALGORITHM WHICH CONSTRUCTS A PERFECT PERMUTATION OF CLIQUES (ASSUMING THAT IT IS POSSIBLE). function Wmat constructs the weighted coincidence matrix for the set	of cliques CLI with weight w(i,j) = number of vertices in C_i cap C_j
#' @param A Description of the parameter
#' @keywords A
#' @export
#' @examples
#' Wmat(CLIQUES)
Wmat=function(CLI)
{
    n=length(CLI)
    W=array(0,c(n,n))
    for (i in 1:(n-1))
    {
        for (j in (i+1):n)
        {
            W[i,j]=length(intersect(CLI[[i]],CLI[[j]]))
        }
    }
    W=W+t(W)
    return(W)
}
#' matmod_new
#'
#' FUNCTIONS Wmat AND matmod FORM THE TECHNICAL CORE OF THE ALGORITHM WHICH CONSTRUCTS A PERFECT PERMUTATION OF CLIQUES (ASSUMING THAT IT IS POSSIBLE). function matmod puts zeros in the weighted	incidence matrix A	for those edges connecting vertices from R to vertices from R^c which are removable with respect to R
#' @param A Description of the parameter
#' @keywords W
#' @export 
#' @examples
#' matmod(W,R,CLI)
matmod=function(W,A,CLI=CLI)
{
    n=length(CLI)
#    AC=setdiff(1:n, A)
    W.not_in.A=array(TRUE, c(n,n)) #initialize as if all edges are not in A
    if (length(A)>1) #Find edges in A and mark with "FALSE"
    {
        all.A.pairs=t(combn(A,2))
        W.not_in.A[all.A.pairs]=FALSE 
        W.not_in.A=W.not_in.A & t(W.not_in.A) #Make it symmetric
    }
    sorted.unique.weights=sort(unique(as.vector(W)))[-1] #Only iterate over actually present weights
    for (k in sorted.unique.weights)
    {
        W.less_k=(W < k)
        W2=W.less_k & W.not_in.A #Removing candidates on level 'k'
        W1=(W2==FALSE) #Not removable
        W1_m=W1
        m=1
        while (m  < (n-1))
        {
            m=m+1
            W1_m=W1_m%*%W1
            W[W2*W1_m>0]=0 #remove edge if there exists a path of length m that goes through non-removable edges on level k
        }
    }
    return(W)
}

#' matmod_old
#'
#' FUNCTIONS Wmat AND matmod FORM THE TECHNICAL CORE OF THE ALGORITHM WHICH CONSTRUCTS A PERFECT PERMUTATION OF CLIQUES (ASSUMING THAT IT IS POSSIBLE). function matmod puts zeros in the weighted	incidence matrix A	for those edges connecting vertices from R to vertices from R^c which are removable with respect to R
#' @param A Description of the parameter
#' @keywords W
#' @export 
#' @examples
#' matmod(W,R,CLI)
matmod_old=function(W,A,CLI=CLI) {
	K1=max(W)
	K2=min(W[W>0])
	
	W3=W
	for (k in (K2+1):K1) 
	{
		W1=W3
		W2=W3
		for (i in 1:length(CLI)) 
		{
			for (j in 1:length(CLI)) {
				if (W3[i,j]<k & (length(setdiff(c(i,j),A))>0)) W1[i,j]=0 
				if ((W3[i,j]>=k) | (length(setdiff(c(i,j),A))==0)) W2[i,j]=0
			}
		}

		W12=W1
		B<-matrix(1,ncol=length(CLI),nrow=length(CLI))
		for (m in (1:(length(CLI)-1)))
		{
				W12=W12%*%W1
				for (i in A){
						for (j in 1:length(CLI)) {
								if (W2[i,j]*W12[i,j]>0) {
										B[i,j]=0
										B[j,i]=0
								}
						}
				}
		}

		W3=W3*B

	}
    return(W3) 
}


#' RECURRENT.POC.SEARCH
#'
#' THE CORE OF THE RECURSION WHICH PRINTS ALL PERFECT ORDERING OF "CLIQUES" IS FUNCTION: "RECURRENT.POC.SEARCH" THE FUNCTION "RECURRENT.POC.SEARCH.WRAPPER"	 JUST SETS THE STAGE	
#' @param W Description of the parameter
#' @keywords W
#' @export 
#' @examples
#' RECURRENT.POC.SEARCH(...)
RECURRENT.POC.SEARCH <- function(W, W.labels, o.k, CLIQUES)
{
	n=dim(W)[1]
	B=W
	o.k.c=setdiff(1:n, o.k)
	if (length(o.k)>0 & length(o.k.c)<n)
	{
		W=matmod(W,o.k, CLI=CLIQUES) # killing the unwanted choices of cliques with respect to o.k 
		B=W[o.k,o.k.c, drop=F] # creates the matrix with columns labels of cliques not-chosen yet
	}
	B.ColSum=colSums(B) # sums of columns in B
	B.labels=W.labels[o.k.c] # labels of cliques which have not been chosen yet
	if (n==2)
	{
      ## Using the fact that permutation of initial two elements of POC yields also a POC
      assign("all.POCs", list(c(2,1), c(1,2)),envir = .GlobalEnv)
	}
	for (j in B.labels[B.ColSum>0])
	{
			o.k.j=c(o.k,j) # vector of labels of already chosen cliques
			if (length(o.k.j)==n)
			{
          ## Using the fact that permutation of initial two elements of POC yields also a POC
          assign("all.POCs", c(all.POCs, o.k.j),envir = .GlobalEnv)
          assign("all.POCs", c(all.POCs, o.k.j[c(2,1,3:length(o.k.j))]),envir = .GlobalEnv)
			}
      ## Using the fact that permutation of initial two elements of POC yields also a POC
			if (length(o.k.j)<2) RECURRENT.POC.SEARCH(W,W.labels, o.k=o.k.j, CLIQUES); # recursion step if o.k.j has single clique
			if (length(o.k.j)>=2 & o.k.j[1] < o.k.j[2]) RECURRENT.POC.SEARCH(W,W.labels, o.k=o.k.j, CLIQUES); # recursion step. ensuring we will not double the outcomes with investigating only one order of initial two elements of a POC

	}
}

#' RECURRENT.POC.SEARCH.WRAPPER
#'
#' THE CORE OF THE RECURSION WHICH PRINTS ALL PERFECT ORDERING OF "CLIQUES" IS FUNCTION: "RECURRENT.POC.SEARCH" THE FUNCTION "RECURRENT.POC.SEARCH.WRAPPER"	 JUST SETS THE STAGE	
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' RECURRENT.POC.SEARCH.WRAPPER(...)
RECURRENT.POC.SEARCH.WRAPPER <- function(CLIQUES)
{
    assign("all.POCs", c(), envir = .GlobalEnv)
    W=Wmat(CLIQUES)
    W.labels=col(W)[1,]
    o.k=c()
    out=RECURRENT.POC.SEARCH(W, W.labels, o.k, CLIQUES)
    POCs=matrix(all.POCs, ncol=length(CLIQUES), byrow=TRUE)
    return(POCs)
}

#' get.all.POCs
#'
#' get.all.POCs descroption
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' get.all.POCs(...)
get.all.POCs <- function(CLIQUES, method)
{
    if (method=='recurrent')
    {
        POCs=RECURRENT.POC.SEARCH.WRAPPER(CLIQUES)
    }
    if (method=='check.all')
    {
        POCs=c()
        all.permutations=permn(1:length(CLIQUES))
        if (length(all.permutations)==2) POCs=list(c(1,2), c(2,1))
        else
            {
                for (i in 1:length(all.permutations))
                {
                    perm.i=all.permutations[[i]]
                    if (perm.i[1] < perm.i[2] & length(perm.i)>2)
                    {
                        if (CHECK1(CLIQUES[perm.i])) POCs=rbind(POCs,perm.i, perm.i[c(2,1,3:length(perm.i))]);
                    }
                }
            }
    }
    return(POCs)
}


#' single.POC
#'
#' single.POC descroption
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' single.POC(...)
single.POC <- function(W, W.labels, o.k, CLIQUES)
{
		n=dim(W)[1]
		B=W
		o.k.c=setdiff(1:n, o.k)
		if (length(o.k)>0 & length(o.k.c)<n)
		{
				W=matmod(W,o.k, CLI=CLIQUES) # killing the unwanted choices of cliques with respect to o.k 
				B=W[o.k,o.k.c, drop=F] # creates the matrix with columns labels of cliques not-chosen yet
		}
		B.ColSum=colSums(B) # sums of columns in B
		B.labels=W.labels[o.k.c] # labels of cliques which have not been chosen yet
    if (length(B.labels[B.ColSum>0])==1) j=B.labels[B.ColSum>0];
    if (length(B.labels[B.ColSum>0])>1)  j = sample(B.labels[B.ColSum>0],1);
    ##j = B.labels[B.ColSum>0][1]
    o.k.j=c(o.k,j) # vector of labels of already chosen cliques
    if (length(o.k.j)==n)
    {
        #print(o.k.j)
        POC<<-o.k.j
    }
    if (length(o.k.j)<n) single.POC(W,W.labels, o.k=o.k.j, CLIQUES); # recursion step
}

#' get.single.POC
#'
#' get.single.POC descroption
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' get.single.POC(...)
get.single.POC <- function(CLIQUES)
{
		W=Wmat(CLIQUES)
		W.labels=col(W)[1,]
		o.k=c()
		single.POC(W, W.labels, o.k, CLIQUES)
		return(POC)
}

#' is.decomposable
#'
#' get.single.POC descroption
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' is.decomposable(CLIQUES)
is.decomposable <- function(CLIQUES)
{
		POC.candidate=get.single.POC.wrapper(CLIQUES)
		return(CHECK1(CLIQUES[POC.candidate]))
}

#' rPOC
#'
#' rPOC description
#' @param n Description of the parameter
#' @keywords random sample
#' @export 
#' @examples
#' rPOC.unif(1,CLIQUES)
rPOC_recurrent <- function(n,CLIQUES,method='recurrent', replace=TRUE)
{
    nPOCs=list()
    if (method=='recurrent')
    {
        if (replace)
        {
            nPOCs=replicate(get.single.POC(CLIQUES),n)
        }
        if (!(replace))
        {
            if (n>factorial(length(CLIQUES)))
            {
                return(nPOCs)
            }
            counter=1
            while (counter<=n)
            {
                POC.candidate=get.single.POC(CLIQUES)
                if (!(list(POC.candidate) %in% nPOCs))
                {
                    nPOCs[[counter]]=POC.candidate
                    counter=counter+1
                }
            }

        }

    }
    return(nPOCs)
}

#' rPOC_unif
#'
#' rPOC description
#' @param n Description of the parameter
#' @keywords random sample
#' @export 
#' @examples
#' rPOC.unif(1,CLIQUES)
rPOC_unif<-function(n,CLIQUES,method='recurrent', replace=TRUE)
{
    ## it should return an array of POCs of size (n, len(Cliques))
    nPOCs=array(0,c(n, length(CLIQUES)))
    if (method=='rejection')
    {
        if (replace)
        {
            counter=1
            while (counter<=n)
            {
                POC.candidate=sample(1:length(CLIQUES))
                if (CHECK1(CLIQUES[POC.candidate]))
                {
                    nPOCs[counter,]=POC.candidate
                    counter=counter+1
                }
            }
        }
        if (!(replace))
        {
            if (n>factorial(length(CLIQUES)))
            {
                print('Warning: sampling space smaller than sample size')
                return(nPOCs)
            }
            all.permutations=permn(1:length(CLIQUES))
            counter=1
            while (counter<=n)
            {
                POC.candidate=unlist(sample(all.permutations,1))
                if (!(list(POC.candidate) %in% nPOCs) && CHECK1(CLIQUES[POC.candidate]))
                {
                    nPOCs[counter,]=POC.candidate
                    counter=counter+1
                }
            }

        }

    }
    if (method=='recurrent')
    {
        POCs=RECURRENT.POC.SEARCH.WRAPPER(CLIQUES)
        if (replace)
        {
            nPOCs=POCs[sample(1:dim(POCs)[1],n, replace=TRUE),,drop=FALSE]
        }
        if (!(replace))
        {
            nPOCs=POCs[sample(1:dim(POCs)[1],n, replace=FALSE),,drop=FALSE]
        }
    }

    return(nPOCs)
}

