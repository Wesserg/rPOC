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
CHECK1=function(CLIB)
{
		UCLI=list(CLIB[[1]])
		kk=length(CLIB)
		for (k in 2:kk) UCLI=c(UCLI,list(union(UCLI[[k-1]],CLIB[[k]])))
		UCLI
		UCLIS=list(c(1))
		for (k in 2:kk) UCLIS=c(UCLIS,list(intersect(CLIB[[k]],UCLI[[k-1]])))
		UCLIS

		bb=0
		for (k in 2:kk) { 
				b=1
				for (j in 1:(k-1)) {
						b=min(b,length(setdiff(UCLIS[[k]],CLIB[[j]])))
						if (b==0) {
								break
						}
				}
				bb=b+bb
		}
		bb
		if (bb==0) {bbb=1}
		if (bb>0) {bbb=0}
		return(list(bbb=bbb))
}



#' rCLI
#'
#' We draw a random permutation of cliques from CLI and accept it only if it is perfect
#' @param CLIQUES Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' rancli(CLIQUES)
rCLI=function(CLI) {
		UU=0
		while (UU==0)
		{
				nn=length(CLI)
				v=sample(1:nn,nn) # random permutation
				CCC=list()
				for (i in v) CCC=c(CCC,list(CLI[[i]]))
				UU=CHECK1(CCC)$bbb
		}
		return(list(CCC=CCC))
}

#' Wmat
#'
#' FUNCTIONS Wmat AND matmod FORM THE TECHNICAL CORE OF THE ALGORITHM WHICH CONSTRUCTS A PERFECT PERMUTATION OF CLIQUES (ASSUMING THAT IT IS POSSIBLE). function Wmat constructs the weighted coincidence matrix for the set	of cliques CLI with weight w(i,j) = number of vertices in C_i\cap C_j
#' @param A Description of the parameter
#' @keywords A
#' @export
#' @examples
#' Wmat(CLIQUES)

Wmat=function(CLI) {	
	A=matrix(,nrow=length(CLI),ncol=length(CLI))
	for (i in 1:length(CLI)) {
		for (j in 1:length(CLI)) {
			A[i,j]=length(intersect(CLI[[i]],CLI[[j]]))
		}
	}
	
	for (i in 1:length(CLI)) A[i,i]=0
	return(list(A=A))
}

##########################################################################################################
#' matmod
#'
#' FUNCTIONS Wmat AND matmod FORM THE TECHNICAL CORE OF THE ALGORITHM WHICH CONSTRUCTS A PERFECT PERMUTATION OF CLIQUES (ASSUMING THAT IT IS POSSIBLE). function matmod puts zeros in the weighted	incidence matrix A	for those edges connecting vertices from R to vertices from R^c which are removable with respect to R
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' matmod(A,R,CLI)
matmod=function(A,R,CLI=CLI) {
	K1=max(A)
	K2=min(A[A>0])
	
	A3=A
	for (k in (K2+1):K1)
	{
		A1=A3
		A2=A3
		
		for (i in 1:length(CLI)) 
		{
			for (j in 1:length(CLI)) {
				if (A3[i,j]<k & (length(setdiff(c(i,j),R))>0)) A1[i,j]=0 
				if ((A3[i,j]>=k) | (length(setdiff(c(i,j),R))==0)) A2[i,j]=0
			}
		}

		A12=A1
		B<-matrix(1,ncol=length(CLI),nrow=length(CLI))
		for (m in (1:(length(CLI)-1)))
		{
				A12=A12%*%A1
				for (i in R){
						for (j in 1:length(CLI)) {
								if (A2[i,j]*A12[i,j]>0) {
										B[i,j]=0
										B[j,i]=0
								}
						}
				}
		}
		
		A3=A3*B
	}
	
	return(list(W=A3)) 
}


#' RECURRENT.POC.SEARCH
#'
#' THE CORE OF THE RECURSION WHICH PRINTS ALL PERFECT ORDERING OF "CLIQUES" IS FUNCTION: "RECURRENT.POC.SEARCH" THE FUNCTION "RECURRENT.POC.SEARCH.WRAPPER"	 JUST SETS THE STAGE	
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' RECURRENT.POC.SEARCH(...)
RECURRENT.POC.SEARCH <- function(W, W.labels, o.k, CLIQUES, filename)
{
	n=dim(W)[1]
	B=W
	o.k.c=setdiff(1:n, o.k)
	if (length(o.k)>0 & length(o.k.c)<n)
	{
		W=matmod(W,o.k, CLI=CLIQUES)$W # killing the unwanted choices of cliques with respect to o.k 
		B=W[o.k,o.k.c, drop=F] # creates the matrix with columns labels of cliques not-chosen yet
	}
	B.ColSum=colSums(B) # sums of columns in B
	B.labels=W.labels[o.k.c] # labels of cliques which have not been chosen yet
	if (n==2)
	{
			write(c(1,2), file=filename, append=TRUE, sep=",", ncolumns=n) 
			write(c(2,1), file=filename, append=TRUE, sep=",", ncolumns=n)
			return(NULL)
	}
	for (j in B.labels[B.ColSum>0])
	{
			o.k.j=c(o.k,j) # vector of labels of already chosen cliques
			if (length(o.k.j)==n)
			{
					write(o.k.j, file=filename, append=TRUE, sep=",", ncolumns=length(o.k.j)) 
					write(o.k.j[c(2,1,3:length(o.k.j))], file=filename, append=TRUE, sep=",", ncolumns=length(o.k.j)) 
			}
			
			if (length(o.k.j)<2) RECURRENT.POC.SEARCH(W,W.labels, o.k=o.k.j, CLIQUES, filename); # recursion step if o.k.j has single clique
			if (length(o.k.j)>=2 & o.k.j[1] < o.k.j[2]) RECURRENT.POC.SEARCH(W,W.labels, o.k=o.k.j, CLIQUES, filename); # recursion step
		
	}
	return(NULL)
}

#' RECURRENT.POC.SEARCH.WRAPPER
#'
#' THE CORE OF THE RECURSION WHICH PRINTS ALL PERFECT ORDERING OF "CLIQUES" IS FUNCTION: "RECURRENT.POC.SEARCH" THE FUNCTION "RECURRENT.POC.SEARCH.WRAPPER"	 JUST SETS THE STAGE	
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' RECURRENT.POC.SEARCH.WRAPPER(...)
RECURRENT.POC.SEARCH.WRAPPER <- function(CLIQUES, filename='POC_list.csv')
{
	W=Wmat(CLIQUES)$A
	W.labels=col(W)[1,]
	o.k=c()
	if (file.exists(filename)) file.remove(filename)
	out=RECURRENT.POC.SEARCH(W, W.labels, o.k, CLIQUES, filename)
	POCs=matrix()
	SymmGroup.on.C=c()
	if (file.exists(filename))
	{
			POCs=read.table(file=filename, sep=',')
			SymmGroup.on.C=lapply(1:dim(POCs)[1], function(POCi) lapply(POCs[POCi,], function(i) CLIQUES[[i]]))
	}
	return(list(POCs, SymmGroup.on.C))
}

#' get.single.POC
#'
#' get.single.POC descroption
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' get.single.POC(...)
get.single.POC <- function(W, W.labels, o.k, CLIQUES)
{
		n=dim(W)[1]
		B=W
		o.k.c=setdiff(1:n, o.k)
		if (length(o.k)>0 & length(o.k.c)<n)
		{
				W=matmod(W,o.k, CLI=CLIQUES)$W # killing the unwanted choices of cliques with respect to o.k 
				B=W[o.k,o.k.c, drop=F] # creates the matrix with columns labels of cliques not-chosen yet
		}
		B.ColSum=colSums(B) # sums of columns in B
		B.labels=W.labels[o.k.c] # labels of cliques which have not been chosen yet
		for (j in B.labels[B.ColSum>0])
		{
				o.k.j=c(o.k,j) # vector of labels of already chosen cliques
				if (length(o.k.j)==n)
				{
						POC<<-o.k.j
				}
				get.single.POC(W,W.labels, o.k=o.k.j, CLIQUES); # recursion step
		}
}

#' get.single.POC.wrapper
#'
#' get.single.POC descroption
#' @param A Description of the parameter
#' @keywords A
#' @export 
#' @examples
#' get.single.POC.wrapper(...)
get.single.POC.wrapper <- function(CLIQUES)
{
		W=Wmat(CLIQUES)$A
		W.labels=col(W)[1,]
		o.k=c()
		get.single.POC(W, W.labels, o.k, CLIQUES)
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
		return(CHECK1(CLIQUES[POC.candidate])$bbb==1)
}

