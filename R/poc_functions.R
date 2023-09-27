#' Function lp0
#' 
#' Draws (by rejection) randomly any proper subset of the set 1:n with the same probability.
#' @param A A vector of elements to choose from
#' @param p The probability of including each element in the subset
#' @keywords random subset
#' @export 
#' @return A randomly chosen proper subset of A
#' @examples
#' lp1(1:10)
lp0=function(A, p=0.5)
{
    n=length(A)
    w=rep(1,n)
    while (sum(w)==n) w=rbinom(n,1,p)
    return(A[w==1]) # choose from A only elements at positions in which 1 was generated
}


#' Function lp1
#' 
#' Draws (by rejection) randomly any non-empty subset of the set 1:n with the same probability.
#' @param A A vector of elements to choose from
#' @param p The probability of including each element in the subset
#' @keywords random subset
#' @export 
#' @return A randomly chosen non-empty subset of A
#' @examples
#' lp1(1:10)
lp1=function(A, p=0.5)
{
    n=length(A)
    if (n==1) return(A);
    w=c()
    while (sum(w)==0) w=rbinom(n,1,p);
    return(A[w==1]) # choose from A only elements at positions in which 1 was generated
}

#' Function lp2
#' 
#' Draws (by rejection) randomly any non-empty proper subset of the set 1:n with the same probability.
#' @param A A vector of elements to choose from
#' @param p The probability of including each element in the subset
#' @keywords random subset
#' @export 
#' @return A randomly chosen non-empty proper subset of A
#' @examples
#' lp2(1:10)
lp2=function(A, p = 0.5)
{
    n=length(A)
    if (n==1) return(A);
    w=c()
    while (sum(w)==0 | sum(w)==n) w=rbinom(n,1,p);
    return(A[w==1]) # choose from A only elements at positions in which 1 was generated
}


#' Function split_into_3
#' 
#' Splits a set into 3 non-empty subsets..
#' @param A A vector of elements to choose from
#' @keywords random subset
#' @export 
#' @return A list of non-empty subsets of a set contaning all elements
#' @examples
#' split_into_3(1:10)
split_into_3=function(A)
{
    n=length(A)
	assign_first_elements = sample(1:n,3, replace = FALSE)
	w = sample(c(1,2,3), n, replace=TRUE)
	w[assign_first_elements] = c(1,2,3)
    return(split(A, w)) # choose from A only elements at positions in which 1 was generated
}


#' Function lp3
#' 
#' Draws (by rejection) randomly any subset of size >1 of the set 1:n with the same probability.
#' @param A A vector of elements to choose from
#' @param p The probability of including each element in the subset
#' @keywords random subset
#' @export 
#' @return A randomly chosen subset of size >1 of A
#' @examples
#' split_into_3(1:10)
lp3=function(A, p=0.5)
{
    n=length(A)
    w=c()
    while (sum(w)==1 | sum(w)==0) w=rbinom(n,1,p);
    return(A[w==1]) # choose from A only elements at positions in which 1 was generated
}


#' Function poc1
#' 
#' Draws (randomly, non-uniformly) a perfect ordering of cliques of n nodes (and thus a decomposable graph). Its output CLI is the input to the main procedure.
#' @param n The number of nodes in the graph
#' @param p The probability of including each element in the subsets
#' @param structure The desired graph structure: 'any', 'tree', or 'chain'
#' @keywords random graph, perfect ordering
#' @export 
#' @return A list containing the perfect ordering of cliques and the separators
#' @examples
#' poc1(10, structure = 'tree')
poc1=function(n, p=0.5, structure='any')
{
    B=1:n

    C=lp3(B, p) # draw the first clique
    CLI=list(C)

    B=setdiff(B,C) # what remains
    I=length(B)

    if (I==0)
    {
        SEP=list()
    } # in case the first clique is the whole set the set of separators is empty
    else
    {
        S=lp2(C,p) # otherwise draw the separator S2 from the first clique
        SEP=list(S) # and put it on the list of separators
    }	 
    while (I>0)
    { #	 proceed further only if what remains is non-empty
        R=lp1(B,p) # draw subsequent residual from what remained
        C=sort(union(S,R)) # create the subsequent clique 
        CLI=c(CLI,list(C)) # put in on the list of cliques

        B=setdiff(B,C) # what remains
        I=length(B) # size of what remains
        if (I==0) break  # stop if empty set remains
        ## pick up at random one of the cliques which are on the list at present
        if (structure=='any')
        {
            if (length(CLI)==1)
            {
                L=CLI[[1]]
            }
            else
            {
                L=sample(CLI,1)[[1]]
            }
        }
        if (structure=='tree')
        {
         
            L_for_SEP=lapply(1:length(CLI), function(i) setdiff(CLI[[i]], unlist(SEP)))
            L=c()
            if (length(L_for_SEP)==1)
            {

                L=L_for_SEP
            }
            else
            {
                while (length(L)==0) L=sample(L_for_SEP,1)[[1]];
            }

        }

        if (structure=='chain')
        {
            L_for_SEP=lapply(1:length(CLI), function(i) setdiff(CLI[[i]], unlist(SEP)))
            L=c()
            if (length(L_for_SEP)==1)
            {
                L=L_for_SEP
            }
            else (length(L_for_SEP)==1)
            {
                while (length(L)==0) L=sample(list(L_for_SEP[[1]], L_for_SEP[[length(L_for_SEP)]]),1)[[1]];
            }

        }
        S=lp2(L,p); # draw the subsequent separator from the chosen clique
        SEP=c(SEP,list(S)) # put it on the list of separators
    }
    return(list(CLI=CLI,SEP=SEP))
}

#' is.POC
#' 
#' Checks if CLIQUES is a perfect ordering of cliques
#' @param CLIQUES List of cliques
#' @keywords CLIQUES
#' @export 
#' @examples
#' is.POC(CLIQUES)
is.POC=function(CLIQUES)
{
		K=length(CLIQUES)
		UCLI=lapply(1:K, function(k) unique(unlist(CLIQUES[1:k])))
		UCLIS=lapply(2:K, function(k) intersect(CLIQUES[[k]], UCLI[[k-1]]))

		for (j in 1:(K-1))
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

#' Wmat
#'
#' Constructs the weighted coincidence matrix for the set of cliques CLI with weight w(i,j) = number of vertices in C_i cap C_j
#' @param CLI A list of cliques
#' @keywords CLI
#' @export
#' @examples
#' Wmat(CLI)
Wmat=function(CLI)
{
		K=length(CLI)
		W=array(0,c(K,K))
		for (i in 1:(K-1))
		{
				for (j in (i+1):K)
				{
						W[i,j]=length(intersect(CLI[[i]],CLI[[j]]))
				}
		}
		W=W+t(W)
		return(W)
}
#' matmod
#' 
#' Puts zeros in the weighted incidence matrix A for those edges connecting vertices from R to vertices from R^c which are removable with respect to R
#' @param W The weighted coincidence matrix
#' @param A The set of vertices in R
#' @keywords W, A
#' @export 
#' @examples
#' matmod(W, A)
matmod=function(W,A)
{
		K=dim(W)[1]
		W.not_in.A=array(TRUE, c(K,K)) #initialize as if all edges are not in A
		if (length(A)>1) #Find edges in A and mark with "FALSE"
		{
				all.A.pairs=t(combn(A,2))
				W.not_in.A[all.A.pairs]=FALSE 
				W.not_in.A=W.not_in.A & t(W.not_in.A) #Make it symmetric
		}
		sorted.unique.weights=sort(unique(as.vector(W)))[-1] #Only iterate over actually present weights
		for (w in sorted.unique.weights)
		{
				W.less_w=(W < w)
				W2=W.less_w & W.not_in.A #Removing candidates on level 'k'
				W1=(W2==FALSE) #Not removable
				W1_m=W1
				m=1
				while (m	< (K-1))
				{
						m=m+1
						W1_m=W1_m%*%W1
						W[W2*W1_m>0]=0 #remove edge if there exists a path of length m that goes through non-removable edges on level k
				}
		}
		return(W)
}

#' RCM
#'
#' Prints all perfect ordering of cliques
#' @param W The weighted coincidence matrix
#' @param o.k A vector of the indices of the ordered cliques in the perfect ordering
#' @keywords W, o.k
#' @export 
#' @examples
#' RCM(W, o.k)
RCM <- function(W, o.k)
{
		K=dim(W)[1]
		if (K==2)
		{
				## Using the fact that permutation of initial two elements of POC yields also a POC
				o.k.final=unlist(lapply(strsplit(paste(colnames(W), collapse="_"),"_")[[1]], as.numeric))
				assign("all.POCs", c(all.POCs, o.k.final),envir = .GlobalEnv) 
				assign("all.POCs", c(all.POCs, o.k.final[c(2,1,3:length(o.k.final))]),envir = .GlobalEnv)
		}
		## Using the fact that permutation of initial two elements of POC yields also a POC
		B=W
		B.labels=colnames(W)
		if (length(o.k)>0)
		{
				W=matmod(W,o.k) # killing the unwanted choices of cliques with respect to o.k 
				B=W[o.k,-o.k, drop=F] # creates the matrix with columns labels of cliques not-chosen yet
				B.labels=B.labels[-o.k, drop=F] # labels of cliques which have not been chosen yet
		}
		if (length(o.k)>=2)
		{
				collapsed_weights=apply(W[o.k, -o.k, drop=F],2, max)
				collapsed_names = c(paste(colnames(W)[o.k], collapse=""),colnames(W)[-o.k])
				W=rbind(c(0, collapsed_weights), cbind( collapsed_weights, W[-o.k, -o.k, drop=F]))
				colnames(W) = collapsed_names
				rownames(W) = collapsed_names
				o.k=c(1)
		}
		B.ColSum=colSums(B) # sums of columns in B
		for (j in B.labels[B.ColSum>0])
		{
				o.k.j=c(o.k,which(colnames(W)==j)) # vector of labels of already chosen cliques inlcuding possibly collapsed clique
				if (length(o.k.j)<2 & K >1) RCM(W, o.k=o.k.j); # RCM step if o.k.j has single clique
				if (length(o.k.j)>1 & o.k.j[1] < o.k.j[2] & K>1) RCM(W, o.k=o.k.j); # RCM step. ensuring we will not double the outcomes with investigating only one order of initial two elements of a POC
		}
}

#' RCM.WRAPPER
#'
#' The core of the RCM which prints all perfect orderings of "cliques" is function: "RCM". The function "RCM.WRAPPER" sets the stage.
#' 
#' @param CLIQUES a list of cliques 
#' @return a matrix of all perfect orderings of cliques
#' @keywords RCM, wrapper function
#' @export 
#' @examples
#' RCM.WRAPPER(list(1:3, 3:5, 2:4))
RCM.WRAPPER <- function(CLIQUES)
{
		assign("all.POCs", c(), envir = .GlobalEnv)
    	assign("all.PROBSs", c(), envir = .GlobalEnv)

		K=length(CLIQUES)
		if (K==1) return(matrix(c(1),1,1));
		if (K==2) return(matrix(c(1,2,2,1),2,2));
		W=Wmat(CLIQUES)
		colnames(W)=1:dim(W)[1]
		out=RCM(W=W, o.k=c())
		POCs=matrix(all.POCs, ncol=K, byrow=TRUE)

		return(POCs)
}

#' get.all.POCs
#'
#' Gets all possible perfect orderings of cliques using RCM or rejection sampling method.
#' 
#' @param CLIQUES a list of cliques 
#' @param method the sampling method to use (either "RCM" or "rejection")
#' @return a matrix of all perfect orderings of cliques
#' @keywords RCM, rejection sampling, cliques
#' @export 
#' @examples
#' get.all.POCs(list(1:3, 3:5, 2:4), "RCM")
get.all.POCs <- function(CLIQUES, method)
{
		K=length(CLIQUES)
		if (method=='RCM')
		{
				POCs=RCM.WRAPPER(CLIQUES)
		}
		if (method=='rejection')
		{
				POCs=c()
				all.permutations=permn(1:K)
				if (length(all.permutations)==2) POCs=list(c(1,2), c(2,1))
				else
						{
								for (i in 1:factorial(K))
								{
										perm.i=all.permutations[[i]]
										if (perm.i[1] < perm.i[2] & length(perm.i)>2)
										{
												if (is.POC(CLIQUES[perm.i])) POCs=rbind(POCs,perm.i, perm.i[c(2,1,3:length(perm.i))]);
										}
								}
						}
		}
		return(POCs)
}


#' single.RCM
#'
#' The single RCM function uses recursive coalescence modeling to compute one perfect ordering of cliques.
#' 
#' @param W a matrix of weights
#' @param o.k a vector of indices representing cliques already chosen
#' @param prob the probability of the perfect ordering
#' @return a list with two items: a vector representing a perfect ordering of cliques and a probability
#' @keywords RCM, cliques, probability
#' @export 
single.RCM <- function(W, o.k,prob)
{
    K=dim(W)[1]
		if (K==2)
		{
      ## Using the fact that permutation of initial two elements of POC yields also a POC
      poc=unlist(lapply(strsplit(paste(colnames(W), collapse="_"),"_")[[1]], as.numeric))
      ##assign("POC",poc, envir=.GlobalEnv)
      ##assign("PROB",c(prod(prob),prob), envir=.GlobalEnv)
      return(list(poc, c(prod(prob), prob)))
		}

    		## Using the fact that permutation of initial two elements of POC yields also a POC
		B=W
		B.labels=colnames(W)
		if (length(o.k)>0)
		{
				W=matmod(W,o.k) # killing the unwanted choices of cliques with respect to o.k 
				B=W[o.k,-o.k, drop=F] # creates the matrix with columns labels of cliques not-chosen yet
				B.labels=B.labels[-o.k, drop=F] # labels of cliques which have not been chosen yet
		}
		if (length(o.k)>=2)
		{
				collapsed_weights=apply(W[o.k, -o.k, drop=F],2, max)
				collapsed_names = c(paste(colnames(W)[o.k], collapse="_"),colnames(W)[-o.k])
				W=rbind(c(0, collapsed_weights), cbind( collapsed_weights, W[-o.k, -o.k, drop=F]))
				colnames(W) = collapsed_names
				rownames(W) = collapsed_names
				o.k=c(1)
		}
	B.ColSum=colSums(B) # sums of columns in B
    B0.len=length(B.labels[B.ColSum>0])
	
    if (B0.len>= 1)
    {
        j=sample(B.labels[B.ColSum>0],1)
        prob = c(prob, 1/B0.len)
        o.k.j=c(o.k,which(colnames(W)==j)) # vector of labels of already chosen cliques inlcuding possibly collapsed clique

        if (length(o.k.j)<2 & K >1) return(single.RCM(W, o.k=o.k.j, prob=prob)); # RCM step if o.k.j has single clique
        if (length(o.k.j)>1 & K>1) return(single.RCM(W, o.k=o.k.j, prob=prob)); # RCM step. ensuring we will not double the outcomes with investigating only one order of initial two elements of a POC
    }
		}

#' single.RCM.WRAPPER
#'
#' The function single.RCM.WRAPPER returns a perfect ordering of cliques (POC) of the list of cliques.
#' @param CLIQUES A list of cliques.
#' @return Returns a list containing the perfect ordering of cliques (POC) and the probability of the POC.
#' @examples
#' single.RCM.WRAPPER(CLIQUES)
single.RCM.WRAPPER <- function(CLIQUES)
{
  ##assign("POC", c(), envir = .GlobalEnv)
  ##assign("PROB", 1, envir = .GlobalEnv)
  K=length(CLIQUES)
  if (K==1)
  {
	 return(c(c(1),1))
  } 
  if (K==2) 
  {
	return(list(list(c(1,2),c(2,1))[rbinom(1,1,prob=0.5)+1],0.5))
  }
  if (K>2)
  {
	W=Wmat(CLIQUES)
	colnames(W)=1:dim(W)[1]
	poc.prob = single.RCM(W=W, o.k=c(), prob=c())
	return(list(poc.prob[[1]], poc.prob[[2]]))
  }

}

#' rPOC_RCM
#'
#' The function rPOC_RCM generates random samples of perfect orderings of cliques (POCs) from the list of cliques using the RCM method.
#' @param nsamples Number of POCs to be generated.
#' @param CLIQUES A list of cliques.
#' @param replace Boolean indicating whether to sample with or without replacement.
#' @return Returns a list of generated POCs.
#' @examples
#' rPOC_RCM(100, CLIQUES)
rPOC_RCM <- function(nsamples,CLIQUES,method='RCM', replace=TRUE)
{
		nPOCs=list()
		if (method=='RCM')
		{
				if (replace)
				{
						nPOCs=replicate(nsamples,single.RCM.WRAPPER(CLIQUES)[[1]])
				}
				if (!(replace))
				{
						if (nsamples>factorial(length(CLIQUES)))
						{
								return(nPOCs)
						}
						counter=1
						while (counter<=nsamples)
						{
								POC.candidate=single.RCM.WRAPPER(CLIQUES)[[1]]
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
#' The function rPOC_unif generates random samples of perfect orderings of cliques (POCs) from the list of cliques using the rejection method.
#' @param n Number of POCs to be generated.
#' @param CLIQUES A list of cliques.
#' @param method Method to be used for sampling. Either 'RCM' or 'rejection'.
#' @param replace Boolean indicating whether to sample with or without replacement.
#' @param mh_burn Integer indicating the number of Metropolis-Hastings burn-in iterations to be used if method='M-H'.
#' @return Returns a list of generated POCs.
#' @examples
#' rPOC_unif(100, CLIQUES)
rPOC_unif<-function(n,CLIQUES,method='RCM', replace=TRUE, mh_burn=100)
{
		## it should return an array of POCs of size (n, len(Cliques))
		K=length(CLIQUES)
		nPOCs=array(0,c(n, K))
		if (method=='rejection')
		{
				if (replace)
				{
						counter=1
						while (counter<=n)
						{
								POC.candidate=sample(1:K)
								if (is.POC(CLIQUES[POC.candidate]))
								{
										nPOCs[counter,]=POC.candidate
										counter=counter+1
								}
						}
				}
				if (!(replace))
				{
						if (n>factorial(K))
						{
								print('Warning: sampling space smaller than sample size')
								return(nPOCs)
						}
						all.permutations=permn(1:K)
						counter=1
						while (counter<=n)
						{
								POC.candidate=unlist(sample(all.permutations,1))
								if (!(list(POC.candidate) %in% nPOCs) && is.POC(CLIQUES[POC.candidate]))
								{
										nPOCs[counter,]=POC.candidate
										counter=counter+1
								}
						}

				}

		}
		if (method=='RCM')
		{
				POCs=RCM.WRAPPER(CLIQUES)
				if (replace)
				{
						nPOCs=POCs[sample(1:dim(POCs)[1],n, replace=TRUE),,drop=FALSE]
				}
				if (!(replace))
				{
						nPOCs=POCs[sample(1:dim(POCs)[1],n, replace=FALSE),,drop=FALSE]
				}
		}

		if (method=="M-H")
		{
				if (replace)
				{
						nPOCs=t(replicate(n,get.MH.poc(CLIQUES,mh_burn)))

				}
				if(!(replace))
				{
						while (counter<=n)
				{
						POC.candidate=get.MH.poc(CLIQUES,mh_burn)
						if (!(list(POC.candidate) %in% nPOCs))
						{
								nPOCs[counter,]=POC.candidate
								counter=counter+1
						}
				}
				}

		}

		return(nPOCs)
}

#' get.MH.poc
#'
#' The function get.MH.poc generates a perfect ordering of cliques (POC) using the Metropolis-Hastings method.
#' @param CLIQUES A list of cliques.
#' @param mh_burn Integer indicating the number of Metropolis-Hastings burn-in iterations to be used.
#' @return Returns the generated POC.
#' @examples
#' get.MH.poc(CLIQUES,1000)
get.MH.poc<-function(CLIQUES, mh_burn=10)
{
		sigma_and_prob=single.RCM.WRAPPER(CLIQUES)
		sigma=sigma_and_prob[[1]]
		prob_sigma=sigma_and_prob[[2]][1]
		for (i in 1:mh_burn)
		{
				tau_and_prob=single.RCM.WRAPPER(CLIQUES)
				tau=tau_and_prob[[1]]
				prob_tau=tau_and_prob[[2]]
				alpha=min(1, prob_sigma/prob_tau)
				sigma=list(sigma, tau)[[sample(1:2, 1,prob=c(1-alpha, alpha))]]
		}
		return(sigma)
}


#' is.decomposable
#'
#' The function is.decomposable determines if a list of cliques is decomposable by checking if a perfect ordering of cliques (POC) can be generated using the provided method.
#' @param CLIQUES A list of cliques.
#' @param method Method to be used for generating POCs. Either 'RCM' or 'rejection'.
#' @return Returns a boolean indicating if the list of cliques is decomposable.
#' @examples
#' is.decomposable(CLIQUES, method='RCM')
is.decomposable <- function(CLIQUES, method='RCM')
{
  K=length(CLIQUES)
  found_POC=FALSE
  if (method=='RCM')
  {
    while (!found_POC)
  {
    found_POC=is.POC(CLIQUES[rPOC_RCM(1, CLIQUES)])
  }
  }
  if (method=="rejection")
  {
    while(!found_POC)
  {
    found_POC=is.POC(CLIQUES[sample(1:K)])
  }
  }
  return(found_POC)
}
