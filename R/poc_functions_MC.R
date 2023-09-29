#' shuffle_lists
#'
#' Merge and shuffle lists while preserving order of objects within each list
#'
#' This function takes a list of lists and returns a new list containing all objects from the input lists,
#' with the order of the objects shuffled, but preserving the order of objects within each original list.
#'
#' @param lists A list of lists containing the objects to be merged and shuffled.
#' @return A new list containing all objects from the input lists, shuffled but with the order of objects preserved
#' within each original list.
#' @examples
#' list1 <- list("a", "b", "c")
#' list2 <- list("d", "e")
#' shuffled_lists(list(list1, list2))
  shuffle_lists <- function(lists)
  {
  # Get counts on various levels
  components_count <- length(lists)
  components_sizes <- lengths(lists)
  K <- sum(components_sizes)
  proba = prod(sapply(components_sizes, function(Ki) factorial(Ki)))/factorial(K)
# Initialize a new empty list to hold the merged and shuffled objects
  shuffled_list <- vector("list", length = K)
  available_indices <- seq_along(shuffled_list)

# For each component list, shuffle it into the shuffled_list, preserving the order of appearence in original list
  for (i in seq_len(components_count)) 
  {
    component <- lists[[i]]
    n <- components_sizes[i]
    if (length(available_indices)<=1)
    {
      shuffled_list[available_indices] = component
      break()
    }
    if (length(available_indices)>1)
    {
    indices <- sort(sample(c(available_indices), n, replace = FALSE))
    shuffled_list[indices] <- component
    available_indices <- c(setdiff(available_indices, indices))
    }
  }

# Return the merged and shuffled list
  return(list(shuffled_list, proba))
# add proba - ration. K1!*...*Kn! n  -skladowe / (sum(Ki))!
}


#' get_CLI_next_length
#'
#' Generates a random length for the next clique to be added to a graph.
#'
#' This function takes the remaining vertices and the size of the current clique and returns the length of the next clique to be added to the graph. The length of the next clique is chosen randomly, either K-1, K, or K+1, with equal probaability.
#'
#' @param rmng_V A vector containing the vertices that have not yet been added to the graph.
#' @param K The size of the current clique.
#' @return A numeric value indicating the length of the next clique to be added.
#' @examples
#' get_CLI_next_length(c(1, 2, 3), 4)
#' @export
get_CLI_next_length <- function(rmng_V, K)
{
  if (length(rmng_V) >=3)
  {
    L = sample(c(K-1, K, K+1),1)
    proba=1/3

  }
  else
  {
    L = sample(c(K-1, K),1)
    proba=1/2
  }
  return(list(L, proba))
}

#' markov_chain
#'
#' Generates a Markov chain of graphs using the POC algorithm.
#'
#' This function takes an initial graph and generates a Markov chain of graphs using the POC algorithm. The maximum number of iterations is specified by the max_iter argument. The output is a list of graphs generated during the Markov chain.
#'
#' @param CLI_init The initial graph, represented as a list of cliques.
#' @param max_iter The maximum number of iterations for the Markov chain. Default is 100.
#' @return A list of graphs generated during the Markov chain.
#' @examples
#' markov_chain(list(c(1, 2), c(3, 4)))
#' @export
markov_chain <- function(CLI_init, max_iter = 100, poc_shuffle = TRUE)
{
  markov_chain = list()
  for (i in 1:max_iter)
  {
    K = length(CLI_init)
    markov_chain[[i]]=CLI_init
    CLI_next = POC_jump(CLI_init, poc_shuffle = poc_shuffle)
    N = length(CLI_next)

    CLI_init = CLI_next
  }
  return(markov_chain)
}


#' POC_jump
#'
#' Generates a new graph by adding a new clique to an existing graph using the POC algorithm.
#'
#' This function takes an existing graph and generates a new graph by adding a new clique using the POC algorithm. The probaability of adding each vertex to the new clique is specified by the p argument. The output is a list of cliques representing the new graph.
#'
#' @param CLI_init The initial graph, represented as a list of cliques.
#' @param p The probaability of adding each vertex to the new clique. Default is 0.5.
#' @param poc_shuffle Whether to shuffle the vertices before adding the new clique. Default is TRUE.
#' @return A list of cliques representing the new graph.
#' @examples
#' POC_jump(list(c(1, 2), c(2, 3)), p = 0.8)
#' @export
POC_jump = function(CLI_init, p=0.5, poc_shuffle = TRUE)
{
  proba = 1
  all_vertices = unique(unlist(CLI_init))
  n=length(all_vertices)
  K = length(CLI_init)
  
  if (K>=3)
  {
    CLI_next = CLI_init[1:(K-2)] # remove the last two cliques

    rmng_V = setdiff(all_vertices,unique(unlist(CLI_next)))
    L.proba = get_CLI_next_length(rmng_V, K)
    L = L.proba[[1]]
    proba = proba*L.proba[[2]]
    if (L==(K-1))
    {
      separator_clique = sample(CLI_next,1)
      proba = proba*1/length(CLI_next)

      sep = lp0(unlist(separator_clique), p)
      proba = proba*1/(2^length(unlist(separator_clique))-1)

      CLI_next = c(CLI_next, list(union(sep, rmng_V)))
    }
    if (L==K)
    {
      separator_clique = sample(CLI_next,1)
      proba = proba*1/length(CLI_next)
      
      sep = lp0(unlist(separator_clique),p)
      proba = proba*1/(2^length(unlist(separator_clique))-1)

      add_V = lp2(rmng_V,p)
      proba = proba*1/(2^length(rmng_V)-2)

      rmng_V = setdiff(rmng_V, add_V)
      CLI_next = c(CLI_next, list(union(sep, add_V)))
      separator_clique = sample(CLI_next,1)
      proba = proba*1/length(CLI_next)

      sep = lp0(unlist(separator_clique),p)
      proba = proba*1/(2^length(unlist(separator_clique))-1)
      
      CLI_next = c(CLI_next, list(union(sep, rmng_V)))
    }
    if (L==(K+1))
    {
      rmng_V_list = split_into_3(rmng_V)
      # Stirling numbers
      proba = proba*1/Stirling2(length(rmng_V),3)
      for (i in 1:length(rmng_V_list)) #This is always
        {
          rmng_V_current = rmng_V_list[[i]]

          separator_clique = sample(CLI_next,1)
          proba = proba*1/length(CLI_next)

          sep = lp0(unlist(separator_clique),p)
          proba = proba*1/(2^length(unlist(separator_clique)) -1)

          CLI_next = c(CLI_next, list(union(sep, rmng_V_current)))
        }

    }
  }
  if (K==2 && n>=3)
  {
    L.proba = get_CLI_next_length(all_vertices, K)
    L = L.proba[[1]]
    proba = proba*L.proba[[2]]
    if (L==1) CLI_next = list(all_vertices);

    if (L==2)
    {
      CLI_next = list(lp2(all_vertices,p))
      proba = proba*1/(2^length(all_vertices) -2)
      rmng_V = setdiff(all_vertices, unlist(CLI_next))
      sep = lp0(unlist(CLI_next),p)
      proba = proba*1/(2^length(unlist(CLI_next)) -1)

      CLI_next = c(CLI_next, list(union(sep, rmng_V)))
    }
    if (L==3)
    {
      rmng_V_list = split_into_3(all_vertices)
      print(rmng_V_list)
      # Stirling numbers
      proba = proba*1/Stirling2(length(all_vertices),3)
      CLI_next = rmng_V_list[1]
      print(CLI_next)
      for (i in 2:length(rmng_V_list))
      {
        rmng_V_current = rmng_V_list[[i]]
        separator_clique = sample(CLI_next,1)
        proba = proba*1/length(CLI_next)
        sep = lp0(unlist(separator_clique),p)
        proba = proba*1/(2^length(unlist(separator_clique)) -1)
        
        CLI_next = c(CLI_next, list(union(sep, rmng_V_current)))
      }
    }
  }
  if ((K==2 && n==2) || K==1)
  {
    L = sample(c(1,2),1)
    proba = proba*1/2

    if (L==1) CLI_next = list(all_vertices);
    if (L==2)
    {

      CLI_next = list(lp2(all_vertices,p))
      proba = proba*1/(2^length(all_vertices) -2)
      rmng_V = setdiff(all_vertices, unlist(CLI_next))
      sep = lp0(unlist(CLI_next),p)
      proba = proba*1/(2^length(unlist(CLI_next)) -1)
      
      CLI_next = c(CLI_next, list(union(unlist(sep), unlist(rmng_V))))
    }
  }
 # Depending on when in Markov Chain step we do the shuffle, we obtain potentially different MC
  if (poc_shuffle) 
  {
    if (length(CLI_next)==1)
    {
      CLI_next = CLI_next
    }
    else
    {
      shuffle_proba=1
      new_poc_proba = list()
      # Old shuffle option (single connected components)
      # new_poc_order = rPOC_RCM(1,CLI_next, replace=TRUE)
      # CLI_next = CLI_next[new_poc_order]
      
      # New shuffle option (multiple connected components):
      g = Get.weighted.clique.graph(CLI_next)
      cl <- components(g) # get connected components data
      if (cl$no == 1)
      {
        new_poc_order.proba = single.RCM.WRAPPER(CLI_next)
        new_poc_order  = new_poc_order.proba[[1]]
        new_poc_proba[[1]] = new_poc_order.proba[[2]][1]
        CLI_next = CLI_next[unlist(new_poc_order)]
      }
      if (cl$no > 1)
      {
        CLI_nexts = lapply(1:cl$no, function(comp_id) CLI_next[which(cl$membership == comp_id)])
        
        for (i in 1:length(CLI_nexts))
        {
          CLI_component = CLI_nexts[[i]]
          new_poc_order.proba = single.RCM.WRAPPER(CLI_component)

          new_poc_order  = new_poc_order.proba[[1]]
          new_poc_proba[[i]]  = new_poc_order.proba[[2]][1]
          CLI_nexts[[i]] = CLI_component[unlist(new_poc_order)]
        }
        #CLI_nexts_new_order = lapply(CLI_nexts, function(CLI_component) CLI_component[unlist(rPOC_RCM(1,CLI_component, replace=TRUE))])
        CLI_next.proba = shuffle_lists(CLI_nexts)
        CLI_next = CLI_next.proba[[1]]
        shuffle_proba = CLI_next.proba[[2]]
      }
      proba = shuffle_proba*proba*prod(unlist(new_poc_proba))
    }
  }
  print(proba)
  return(list(CLI_next, proba))
}
