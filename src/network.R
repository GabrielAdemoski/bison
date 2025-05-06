library(SpiecEasi, warn.conflicts=FALSE)
library(Matrix, warn.conflicts=FALSE)
library(igraph, warn.conflicts=FALSE)
library(dplyr, warn.conflicts=FALSE)
library(clustAnalytics, warn.conflicts=FALSE)
library(stringr, warn.conflicts=FALSE)
library(docstring, warn.conflicts=FALSE)


network <- function(group, func = "spieceasi", bin = FALSE, thresh=.5, verbose = TRUE){
  #' Generate optimal network
  #'
  #' @description Run either SpiecEasi or SparCC using a phyloseq object to generate an optimal network. From this
  #' network retrieve the correlation matrix.
  #' @details
  #' corr: A correlation matrix where row names and column names are taxa names and the values of the matrix are OTU
  #' amounts. The maxrix shape is NxN where N is the number of taxa.
  #' lambda: The optimal lambda double value that SpiecEasi / StARS has found from the lambda range given.
  #' stab: The network stability double value. A stability of .05 is recomended by the authors of SpiecEasi. See the
  #' Troubleshooting section of the GitHub page (https://github.com/zdk123/SpiecEasi) for more information.
  #'
  #' @param group. A phyloseq object containing taxa names, OTU amounts, and arbitrary OTU keys.
  #' @param func. A string representing the method of determining the correlation matrix. Currently, only
  #' SpiecEasi and SparCC are implemented.
  #' @param bin. A boolean parameter. When TRUE, weights become either 1 or 0.
  #' @param thresh. A double value and is only used when sparcc is chosen as a function option. The value is
  #' the arbitrary user set threshold for correlations to be considered important.
  #' @param verbose. A boolean parameter. When TRUE, troubleshooting and spiec.easi text will be printed.
  #'
  #' @return A list containing corr, lambda, and stab
  # returns a correlation matrix of shape

  tax = data.frame(tax_table(group))
  last2 = tax[tail(colnames(tax), n=2)]
  # combine the first letter of the second to last column with the final column and remove from the final column
  tax.names = paste0(str_sub(last2[[1]], start=4, end=4), ". ", str_sub(last2[[2]], start=4))

  if (tolower(func) == "spieceasi"){
    # used mb beacuase it gave better results in the original SpiecEasi paper
    se.mb <- spiec.easi(group, method = 'mb', nlambda=1000, lambda.min.ratio=.1,
                          lambda.log=FALSE, verbose=verbose,
                          pulsar.params = list(rep.num = 100, ncores = 16, seed = 42),
                          sel.criterion='stars', pulsar.select=TRUE)
    # The optimal lambda double value that SpiecEasi / StARS has found from the lambda range given.
    lambda = getOptLambda(se.mb)
    # The network stability double value. A stability of .05 is recomended by the authors of SpiecEasi. See the
    # Troubleshooting section of the GitHub page (https://github.com/zdk123/SpiecEasi) for more information.
    stab = getStability(se.mb)

    se.mb.corr = symBeta(getOptBeta(se.mb), mode = "maxabs")
    diag(se.mb.corr) <- 0
    if (bin == TRUE){corr = as.matrix(getRefit(se.mb))}
    else{corr = as.matrix(se.mb.corr)}
    if (verbose == TRUE){
      cat("Optimal lambda: ", lambda, "\n")
      cat("Stability: ", stab, ", .05 is the target. \n")
    }

  }else if(tolower(func) == "sparcc"){
    OTU = t(data.frame(otu_table(group)))

    sparcc.mat <- sparcc(OTU)
    lambda = -1
    stab = -1

    diag(sparcc.mat$corr) <- 0
    if (bin == TRUE){corr = as.matrix(abs(sparcc.mat$Cor) >= thresh)}
    else{corr = as.matrix((abs(sparcc.mat$Cor) >= thresh)*(sparcc.mat$corr))}
  }
  colnames(corr) = tax.names
  rownames(corr) = tax.names

  return (list(corr, stab, lambda))
}

network_stats <- function (corr){
  #' Calculate network statistics
  #'
  #' @description Create an igraph object and calculate network statistics on it.
  #' @details The statistics include:
  #' - number of positive edges
  #' - number of negative edges
  #' - total number of verticies (nodes)
  #' - total number of edges
  #' - average degree of the graph
  #' - average path length of the graph
  #' - clustering coefficient of the graph
  #' - number of lone nodes
  #' - number of clusters
  #' - number of nodes in the largest cluster
  #' - number of nodes in the smallest cluster
  #'
  #' @param corr. A correlation matrix where row names and column names are taxa names and the values of the matrix
  #' are OTU amounts.
  #' @return A list containing the network statistics listed above.


  edge.pos = sum(corr > 0) / 2
  edge.neg = sum(corr < 0) / 2

  # some of the following algorithms require only positive edge weights
  graph = adj2igraph(abs(corr))
  vert = V(graph)
  edge = E(graph)

  vert.total = length(vert)
  edge.total = length(edge)
  deg = degree(graph, vert)
  deg.avg = sum(deg) / length(deg)
  path.len = mean_distance(graph)
  clust.coeff = transitivity(graph, type = "global")


  all.clust = components(graph)
  # delete node with 0 connections
  graph <- delete_vertices(graph, which(degree(graph) < 1))
  trimmed.clust = components(graph)
  lone.nodes = all.clust$no - trimmed.clust$no

  return(list(NumVert = c(vert.total), NumEdges = c(edge.total), AvgDegree = c(deg.avg),
            AvgPathLen = c(path.len), ClusteringCoeff = c(clust.coeff),
            PositiveEdges = c(edge.pos), NegativeEdges = c(edge.neg),
            NumLoneNodes = c(lone.nodes), NumCluster = c(trimmed.clust$no),
            LargestCluster = c(max(trimmed.clust$csize)), SmallestCluster = c(min(trimmed.clust$csize))))

}





