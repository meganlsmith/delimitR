### This is a function that checks the input provided by the user for several things:
### 1) an incorrectly specified tree file
### 2) An inaccurate number of species
### 3) The wrong number of population size, sample size, or divergence times given in the prior.

#' Check the user input for issues
#' @param tree A guide tree specified in Newick format. Tip labels should be 0 to (nspec-1).
#' @param samplesizes A vector of sample sizes.
#' @param nspec The number of potential species.
#' @param popsizeprior Prior for population sizes.
#' @param divtimeprior Prior on divergence times.

checkinput <- function(tree, nspec, samplesizes, popsizeprior,divtimeprior){
  store <- gregexpr('[0-9]', tree, perl = TRUE)[[1]]
  storesub <- gsub(',','',tree)
  storesub <- gsub(')','',storesub)
  storesub <- gsub('\\(','',storesub)
  storesub <- gsub(';','',storesub)
  storesub <- strsplit(storesub,'')
  storesub2 <- lapply(storesub, function(x) as.numeric(x))
  storesub2 <- c(storesub2)
  storesub <- sort(storesub2[[1]], decreasing = FALSE)
  storesub <- paste(storesub, collapse = "")
  shouldname <- paste(0:(nspec-1), collapse = "")
  if (length(store)!= nspec){
    stop("The number of species in the guide tree does not match the provided number of species.")
  }
  if (length(samplesizes) != nspec){
    stop("The number of sample sizes provided does not match the provided number of species.")
  }
  if (length(popsizeprior) != nspec){
    stop("The number of population size priors provided does not match the provided number of species.")
  }
  if (length(divtimeprior) != nspec-1){
    stop("The number of divergence time priors provided is not consistent the provided number of species.")
  }
  if (storesub != shouldname){
    stop("The guide tree was specified incorrectly. Make sure you are following the naming conventions discussed in the tutorial.")
  }
  if(length(grep(';',tree))==0){
    stop("The guide tree was specified incorrectly. Ensure that the tree is in Newick format.")
  }
}


