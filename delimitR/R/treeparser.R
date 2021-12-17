# This is a function to return a list of matche based on a guide tree.
# There are two helper functions: countmatches() and getfirstmatches()
# The input is a guide tree (e.g. ((1,2),3))
# The output is a list of lists:
# Each outer list encompasses a coalescent interval.
# Each inner list is a coalescent event within that coalescent interval.

#' Generate a list of lists that containes information about coalescent events from a guide tree.
#' @param tree A guide tree.


parsetree <- function(tree){
  index <- 1 # initiate a counter
  listofmatches <- list(list()) # initiate an empty list of lists
  currenttree <- tree # create a variable and store the tree
  currenttree <- gsub(';','',tree) # remove the trailing semicolon
  thecount <- 2 # generate a counter that starts at 2
  while (length(strsplit(currenttree, split="")[[1]]) >5) { # as long as there are more than 5 matches (meaning more than one split left), run this loop
    match_bi <- gregexpr('[[:digit:]],[[:digit:]]', currenttree) # find the number of matches
    thecount <- countmatches(match_bi) # use the helper function, countmatches, to count the number of coalescent events in the current interval.
    myfirstmatches <- getfirstmatches(currenttree,match_bi,thecount) # get the first match
    listofmatches[[index]] <- myfirstmatches # add this to our list
    newtree <- gsub('\\([[:digit:]],([[:digit:]])\\)',"\\1",currenttree ) # replace this first match one of the two species (the second listed in the tree)
    currenttree <- newtree # replace the current tree with this version
    index <- index + 1 # increase the counter
  }
  match_bi <- gregexpr('[[:digit:]],[[:digit:]]', currenttree) # when we're at the last interval, we find the final match
  thecount <- countmatches(match_bi) # we use the helping function to count
  myfirstmatches <- getfirstmatches(currenttree,match_bi,thecount) # we use the helping function to get the match
  listofmatches[[index]] <- myfirstmatches # we add it to the list of matches.
  return(listofmatches)
}

# This function determines how many matches we had when we run gregexpr() to find the lenght of the coalescent interval.
countmatches <- function(myvector){
  value <- "start"
  count <- 1
  while (is.na (value) == FALSE) {
    value = myvector[[1]][count]
    count <- count + 1
  }
  return(count - 2)
}

# This function pulls the matches out of the tree based on indices,
# so we know which species are involved.
getfirstmatches <- function(tree, matchlist, countofmatches){
  treesplit <- strsplit(tree, split="")
  split1list <- list(list())
  for (i in 1:countofmatches) {
    start <- matchlist[[1]][i]
    stop <- matchlist[[1]][i] + 2
    split1list[[i]] <- treesplit[[1]][start:stop]
  }
  return(split1list)
}

