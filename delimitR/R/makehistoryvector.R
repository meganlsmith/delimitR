#' history vector
#' @keywords internal


# This takes the output from renaming species and creates a list of historic events.
makehistoryvector <- function(parseoutput){
  myhistevents <- list(list()) # create an empty list to store the history
  indexhist <- 1 # initiate a counter at 1
  for (i in 1:length(parseoutput)) { # for each coalescent interval
    thegroup <- parseoutput[[i]] # store the part of the parse output list for this interval
    for (j in 1:length(thegroup)) { # for each coalescent event in this interval
      theevent <- thegroup[[j]] # store the coalescent event
      itemtoadd <- paste(theevent[1],theevent[3]) # get the taxa indices from this
      itemtoadd <- gsub('\\"', "", itemtoadd) # get rid of quotation marks in this
      myhistevents[[indexhist]] <- itemtoadd # add the historic event
      indexhist <- indexhist + 1 # add one to our counter
    }
  }
  return(myhistevents) # return a list of events
}
