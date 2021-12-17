#' Generate history
#' @keywords internal
# This function formats coalescent events as they should be for fsc2.
gethistory <- function(nhistevents, myhistory){
  thehistory <- c() # empty vector
  for (j in 1:nhistevents) { # for each event
    events <- myhistory[[j]] # assign to 'events' that item from my history
    eventssplit <- strsplit(events, " ") # split into a list of populaiton indices
    thehistory[j] <- paste('Tdiv', j, '$', ' ', eventssplit[[1]][1], ' ', eventssplit[[1]][2], ' ', 1, ' ', 1, ' ', 0, ' ', 'keep', sep = '') # create the fsc2-style divergence event
  }
  return(thehistory)
}

# this generates a grid with coalescent and collapse events
geteventgrid <- function(thehistory){
  # generate all binary combos, 0=collapse, 1=don't collapse
  eventgrid <- expand.grid(rep(list(0:1),length(thehistory)))
  isokay <- data.frame(matrix(TRUE, nrow = nrow(eventgrid), ncol = ncol(eventgrid)))
  for (i in 1:nrow(eventgrid)) { ## loop through each potential combo
    for (j in 1:ncol(eventgrid)) { ## go through and check each column, to see if rules are broken
      if (j > 1 & eventgrid[i,j] != 1) {
        pop1 <- strsplit(thehistory[[j]], split=" ")[[1]][2]
        pop2 <- strsplit(thehistory[[j]], split=" ")[[1]][3]
        klistcounts <- list()
        acount = 0
        for (k in 2:j-1) {
          acount = as.numeric(strsplit(thehistory[[k]], split=" ")[[1]][2]==pop1)
          acount = acount+ as.numeric(strsplit(thehistory[[k]], split=" ")[[1]][3]==pop1)
          acount = acount+ as.numeric(strsplit(thehistory[[k]], split=" ")[[1]][2]==pop2)
          acount = acount+ as.numeric(strsplit(thehistory[[k]], split=" ")[[1]][3]==pop2)
          klistcounts[k] <- acount
        }
        if (length(which(klistcounts==0))!=length(klistcounts)) {
          for (k in 1:(j-1)) {
            pop1 <- strsplit(thehistory[[j]], split=" ")[[1]][2]
            pop2 <- strsplit(thehistory[[j]], split=" ")[[1]][3]
            acount1 <- as.numeric(strsplit(thehistory[[k]], split=" ")[[1]][2] == pop1)
            acount2 <- as.numeric(strsplit(thehistory[[k]], split=" ")[[1]][3] == pop1)
            acount3 <- as.numeric(strsplit(thehistory[[k]], split=" ")[[1]][2] == pop2)
            acount4 <- as.numeric(strsplit(thehistory[[k]], split=" ")[[1]][3] == pop2)
            if (max(acount1, acount2, acount3,acount4)>0) {
                if (eventgrid[i,k] == 1) {
                isokay[i,j] = 'FALSE'
              }
            }
          }
        }
      }
    }
  }
  NAvec <- which(isokay==FALSE, arr.ind=TRUE)
  rowstodrop <- NAvec[,1]
  if(length(rowstodrop) > 0) {
    mycollapses <- eventgrid[-c(rowstodrop), ]
  }
  else {
    mycollapses <- eventgrid
  }
  return(mycollapses)
}

## now we need to generate lists of history events with all possible combos of collapses
historieswithcollapse <- function(eventgrid, thehistory){
  historywithcollapse <- list()
  for (i in 1:nrow(eventgrid)) {
    thishistory <- list()
    for (j in 1:ncol(eventgrid)) {
      if (eventgrid[i,j] == 0) {
        search <- paste('Tdiv', j, '\\$', sep = "")
        store <- thehistory[[j]]
        thishistory[j] <- gsub(search, 0, store)
      }
      else {
        thishistory[j] <- thehistory[[j]]
      }
    }
    historywithcollapse <- append(historywithcollapse, list(thishistory))
  }
  return(historywithcollapse)
}
