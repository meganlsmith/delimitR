#' Simulate data in fastsimcoal2
#' @export
#' @param prefix The prefix for the .tpl and .est files in your working directory.
#' @param pathtofsc The path to fastsimcoal2.
#' @param nreps The number of datasets to simulate under each model.
#' @return This will simulate data in directories created in your working directory.

## now it's time to do our fastsimcoal simulations

fastsimcoalsims <- function(prefix,pathtofsc,nreps){
  listoftpl <- list()
  listofest <- list()
  tpllist <- system(paste("ls ", prefix,"*.tpl",sep=""),intern=T)
  estlist <- system(paste("ls ", prefix,"*.est",sep=""),intern=T)
  listoftpl <- c(listoftpl, tpllist)
  listofest <- c(listofest, estlist)
  count <- 1
  for(j in 1:length(listoftpl)){
    print(paste(pathtofsc, " -t ", prefix, "_", count, ".tpl",  " -e ", prefix, "_", count, ".est", " -n 1 --msfs -q --multiSFS -x -E" ,nreps, sep = ""))
    system(paste(pathtofsc, " -t ", prefix, "_", count, ".tpl",  " -e ", prefix, "_", count, ".est", " -n 1 --msfs -q --multiSFS -x -E" ,nreps, sep = ""), ignore.stdout = TRUE)
    count <- count + 1
  }
}
