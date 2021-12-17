#' Remove unnecessary files from the working directory
#' @export
#' @param prefix The prefix used in the delimitR analysis you wish to delete files from.
clean_working <- function(prefix){
  listoffolders <- list()
  folderlist <- system(paste("ls -d ", prefix,"*/",sep=""),intern=T)
  listoffolders <- c(listoffolders,folderlist)
  for( folder in listoffolders){
    system(paste("rm -r ", folder, sep = ""))
  }
  system("rm *.params")
  system(paste("rm ", " *Processed_*.obs", sep = ""))
  system(paste("rm ", prefix, "*_MSFS.obs", sep = ""))
}

