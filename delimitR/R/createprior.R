#' Build a prior
#' @export
#' @import utils
#' @param prefix The prefix for the .tpl and .est files in your working directory.
#' @param nspec The number of putative species.
#' @param nclasses The number of bins to use in the SFS.
#' @param mydir The working directory for the analysis.
#' @param traitsfile A beast-formatted traits file.
#' @param threshold If a downsampling approach was used to create the observed SFS, provide the threshold here.
#' @param ncores The number of cores available for the analysis.
#' @param thefolder A name for the folder to save the prior in.
#' @return This will return the full prior to use when building the RF classifier.


## organize and bin
makeprior <- function(prefix, nspec,nclasses,mydir,traitsfile,threshold, thefolder,ncores){
  listoffolders <- list()
  folderlist <- system(paste("ls -d ", prefix,"*/",sep=""),intern=T)
  listoffolders <- c(listoffolders,folderlist)
  #print(listoffolders)
  for(j in 1:length(listoffolders)){
    #print(paste('cp ',listoffolders[[j]],'*.obs ', mydir, sep=""))
    system(paste('cp ',listoffolders[[j]],'*.obs ', mydir, sep=""))
  }

  #  for(i in 1:length(listoffolders)){
  #    system(paste("awk 'NR == 0 || NR % 3 == 0' ",gsub("/", "", listoffolders[[i]]),"_MSFS.obs > ",gsub("/", "", listoffolders[[i]]),"_processed.obs",sep=""))
  #  }
  maximum <- length(listoffolders)
  for(i in 1:maximum){
    listoffolders[[i]] <- gsub("/", "_MSFS.obs", listoffolders[[i]])
  }
  parallel::mclapply((listoffolders),
                     function(fn) system(paste("awk 'NR == 0 || NR % 3 == 0' ",fn," > ","Processed_",fn,sep="")),
                     mc.cores=ncores, mc.preschedule=F)
  for(i in 1:maximum){
    listoffolders[[i]] <- gsub("^", "Processed_", listoffolders[[i]])
  }
  parallel::mclapply((listoffolders),
                     function(fn) system(paste('python ',system.file(package="delimitR"), '/', 'SFS_CreateBinned_', nspec, 'pops.py ',traitsfile,' ', fn, " ",threshold,' ',nclasses,' ','Binned_',fn,sep="")),
                     mc.cores=ncores, mc.preschedule=F)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  system(paste('mkdir', thefolder, sep = " "))
  for(i in 1:maximum){
    listoffolders[[i]] <- gsub("^", "Binned_", listoffolders[[i]])
  }
  for(j in 1:length(listoffolders)){
    system(paste('cp ',listoffolders[[j]],' ', thefolder, sep=""))
  }

  filelist <- list.files(thefolder)
  newlist <- c()
  for(filenum in 1:length(filelist)){
    newlist<- c(newlist,paste('Binned_Processed_',prefix,'_',filenum,'_MSFS.obs',sep=""))
  }
  filelist <- newlist
  filenames <- list()
  varlist <- list()
  for(i in 1:length(filelist)){
    filename <- filelist[i]
#    print(filename)
    varname <- gsub(".obs","",filename)
    filenames <- c(filenames, filename)
    varlist <- c(varlist, varname)
    assign(gsub(".obs","",filename), utils::read.csv(filename,sep="\t", colClasses=c(rep("numeric",(nclasses**nspec)),rep("NULL",1)),header=F))
  }
  FullPrior <- data.frame()
  for(i in 1: length(varlist)){
    test <- as.name(varlist[[i]])
    test
    tester <- eval(test)
    tester$Model <- as.factor(gsub("Binned_Processed_","",varlist[[i]]))
    tester
    assign(varlist[[i]],tester)
    FullPrior <- rbind(FullPrior,tester)
  }
  return(FullPrior)
}


