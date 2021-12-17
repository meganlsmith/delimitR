#' Generate input files for fastsimcoal2
#' @export
#' @param tree A guide tree specified in Newick format. Tip labels should be 0 to (nspec-1).
#' @param samplesizes A vector of sample sizes.
#' @param nspec The number of potential species.
#' @param popsizeprior Prior for population sizes.
#' @param divtimeprior Prior on divergence times.
#' @param migmatrix Matrix specifying which species migration can happen between.
#' @param nsnps The number of SNPs to simulate.
#' @param prefix A prefix for naming output filess.
#' @param migrateprior Prior on migration rates.
#' @param secondarycontact Boolean indicating whether or not to include secondary contact in models.
#' @param divwgeneflow Boolean indicating whether or not to include divergence with gene flow in models.
#' @param myrules A list of rules for the order of coalescence times.
#' @param maxmigrations The maximum number of migration events to include in any single model.

setup_fsc2 <- function(tree, samplesizes, nspec, popsizeprior,divtimeprior, migrateprior, migmatrix, nsnps, prefix,secondarycontact,divwgeneflow, maxmigrations, myrules = NULL){
  startnum = 0
  for(currenttree in 1:length(tree)){
    thistree <- tree[currenttree]
    mydir <- dir  # set a working directory to ensure all files are saved correctly
    checkinput(thistree, nspec, samplesizes, popsizeprior, divtimeprior) # check the user input for errors
    theparsed <- parsetree(thistree) # this returns a list of lists. The outer layer is a list of coalescent intervals. Within each interval, we list the coalescent events.
    if(  exists('theparsed') == FALSE) { # this will let the user know if there is an issue parsing the guide tree that wsn't caught by checkinput().
      stop("There was an issue parsing your guide tree")
    }
    myhistory <- makehistoryvector(theparsed) # get a simple list of coalescent events
    nhistevents <- length(myhistory) # count the coalescent events
    ntimes <- length(theparsed) # count the coalescent intervals
    theeventhistory <- gethistory(nhistevents,myhistory) # create a list of fsc2 style events for each coalescent event
    eventgrid <- geteventgrid(theeventhistory) # get a matrix of 0s and 1s for when you can and cannot collapse nodes
    mycollapsedhist <- historieswithcollapse(eventgrid, theeventhistory) # create fsc2 style events for all coalescent events with all possible collapses of nodes. This will be a list of lists. The outer list will be a list of models. Within each list we have a list of the relevant divergence events.
    migrationlist <- generatemigratelist(theeventhistory, theparsed, migmatrix) # get a list of migration events
    #print(migrationlist)
    myhistoricinfo <- generatemighistory(migrationlist, mycollapsedhist,secondarycontact,divwgeneflow,nspec,theparsed,maxnmigrate = maxmigrations)
    #print(myhistoricinfo)
    writetpl(samplesizes, nhistevents, nsnps, prefix, myhistoricinfo, startnum)
    endnum = writeest(samplesizes, nhistevents, nsnps, prefix, myhistoricinfo,popsizeprior,divtimeprior,migrateprior,myrules, startnum)
    startnum = endnum
  }
  if(length(tree) > 1){
    identical = checkidentical(prefix)  
  }

}
