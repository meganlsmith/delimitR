#' Generate .tpl files
#' @keywords internal

writetpl <- function(samplesizes, nhistevents, nsnps, prefix, myhistoricinfo,startnum){

  # we need to write a tpl file for each model
  for(model in 1:length(myhistoricinfo[[1]])){
#    cat('\nThis is model ',model,'\n')
    linecount <- 1
    mylines <- list()
    mylines[linecount] <- '//Number of population samples (demes)'
    linecount = linecount + 1
    mylines[linecount] <- as.character(length(samplesizes))
    linecount = linecount + 1
    mylines[linecount] <- '//Population effective sizes (number of genes)'
    linecount <- linecount + 1
    for(i in 1:length(samplesizes)){
      towrite <- paste('N',i-1,'$',sep="")
      mylines[linecount] <- towrite
      linecount <- linecount + 1
    }
    mylines[linecount] <- '//Sample sizes'
    linecount <- linecount + 1
    for(i in 1:length(samplesizes)){
      mylines[linecount] <- as.character(samplesizes[i])
      linecount <- linecount + 1
    }
    mylines[linecount] <- '//Growth rates: negative growth implies population expansion'
    linecount <- linecount + 1
    for(i in 1:length(samplesizes)){
      mylines[linecount] <- 0
      linecount <- linecount + 1
    }
    mylines[linecount] <- '//Number of migration matrices : 0 implies no migration between demes'
    linecount <- linecount + 1
    if(!is.matrix(myhistoricinfo[[2]][[model]])){
      mylines[linecount] <- 0
      linecount <- linecount + 1
    }
    else if(is.matrix(myhistoricinfo[[2]][[model]])){
      mylines[linecount] <- 2
      linecount <- linecount + 1
      mylines[linecount] <- '//migrationmatrix'
      linecount <- linecount + 1
      for (species in 1:length(samplesizes)){
        mylines[linecount] <- paste(as.character(rep(0,length(samplesizes))),sep = ' ',collapse = ' ')
        linecount <- linecount + 1

      }
      mylines[linecount] <- '//migrationmatrix'
      linecount <- linecount + 1
      matrix <- myhistoricinfo[[2]][[model]]
      for(row in 1:nrow(matrix)){
        towrite <- c()
        for(col in 1:ncol(matrix)){
          towrite <- c(towrite,matrix[row,col])
        }
        towrite = paste(towrite, sep = '',collapse = ' ')
        mylines[linecount] <- towrite
        linecount <- linecount + 1

      }
    }
    mylines[linecount] <- '//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix'
    linecount <- linecount + 1
    mylines[linecount] <- paste(length(myhistoricinfo[[1]][[model]]), ' historical event', sep = "")
    linecount <- linecount + 1
    for(historicevent in 1:length(myhistoricinfo[[1]][[model]])){
      mylines[linecount] <- myhistoricinfo[[1]][[model]][[historicevent]]
      linecount <- linecount + 1
    }
    mylines[linecount] <- '//Number of independent loci [chromosome]'
    linecount <- linecount + 1
    mylines[linecount] <- nsnps
    linecount <- linecount + 1
    mylines[linecount] <- '//Per chromosome: Number of linkage blocks'
    linecount <- linecount + 1
    mylines[linecount] <-1
    linecount <- linecount + 1
    mylines[linecount] <- '//per Block: data type, num loc. rec. rate and mut rate + optional paramters'
    linecount <- linecount + 1
    mylines[linecount] <- 'SNP 1 0 0 0.01'

    #print(mylines)
    filename <- paste(prefix,'_',model+startnum,'.tpl',sep = '')
    write(paste(mylines,sep='\n'), file=filename)


  }

}
