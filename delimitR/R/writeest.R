#' Generate .est files
#' @keywords internal

writeest <- function(samplesizes, nhistevents, nsnps, prefix, myhistoricinfo,popsizeprior,divtimeprior,migrateprior,myrules,startnum){
  countforend = 0
  for(model in 1:length(myhistoricinfo[[1]])){
    countforend = countforend + 1
    modelinfo = myhistoricinfo[[1]][[model]]
#    cat('\nThis is model ',model,'\n')
    linecount <- 1
    mylines <- list()
    mylines[linecount] <- '// Search ranges and rules file'
    linecount = linecount + 1
    mylines[linecount] <- '// ****************************'
    linecount = linecount + 1
    mylines[linecount] <- '[PARAMETERS]'
    linecount = linecount + 1
    mylines[linecount] <- '//#isInt? #name   #dist.#min  #max'
    linecount = linecount + 1
    mylines[linecount] <- '//all Ns are in number of haploid individuals'
    linecount = linecount + 1
    for(popsize in 1:length(samplesizes)){
      towrite <- paste('1   N',popsize-1,'$   unif   ',popsizeprior[[popsize]][1],'   ',popsizeprior[[popsize]][2],'   output', sep = '')
      mylines[linecount] <- towrite
      linecount = linecount + 1
    }
    divevents <- list()
    for(event in 1:length(modelinfo)){
      thisevent <- strsplit(modelinfo[[event]][[1]],' ')[[1]][1]
      if(grepl('div',thisevent)){
        mynum <- strsplit(thisevent,'v')[[1]][2]
        mynum <- gsub('$','',mynum,fixed = T)
        towrite <- paste('1   Tdiv',mynum,'$   unif   ',divtimeprior[[as.integer(mynum)]][1], '   ',divtimeprior[[as.integer(mynum)]][2], '   output',sep = '')
        mylines[linecount] <- towrite
        linecount = linecount + 1
        divevents <- c(divevents,paste('Tdiv',mynum,'$', sep = ''))
      }
    }
    migratepriorinfo <- myhistoricinfo[[2]][[model]]
    migratelist <- list(list)
    if(is.matrix(migratepriorinfo)){
      for(row in 1:nrow(migratepriorinfo)){
        for(col in 1:ncol(migratepriorinfo)){
          if(migratepriorinfo[row,col] != 0 & !(migratepriorinfo[row,col] %in% migratelist)){
            migstring <- strsplit(migratepriorinfo[row,col],'G')[[1]][2]
            migstring <- gsub('$','',migstring,fixed=T)
            towrite <- paste('0   ',migratepriorinfo[row,col],'   unif   ',migrateprior[[1]][1],'   ',migrateprior[[1]][2],'   output', sep = '')
            mylines[linecount] <- towrite
            linecount = linecount + 1
            migratelist <- c(migratelist, migratepriorinfo[row,col])
          }
        }
      }
    }
    mylines[linecount] <- ''
    linecount = linecount + 1
    mylines[linecount] <- '[RULES]'
    linecount = linecount + 1
    if(is.null(myrules)){
      mylines[linecount] <- ''
      linecount = linecount + 1
    }
    else if(!is.null(myrules)){
      for(rule in 1:length(myrules)){
        pop1 <- strsplit(x = myrules[[rule]],split = '$', fixed = T)[[1]][1]
        pop2 <- strsplit(myrules[[rule]],'T', fixed = T)[[1]][3]
        pop1 <- paste(pop1,'$',sep = '')
        pop2 <- paste('T',pop2,sep = '')
        if(pop1 %in% divevents & pop2 %in% divevents){

          mylines[linecount] <- myrules[[rule]]
          linecount = linecount + 1
        }
       }
      mylines[linecount] <- ''
      linecount = linecount + 1
    }
    mylines[linecount] <- '[COMPLEX PARAMETERS]'
    linecount = linecount + 1
    estcomplexinfo <- myhistoricinfo[[3]][[model]]
    if(!is.na(estcomplexinfo[1])){
      for(thingtowrite in 1:length(estcomplexinfo)){
        mylines[linecount]<- estcomplexinfo[[thingtowrite]]
        linecount = linecount + 1
      }
    }
    #print(mylines)
    filename <- paste(prefix,'_',model+startnum,'.est',sep = '')
    write(paste(mylines,sep='\n'), file=filename)
  }
  return(countforend + startnum)
}
