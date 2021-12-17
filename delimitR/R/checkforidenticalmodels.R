#' Generate input files for fastsimcoal2
#' @param prefix The prefix used for the analyses.

checkidentical <- function(prefix){
  
  mytplfiles <- list.files('./',pattern='*.tpl') # list of tpl files
  storeresults = list() # place to store hisotric events
  storemigration = list()
  
  # get historical events and file names
  for(file in 1:length(mytplfiles)){
    
    # read in this file
    currentfile = read.table(mytplfiles[file],sep="\n")
    
    # get the starting and ending points that you need for the file.
    number_of_events = currentfile[grep(" historical event", currentfile$V1),  ]
    starting_point = which(grepl(" historical event", currentfile$V1))
    number_of_events = as.numeric(strsplit(as.character(number_of_events),' ')[[1]][1])
    ending_point = starting_point + number_of_events
    starting_point = starting_point + 1
    mylines = currentfile[starting_point:ending_point,]
    mylines = as.character(mylines)
    mylines = c(mytplfiles[file],mylines)
    storeresults[[file]] <- mylines
    
    # get migration matrix if it exists
    migration_matrix = currentfile[grep("matrices", currentfile$V1),  ]
    startmigr = which(grepl(" matrices", currentfile$V1))
    getnumber = currentfile[startmigr:startmigr+1,]
    nummig = as.numeric(as.character(getnumber))
    if(nummig == 0){
      migmatrix = 0
      storemigration[[file]] <- migmatrix
      
    }
    else{
      startmatrix = which(grepl("igrationmatrix", currentfile$V1))[2]
      endmatrix = which(grepl("time, source", currentfile$V1))
      migmatrix= currentfile[(startmatrix+1):(endmatrix-1),]
      storemigration[[file]] <- migmatrix
      
    }
  }
  # check for identical models
  identicalmodels = list()
  for(model in 1:length(storeresults)){
    modellength = length(storeresults[[model]])-1
    modeldata = data.frame(matrix(data=NA, nrow=modellength, ncol=7))
    modelmigration = storemigration[[model]]
    names(modeldata) = c('time','source','sink','migrants','newsize','newgrowth','migmatrix')
    for(event in 2:length(storeresults[[model]])){
      thisevent = storeresults[[model]][[event]]
      thisevent = strsplit(thisevent,' ')
      modeldata[event-1,] <- thisevent[[1]]
    }
    # get rows where divergence happens at time zero, because unless there are at least two such rows, the issue cannot occur
    issuerows = modeldata[modeldata$time==0 & modeldata$migmatrix=='keep',]
    if(nrow(issuerows)>1){
      check=TRUE
    }else{
      check=FALSE
    }
    if(check==TRUE){ # if an issue is possible
      # get the same info for other models
      for(model2 in (model+1):length(storeresults)){
        model2length = length(storeresults[[model2]])-1
        model2migration = storemigration[[model2]]
        model2data = data.frame(matrix(data=NA, nrow=model2length, ncol=7))
        names(model2data) = c('time','source','sink','migrants','newsize','newgrowth','migmatrix')
        for(event2 in 2:length(storeresults[[model2]])){
          thisevent2 = storeresults[[model2]][[event2]]
          thisevent2 = strsplit(thisevent2,' ')
          model2data[event2-1,] <- thisevent2[[1]]
        }
        # again, pull out the rows with zero divergence times
        issuerows2 = model2data[model2data$time==0 & model2data$migmatrix=='keep',]
        if(nrow(issuerows2)>1){
          check2=TRUE
        }else{
          check2=FALSE
        }
        if(check2 == TRUE  & nrow(modeldata) == nrow(model2data)){ # if there are at least two of these, then check to see if the collapsed populations are the same, and if so add to list to remove. Also check to make sure the two have the same number of historic events, otherwise, one has migration and the other doesn't; they aren't identical
          same = c(as.numeric(issuerows$source), as.numeric(issuerows$sink))
          same2 = c(as.numeric(issuerows2$source), as.numeric(issuerows2$sink))
          same = sort(unique(same))
          same2 = sort(unique(same2))
          identical = identical(same, same2)
          if(identical == TRUE){
            migidentical = identical(modelmigration, model2migration)
            if(migidentical == TRUE){
              identicalmodels[[length(identicalmodels)+1]] <- c(storeresults[[model]][1], storeresults[[model2]][1])
              
            }

          }
         }
        
      }
    }
  }
  identicalmodels = unique(identicalmodels)
  removedmodels = c()
  # remove the identical models
  if(length(identicalmodels) > 0){
    for(identicalset in 1:length(identicalmodels)){
      if (identicalmodels[[identicalset]][2] %in% removedmodels){
        next()
      }
      else{
        command = paste('rm ',identicalmodels[[identicalset]][2], sep = '')
        system(command)
        command2 = gsub('tpl','est',command)
        system(command2)
        removedmodels = c(removedmodels, identicalmodels[[identicalset]][2])
        
      }
    }
    # check folder and rename models as appropriate
    
    mytplfiles <- list.files('./',pattern='*.tpl') # list of tpl files
    
    numbermodels = length(mytplfiles)
    
    mytplnums = sapply(mytplfiles, function(i) gsub(prefix,'',i))
    mytplnums = sapply(mytplnums, function(i) gsub('_','',i))
    mytplnums = sapply(mytplnums, function(i) gsub('.tpl','',i))
    mytplnums = sapply(mytplnums, function(i) as.numeric(i))
    currentsort = sort(mytplnums)
    desiredsort = 1:numbermodels
    for(modelname in 1:length(currentsort)){
      currentname = names(currentsort[modelname])
      desiredname = paste(prefix,'_',desiredsort[modelname], sep = '')
      commandmv = paste('mv ', currentname, ' ', desiredname,'_test.tpl', sep = '')
      commandmv2 = gsub('tpl','est',commandmv)
      system(commandmv)
      system(commandmv2)
      commandmv3 = paste('mv ', desiredname,'_test.tpl ', desiredname, '.tpl',sep = '')
      commandmv4 = gsub('tpl','est',commandmv3)
      system(commandmv3)
      system(commandmv4)
    }
    
  }
  
}