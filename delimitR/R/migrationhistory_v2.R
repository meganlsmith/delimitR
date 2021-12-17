#' Generate history
#' @keywords internal
# This function will create a list of migration events. For each event, the list will have: population1, population2, tstart, and tstop
generatemigratelist <- function(history, parsedtree,migrationmatrix){
  migchecked <- list() # create a list to hold which migration events have already been added (so if we add 1<->2 we won't add 2<->1)
  migrationlist <- list(list(list())) # create a list of migration events that has population1, poulation2, NA, NA (we'll add start and stop times next)
  i = 1 # create a counter for adding to migration list
  for(row in 1:nrow(migrationmatrix)){ # loop through rows
    for(col in 1:ncol(migrationmatrix)){ # and columns of the migration matrix
      if((row != col & !(paste(col,row,sep='_') %in% migchecked)) & (migrationmatrix[row,col] == TRUE)){ # check to see if the row matches the column (avoid migration between 1<->1, etc), and if col, row is already in the migchecked list.
        addtocheck <- paste(row,col,sep='_') # if not, then add row, col to the migchecked list as a string
        migchecked <- c(migchecked,addtocheck) # add the string
        population1 <- row-1 # put the row population as pop1 (subtract 1 for 0 index)
        population2 <- col-1 # and the column population as pop2 (subtract 1 for 0 index)
        migrationlistadd <- list(population1,population2) # create vector
        migrationlist[[i]] <- migrationlistadd # add vector to list
        i<- i+1 # increment counter
      }
    }
  }
  return(migrationlist) # return the migration list
}







generatemighistory <- function(migrationlist, collapsedhistory,SC,DWG,nspec,theparsed, maxnmigrate){
  #This function will get our full history for our tpl file.

  i = 1 # set a counter to i for our historic events
  z = 1
  allhistoricevents <- list(list())# and create an empty list
  allmigmatrices <- list()
  estmig <- list(list())
  # loop through our models of collpased history
  for(model in 1:length(collapsedhistory)){
    allhistoricevents[[i]] <- collapsedhistory[[model]] # first, add the model without migration to the list
    i <- i+1 # and increment the counter
    allmigmatrices[[z]] <- NA
    estmig[[z]] <- NA
    z <- z+1
    currenthistory <- collapsedhistory[[model]] # now get the collapsed history

    # for each collapsed model, get the correct stop population for each migration event
    for(migration in 1:length(migrationlist)){
      thismigpop1 <- migrationlist[[migration]][[1]]
      thismigpop2 <- migrationlist[[migration]][[2]]
      FOUND <- FALSE
      for(divergenceevent in 1:length(currenthistory)){
        thisdivtime <- strsplit(currenthistory[[divergenceevent]],' ')[[1]][1]
        thisdivpop1 <- strsplit(currenthistory[[divergenceevent]],' ')[[1]][2]
        thisdivpop2 <- strsplit(currenthistory[[divergenceevent]],' ')[[1]][3]
        if((thisdivpop1 == thismigpop1 & thisdivpop2 == thismigpop2) | (thisdivpop1 == thismigpop2 & thisdivpop2 == thismigpop1)){
          if(thisdivtime != 0){
            migrationlist[[migration]][[3]] <- thisdivtime
            FOUND <- TRUE
          }
          else if(thisdivtime == 0){
            migrationlist[[migration]][[3]] <- NA
          }
        }
      } # end divergence for loop
      counter = 0
      while(FOUND == FALSE & counter <= length(currenthistory)){
        counter = counter + 1
        for(divergenceevent in 1:length(currenthistory)){
          thisdivtime <- strsplit(currenthistory[[divergenceevent]],' ')[[1]][1]
          thisdivpop1 <- strsplit(currenthistory[[divergenceevent]],' ')[[1]][2]
          thisdivpop2 <- strsplit(currenthistory[[divergenceevent]],' ')[[1]][3]
          if((thisdivpop1 == thismigpop1 & thisdivpop2 == thismigpop2) | (thisdivpop1 == thismigpop2 & thisdivpop2 == thismigpop1)){
            if(thisdivtime != 0){
              migrationlist[[migration]][[3]] <- thisdivtime
              FOUND <- TRUE
            }
            else if(thisdivtime == 0){
              migrationlist[[migration]][[3]] <- NA
              FOUND <- TRUE
            }
          }
          if(thisdivpop1 == thismigpop1){
            thismigpop1 = thisdivpop2
            break()
          }
          if(thisdivpop1 == thismigpop2){
            thismigpop2 = thisdivpop2
            break()
          }

        } # end divergence for loop

      }

    } # end migration for loop
    #print(migrationlist)


    for(scdwg in range(1:2)){ # we need to consider secondary contact & divergence with gene flow

      if(scdwg == 1 & SC == TRUE){ # if the user wants secondary contact, and we're on those models

        # we need to loop through the history and get a list of divergence events present, to help us decide which migration events to include
        tdivlist <- list() # create an empty list for storing divergence times
        nonzerotdivlist <- list () # create a list for storing non zero div times
        for(divevent in 1:length(currenthistory)){ # for each divergence event in the current history
          divtime <- strsplit(currenthistory[[divevent]],' ')[[1]][1] # get the name of the event
          tdivlist <- c(tdivlist,divtime) # and add it to our list
          if(divtime != 0){ # if the divergence time isn't zero, then add it to our non-zero list
            nonzerotdivlist <- c(nonzerotdivlist,divtime)
          }
        }
        j = 1 # initiate a counter for migrations to include
        migrationstoinclude <- list() # intiate a list

        # loop through potential migration movements, and see if they are appropriate for the current divergence model
        migpoplist <- c()
        for(potentialmigrationevent in 1:length(migrationlist)){
          PRESENT <- FALSE # an indicator for whether or not a migration is already included
          migration <- migrationlist[[potentialmigrationevent]] # get the event
          migpop1 <- migration[[1]] # get pop 1
          migpop2 <- migration[[2]] # get pop 2

          # if either of the involved populations merges into another at time zero, then we should change the population name to the one it merged with (unless they merge into teh same, in which case exclude)
          notinclude = FALSE # indicator for excluding
          for(divevent in 1:length(currenthistory)){
            divtime <- strsplit(currenthistory[[divevent]],' ')[[1]][1] # get div event name
            divpop1 <- strsplit(currenthistory[[divevent]],' ')[[1]][2] # get pop 1
            divpop2 <- strsplit(currenthistory[[divevent]],' ')[[1]][3] # get pop 2
            if(divtime==0 & ((divpop1 == migpop1 & divpop2 == migpop2) | (divpop1 == migpop2 & divpop2 == migpop1))){ # see if the divtime is zero, and if our two pops merged at time zero. Then, we should exclude the event
              notinclude <- TRUE # by changing the indicator flag to TRUE
            }

            if(divtime==0 & divpop1 == migpop1){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we change population
              migpop1 <- divpop2
            }
            if(divtime==0 & divpop1 == migpop2){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we change population
              migpop2 <- divpop2
            }
          }

          # now, we can look at our migraiton stop time, and, if it's in our divergence time list, we will add it to our migraitons to include
          migstop <- migration[[3]] # get stop times
          # loop through potential stop times
          for(potentialstop in 1:length(migstop)){
            thispotentialstop <- migstop[[potentialstop]] # get the potential stop time
            if(thispotentialstop %in% tdivlist & PRESENT == FALSE & notinclude == FALSE){ # if neither of our indicators are on, and the stop time is in the divergence time list
              checkpop1 = migpop1
              checkpop2 = migpop2
              mycheckstring = paste(checkpop1,'_',checkpop2, sep = '')
              mycheckstring2 = paste(checkpop2,'_',checkpop1, sep = '')
              if(mycheckstring %in% migpoplist | mycheckstring2 %in% migpoplist){
                PRESENT <- TRUE
              }
              else{
                migrationstoinclude[[j]] <- migration # add teh migraiton event
                j = j+1 # increment the counter

              }
              migpoplist <- c(migpoplist,mycheckstring)

            }
          }
        }
        if(length(migrationstoinclude) == 0){ # if we have no migration events to consider for the model
          next() # then move to the next model
        }

        # if we do have migraiton events to consider, then we need to consider all potential combinations of these events
        allcombos <- list() #create a list to hold these combinations
        k=1 # begin a counter for our list of combos
        for(potentiallength in 1:length(migrationstoinclude)){ # loop through 1 to the nubmer of migraiton to include and create combos with 1 to n events
          mycombos <- utils::combn(1:length(migrationstoinclude),m = potentiallength) # create combos
          allcombos[[k]] <- mycombos # add combos to combos list
          k=k+1 # increment the counter
        }

        # now we need to loop through these combos and create some models
        for(potentialcombonum in 1:length(allcombos)){
          if(potentialcombonum < (maxnmigrate+1)){
            thiscombonum <- allcombos[[potentialcombonum]] # store the current combo column
            for(potentialcombo in 1:ncol(thiscombonum)){ # look at each combination in the column, which will correspond to a model
              currentcombo <- thiscombonum[,potentialcombo] # get the combination
              currentmigrations <- list() #  make a list for the migration events this includes
              l = 1 # initiate a counter for this list
              for(migcombo in 1:length(currentcombo)){ # loop through the events included in teh current combo
                toget = currentcombo[migcombo] # get the event number
                currentmigrations[[l]]<- migrationstoinclude[toget] # and grab the migration list for the event & add it to our current migrations
                l = l+1
              }

              # now, we have a list of which migration events are included in the model. We need to make events and matrices for these. For SC scenarios, we will have only one matrix per model, and two events per model (one at zero, one halfway to the first coalescent)
              MYMATRIX <- matrix(rep(0,nspec*nspec),nrow=nspec,ncol=nspec,byrow=T)
              m = 1
              for(migrationitem in 1:length(currentmigrations)){
                pop1 <- currentmigrations[[migrationitem]][[1]][[1]]
                pop2 <- currentmigrations[[migrationitem]][[1]][[2]]
                for(divevent in 1:length(currenthistory)){
                  divtime <- strsplit(currenthistory[[divevent]],' ')[[1]][1] # get div event name
                  divpop1 <- as.integer(strsplit(currenthistory[[divevent]],' ')[[1]][2]) # get pop 1
                  divpop2 <- as.integer(strsplit(currenthistory[[divevent]],' ')[[1]][3]) # get pop 2

                  if(divtime==0 & divpop1 == pop1){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we should exclude the event
                    pop1 <- divpop2
                  }
                  if(divtime==0 & divpop1 == pop2){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we should exclude the event
                    pop2 <- divpop2
                  }
                }
                popadd1 <- pop1+1
                popadd2 <- pop2+1
                MYMATRIX[popadd1,popadd2] <- paste('MIG',m,'$',sep="")
                MYMATRIX[popadd2,popadd1] <- paste('MIG',m,'$',sep="")
                m = m+1
                #cat('pop 1 is ',pop1,', and pop 2 is ',pop2,'\n')

              }
              allmigmatrices[[z]] <- MYMATRIX

              # now that we have the matrix, we need the events
              startingevents <- collapsedhistory[[model]]
              toadd = length(startingevents) + 1
              startingevents[[toadd]] <- paste('0 0 0 0 1 0 1')
              toadd <- toadd+1
              startingevents[[toadd]] <- paste('Tmig1end$ 0 0 0 1 0 0')
              allhistoricevents[[i]] <- startingevents # first, add the model without migration to the list
              i <- i+1 # and increment the counter

              # now we need stop times for the est files
              #print(length(tdivlist))
              allestinfo <- list() # create an empty list to store info for est files

              # only consider divtimes that involve one of the species involved in migration
              divtoinclude <- nonzerotdivlist
              #print(divtoinclude)

              if(length(divtoinclude) == 1){ # if there's only 1 divergence time, then you just need that / 2
                #print('here')
                allestinfo <- paste('1 Tmig1end$ = ',divtoinclude[1],' / 2 output')
              }

              # we only want to consider multiple divergence times if they involve the species involved in migration.
              else if (length(divtoinclude) == 2){ # if there are two, take the minimum
                for(thisdivtime in 1:length(divtoinclude)){
                  estinfo <- paste('1 Tmig1end',thisdivtime,'$ = ',divtoinclude[[thisdivtime]],' / 2 output', sep = '')
                  allestinfo <- c(allestinfo, estinfo)
                }
                #print(estinfo)
                estinfo <- paste('1 Tmig1end$ = Tmig1end1$ %min% Tmig1end2$')
                allestinfo <- c(allestinfo, estinfo)
                #print(estinfo)
              }
              else{ # if there are more, take the minimum also.
                for(thisdivtime in 1:length(divtoinclude)){
                  estinfo <- paste('1 Tmig1end',thisdivtime,'$ = ',divtoinclude[[thisdivtime]],' / 2 output', sep = '')
                  #print(estinfo)
                  allestinfo <- c(allestinfo, estinfo)
                }
                for(thisdivtime in 1:length(divtoinclude)){
                  if(thisdivtime  == 1){
                    estinfo <- paste('1 Tmig1end',thisdivtime+length(divtoinclude),'$ = Tmig1end',thisdivtime,'$ %min% Tmig1end',thisdivtime+1,'$ output',sep = '')
                    #print(estinfo)
                    allestinfo <- c(allestinfo, estinfo)
                  }
                  else if (thisdivtime == 2){
                    next()
                  }
                  else {
                    estinfo <- paste('1 Tmig1end',thisdivtime+length(divtoinclude)-1,'$ = Tmig1end',thisdivtime+length(divtoinclude)-2,'$ %min% Tmig1end',thisdivtime,'$ output',sep = '')
                    allestinfo <- c(allestinfo, estinfo)

                  }
                }
                estinfo <- paste('1 Tmig1end$ = Tmig1end', thisdivtime+length(divtoinclude)-1, '$ * 1 output', sep = '')
                allestinfo <- c(allestinfo, estinfo)
              }
              #print(allestinfo)
              estmig[[z]] <- allestinfo
              #print(z)

              z <- z+1
            }
          }
        }
      } # end if for sc
      if(scdwg == 2 & DWG == TRUE){ # if the user wants dwg, and we're on those models
        currenthistory <- collapsedhistory[[model]] # now get the collapsed history

        # we need to loop through the history and get a list of divergence events present, to help us decide which migration events to include
        tdivlist <- list() # create an empty list for storing divergence times
        nonzerotdivlist <- list () # create a list for storing non zero div times
        for(divevent in 1:length(currenthistory)){ # for each divergence event in the current history
          divtime <- strsplit(currenthistory[[divevent]],' ')[[1]][1] # get the name of the event
          tdivlist <- c(tdivlist,divtime) # and add it to our list
          if(divtime != 0){
            nonzerotdivlist <- c(nonzerotdivlist,divtime)
          }
        }
        j = 1 # initiate a counter for migrations to include
        migrationstoinclude <- list() # intiate a list

        # loop through potential migration movements, and see if they are appropriate for the current divergence model
        migpoplist <- c()
        for(potentialmigrationevent in 1:length(migrationlist)){
          PRESENT <- FALSE # an indicator for whether or not a migration is already included
          migration <- migrationlist[[potentialmigrationevent]] # get the event
          #print(migration[[4]])
          migpop1 <- migration[[1]] # get pop 1
          migpop2 <- migration[[2]] # get pop 2
          additionalstoptimes <- c()
          # now, if our two populations of interest have non-zero div times with each other, we will include them, provided there isn't a coalescent event for one of them in a previous interval
          for(divevent in 1:length(currenthistory)){
            divtime <- strsplit(currenthistory[[divevent]],' ')[[1]][1] # get div event name
            divpop1 <- strsplit(currenthistory[[divevent]],' ')[[1]][2] # get pop 1
            divpop2 <- strsplit(currenthistory[[divevent]],' ')[[1]][3] # get pop 2
            if(divtime==0 & ((divpop1 == migpop1 & divpop2 == migpop2) | (divpop1 == migpop2 & divpop2 == migpop1))){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we should exclude the event
              notinclude <- TRUE # by changing the indicator flag to TRUE
            }
            else if(divtime==0 & divpop1 == migpop1){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we should exclude the event
              migpop1 = divpop2
            }
            else if(divtime==0 & divpop1 == migpop2){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we should exclude the event
              migpop2 = divpop2
            }
            if(divtime!=0 & ((divpop1 == migpop1 & divpop2 == migpop2) | (divpop1 == migpop2 & divpop2 == migpop1)) & PRESENT == FALSE){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we should exclude the event
              currentstops <- migration[[3]]
              #cat('adding ', divtime, '\n')
              migration[[3]] <- c(currentstops,divtime)
              # what is the coalescent interval
              include <- TRUE
              for(coalinterval in 1:length(theparsed)){
                #print(length(theparsed[[coalinterval]]))
                for(coalevent in 1:length(theparsed[[coalinterval]])){
                  #print(theparsed[[coalinterval]][coalevent])
                  coalpop1 <- theparsed[[coalinterval]][[coalevent]][1]
                  coalpop2 <- theparsed[[coalinterval]][[coalevent]][3]
                  if((coalpop1 == divpop1 & coalpop2 == divpop2) | (coalpop1 == divpop2 & coalpop2 == divpop1)){
                    mycoalinterval = coalinterval
                  }
                }
              }
              for(divevent2 in 1:length(currenthistory)){
                divtime2 <- strsplit(currenthistory[[divevent2]],' ')[[1]][1] # get div event name
                div2pop1 <- strsplit(currenthistory[[divevent2]],' ')[[1]][2] # get pop 1
                div2pop2 <- strsplit(currenthistory[[divevent2]],' ')[[1]][3] # get pop 2
                if((divpop1 == div2pop1 | divpop1 == div2pop2 | divpop2 == div2pop1 | divpop2 == div2pop2) & divtime2 != 0){
                  # what is the coalescent interval

                  for(coalinterval in 1:length(theparsed)){
                    #print(length(theparsed[[coalinterval]]))
                    for(coalevent in 1:length(theparsed[[coalinterval]])){
                      #print(theparsed[[coalinterval]][coalevent])
                      coalpop1 <- theparsed[[coalinterval]][[coalevent]][1]
                      coalpop2 <- theparsed[[coalinterval]][[coalevent]][3]
                      if((coalpop1 == div2pop1 & coalpop2 == div2pop2) | (coalpop1 == div2pop2 & coalpop2 == div2pop1)){
                        thiscoalinterval = coalinterval
                        if(thiscoalinterval < mycoalinterval){
                          include <- FALSE
                          next()
                        }
                      }
                    }
                  }


                }

              }

              if (include == TRUE){
                checkpop1 = migpop1
                checkpop2 = migpop2
                mycheckstring = paste(checkpop1,'_',checkpop2, sep = '')
                mycheckstring2 = paste(checkpop2,'_',checkpop1, sep = '')
                if(mycheckstring %in% migpoplist | mycheckstring2 %in% migpoplist){
                  PRESENT <- TRUE
                }
                else{
                  migrationstoinclude[[j]] <- migration # add teh migraiton event
                  #                  print('adding this')
                  migpoplist <- c(migpoplist,mycheckstring)
                  j = j+1 # increment the counter
                  PRESENT <- TRUE # and flip the indicator
                }
              }
            }
          }

          # loop through potential stop times
        }
        #print(migrationstoinclude)
        if(length(migrationstoinclude) == 0){ # if we have no migration events to consider for the model
          next() # then move to the next model
        }

        # if we do have migraiton events to consider, then we need to consider all potential combinations of these events
        allcombos <- list() #create a list to hold these combinations
        k=1 # begin a counter for our list of combos
        for(potentiallength in 1:length(migrationstoinclude)){ # loop through 1 to the nubmer of migraiton to include and create combos with 1 to n events
          mycombos <- utils::combn(1:length(migrationstoinclude),m = potentiallength) # create combos
          allcombos[[k]] <- mycombos # add combos to combos list
          k=k+1 # increment the counter
        }

        # now we need to loop through these combos and create some models
        for(potentialcombonum in 1:length(allcombos)){
          if(potentialcombonum < maxnmigrate+1){
            thiscombonum <- allcombos[[potentialcombonum]] # store the current combo column

            for(potentialcombo in 1:ncol(thiscombonum)){ # look at each combination in the column, which will correspond to a model
              #print('new combo')
              currentcombo <- thiscombonum[,potentialcombo] # get the combination
              currentmigrations <- list() #  make a list for the migration events this includes
              l = 1 # initiate a counter for this list
              for(migcombo in 1:length(currentcombo)){ # loop through the events included in teh current combo
                toget = currentcombo[migcombo] # get the event number
                currentmigrations[[l]]<- migrationstoinclude[toget] # and grab the migration list for the event & add it to our current migrations
                l = l+1
              }
              # now, we have a list of which migration events are included in the model. We need to make events and matrices for these. For SC scenarios, we will have only one matrix per model, and two events per model (one at zero, one halfway to the first coalescent)
              MYMATRIX <- matrix(rep(0,nspec*nspec),nrow=nspec,ncol=nspec,byrow=T)
              m = 1
              for(migrationitem in 1:length(currentmigrations)){
                pop1 <- currentmigrations[[migrationitem]][[1]][[1]]
                pop2 <- currentmigrations[[migrationitem]][[1]][[2]]
                for(divevent in 1:length(currenthistory)){
                  divtime <- strsplit(currenthistory[[divevent]],' ')[[1]][1] # get div event name
                  divpop1 <- as.integer(strsplit(currenthistory[[divevent]],' ')[[1]][2]) # get pop 1
                  divpop2 <- as.integer(strsplit(currenthistory[[divevent]],' ')[[1]][3]) # get pop 2

                  if(divtime==0 & divpop1 == pop1){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we should exclude the event
                    pop1 <- divpop2
                  }
                  if(divtime==0 & divpop1 == pop2){ # see if the divtime is zero, and divpop 1 is either of our migration pops. Then, we should exclude the event
                    pop2 <- divpop2
                  }
                }
                popadd1 <- pop1+1
                popadd2 <- pop2+1
                MYMATRIX[popadd1,popadd2] <- paste('MIG',m,'$',sep="")
                MYMATRIX[popadd2,popadd1] <- paste('MIG',m,'$',sep="")
                m = m+1
                #cat('pop 1 is ',pop1,', and pop 2 is ',pop2,'\n')

              }
              allmigmatrices[[z]] <- MYMATRIX

              # now that we have the matrix, we need the events
              startingevents <- collapsedhistory[[model]]
              toadd = length(startingevents) + 1
              startingevents[[toadd]] <- paste('Tmig1init$ 0 0 0 1 0 1')
              toadd <- toadd+1
              startingevents[[toadd]] <- paste('Tmig1end$ 0 0 0 1 0 0')
              allhistoricevents[[i]] <- startingevents # first, add the model without migration to the list
              i <- i+1 # and increment the counter

              # now we need stop times for the est files
              #print(length(tdivlist))
              allestinfo <- list() # create an empty list to store info for est files

              # okay, now for the est files... these are more complicated in the divw gene flow case (maybe)...
              #print(currentmigrations)
              if(length(currentmigrations) == 1){
                stoptime = currentmigrations[[1]][[1]][3][[1]]
                newstoptime = list()
                for(potentialstoptime in 1:length(stoptime)){
                  if(stoptime[[potentialstoptime]] %in% nonzerotdivlist & !(stoptime[[potentialstoptime]] %in% newstoptime)){
                    newstoptime = c(newstoptime,stoptime[[potentialstoptime]])
                  }
                }
                stoptime = newstoptime
                #cat('mig pop 1 is ',migpop1,' and mig pop 2 is ', migpop2, '\n')
                if(length(stoptime) > 1){
                  for(thiscurrentstoptime in 1:length(stoptime)){
                    estinfo <- paste('1 Tmig1init',thiscurrentstoptime,'$ = ',stoptime[[thiscurrentstoptime]],' / 2 output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)
                    estinfo <- paste('1 Tmig1end',thiscurrentstoptime,'$ = ',stoptime[[thiscurrentstoptime]],' * 1 output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)

                  }
                  if(length(stoptime) == 2){
                    estinfo <- paste('1 Tmig1init','$ = ','Tmig1init1$ %min% Tmig1init2$ output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)
                    estinfo <- paste('1 Tmig1end$ = Tmig1end1$ %min% Tmig1end2$ output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)

                  }
                  if(length(stoptime)>2){
                    for(time in 1:length(stoptime)){
                      if (time == 1){
                        estinfo <- paste('1 Tmig1init',length(stoptime)+time,'$ = ','Tmig1init1$ %min% Tmig1init2$ output',sep = "")
                        allestinfo <- c(allestinfo,estinfo)
                        estinfo <- paste('1 Tmig1end',length(stoptime)+time,'$ = ','Tmig1end1$ %min% Tmig1end2$ output',sep = "")
                        allestinfo <- c(allestinfo,estinfo)
                      }
                      if(time == 2){
                        next()
                      }
                      if(time > 2){
                        estinfo <- paste('1 Tmig1init', length(stoptime)+time-1 ,'$ = Tmig1init',length(stoptime)+time-2,'$ %min% Tmig1init',time,'$ output',sep = "")
                        allestinfo <- c(allestinfo,estinfo)
                        estinfo <- paste('1 Tmig1end', length(stoptime)+time-1 ,'$ = Tmig1end',length(stoptime)+time-2,'$ %min% Tmig1end',time,'$ output',sep = "")
                        allestinfo <- c(allestinfo,estinfo)

                      }
                    }
                    estinfo <- paste('1 Tmig1init$ = Tmig1init', length(stoptime)+time - 1 ,'$ * 1 output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)
                    estinfo <- paste('1 Tmig1end$ = Tmig1end', length(stoptime)+time - 1 ,'$ * 1 output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)

                  }


                }
                else{
                  estinfo <- paste('1 Tmig1init$ = ',stoptime,' / 2 output',sep = "")
                  allestinfo <- c(allestinfo,estinfo)
                  estinfo <- paste('1 Tmig1end$ = ',stoptime,' * 1 output',sep = "")
                  allestinfo <- c(allestinfo,estinfo)
                }
              }
              else if(length(currentmigrations) == 2){
                allstops <- list()
                for(themigration in 1:length(currentmigrations)){
                  stoptime = currentmigrations[[themigration]][[1]][3]
                  newstoptime = list()
                  for(potentialstoptime in 1:length(stoptime[[1]])){
                    if(stoptime[[1]][[potentialstoptime]] %in% nonzerotdivlist & !(stoptime[[1]][[potentialstoptime]] %in% newstoptime)){
                      newstoptime = c(newstoptime,stoptime[[1]][[potentialstoptime]])
                    }
                  }
                  stoptime = newstoptime
                  allstops <- c(allstops, stoptime)
                }
                for(timetostop in 1:length(allstops)){
                  stoptime = allstops[[timetostop]]
                  estinfo <- paste('1 Tmig1init',timetostop,'$ = ',stoptime,' / 2 output',sep = "")
                  allestinfo <- c(allestinfo,estinfo)
                  estinfo <- paste('1 Tmig1end',timetostop,'$ = ',stoptime,' * 1 output',sep = "")
                  allestinfo <- c(allestinfo,estinfo)

                }
                if(length(allstops) == 2){
                  estinfo <- paste('1 Tmig1init$ = Tmig1init1$ %min% Tmig1init2$ output',sep = "")
                  allestinfo <- c(allestinfo,estinfo)
                  estinfo <- paste('1 Tmig1end$ = Tmig1end1$ %min% Tmig1end2$ output',sep = "")
                  allestinfo <- c(allestinfo,estinfo)

                }
                else if (length(allstops) > 2){
                  for(timetostop in 1:length(allstops)){
                    if(timetostop == 1){
                      estinfo <- paste('1 Tmig1init', length(allstops)+1 ,'$ = Tmig1init1$ %min% Tmig1init2$ output',sep = "")
                      allestinfo <- c(allestinfo,estinfo)
                      estinfo <- paste('1 Tmig1end', length(allstops)+1 ,'$ = Tmig1end1$ %min% Tmig1end2$ output',sep = "")
                      allestinfo <- c(allestinfo,estinfo)
                    }
                    else if (timetostop == 2){
                      next()
                    }
                    else{
                      estinfo <- paste('1 Tmig1init', length(allstops)+timetostop-1 ,'$ = Tmig1init',length(allstops)+timetostop-2,'$ %min% Tmig1init',timetostop,'$ output',sep = "")
                      allestinfo <- c(allestinfo,estinfo)
                      estinfo <- paste('1 Tmig1end', length(allstops)+timetostop-1 ,'$ = Tmig1end',length(allstops)+timetostop-2,'$ %min% Tmig1end',timetostop,'$ output',sep = "")
                      allestinfo <- c(allestinfo,estinfo)

                    }
                  }
                }
              }
              else if(length(currentmigrations) > 2){
                for(themigration in 1:length(currentmigrations)){
                  stoptime = currentmigrations[[themigration]][[1]][3]
                  newstoptime = list()
                  for(potentialstoptime in 1:length(stoptime[[1]])){
                    if(stoptime[[1]][[potentialstoptime]] %in% nonzerotdivlist & !(stoptime[[1]][[potentialstoptime]] %in% newstoptime)){
                      newstoptime = c(newstoptime,stoptime[[1]][[potentialstoptime]])
                    }
                  }
                  stoptime = newstoptime
                  estinfo <- paste('1 Tmig1init',themigration,'$ = ',stoptime,' / 2 output',sep = "")
                  allestinfo <- c(allestinfo,estinfo)
                  estinfo <- paste('1 Tmig1end',themigration,'$ = ',stoptime,' * 1 output',sep = "")
                  allestinfo <- c(allestinfo,estinfo)
                }
                for(themigration in 1:length(currentmigrations)){
                  if(themigration == 1){
                    estinfo <- paste('1 Tmig1init',themigration + length(currentmigrations),'$ = Tmig1init1$ %min% Tmig1init2$ output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)
                    estinfo <- paste('1 Tmig1end',themigration + length(currentmigrations),'$ = Tmig1end1$ %min% Tmig1end2$ output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)
                  }
                  if(themigration == 2) {
                    next()
                  }
                  if(themigration > 2) {
                    estinfo <- paste('1 Tmig1init',themigration + length(currentmigrations) - 1 ,'$ = Tmig1init',themigration + length(currentmigrations) - 2,'$ %min% Tmig1init',themigration,'$ output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)
                    estinfo <- paste('1 Tmig1end',themigration + length(currentmigrations) - 1 ,'$ = Tmig1end',themigration + length(currentmigrations) - 2,'$ %min% Tmig1end',themigration,'$ output',sep = "")
                    allestinfo <- c(allestinfo,estinfo)

                  }
                }
                estinfo <- paste('1 Tmig1init$ = Tmig1init',length(currentmigrations) + length(currentmigrations) - 1 ,'$ * 1 output', sep = "")
                allestinfo <- c(allestinfo,estinfo)
                estinfo <- paste('1 Tmig1end$ = Tmig1end',length(currentmigrations) + length(currentmigrations) - 1 ,'$ * 1 output', sep = "")
                allestinfo <- c(allestinfo,estinfo)
              }
              estmig[[z]] <- allestinfo
              #print(z)
              z <- z+1

            }
          }
        }
      } # end loop for dgf
    } #end loop for sc vs dgf
  }
  #print(estmig)
  return(list(allhistoricevents, allmigmatrices,estmig))
}
