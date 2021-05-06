#########Number 1##############
initialpop = c(1)
for(i in 1:999){
  initialpop = c(initialpop, 0)
}
newinput = initialpop

nextGeneration = function(inputpop){
  biasedChancesVector = numeric(length(inputpop))
  for(i in 1:length(inputpop)){ ### creating vector containing biased chances
    if(inputpop[i] == 0){
      biasedChancesVector[i] = 0.97
    }
    else{
      biasedChancesVector[i] = 1
    }
  }
  
  nextgen = sample(inputpop, size = 1000, replace = TRUE, prob = biasedChancesVector)
  return(nextgen)
}

fixSelection = function(newinput){    
  mutationfreq = c(sum(newinput)/1000)    ###initialize vector to store frequency of mutant in population, with the initial mut freq recorded. 
  for(i in 1:100000){
    x = nextGeneration(newinput)
    newinput = x    ###uses temp variable x to allow the output of nextGeneration to become input of next iteration of nextGeneration
    mutationfreq = c(mutationfreq, (sum(newinput)/1000))     ###adds frequency of mutant to vector; should be down after the sampling process. 
    if(sum(x) == 0 | sum(x) > 999){      ###if statement to check for fixation
      return(mutationfreq)
      break
    }
  }
}  
#fixSelection(newinput)
###plots selection graph
plot(0, type='n',xlab="number of generations", ylab="frequency of mutant allele", ylim=c(0,1),xlim=c(1,500), main="Selection")
for(i in 1:10000){
  temp = fixSelection(newinput)
  lines(temp)
}

nextGenNoBias = function(newinput){
  mutantfreq = c(sum(newinput)/1000)
  for(i in 1:10000){
    x = sample(newinput, size = 1000, replace = TRUE)
    newinput = x
    mutantfreq = c(mutantfreq, (sum(newinput)/1000))
    if(sum(x) == 0 | sum(x) == 1000){
      return(mutantfreq)
      break
    }
  }
}

fixGenDrift = function(newinput){    
  mutantfreq = c(sum(newinput)/1000)
  while(sum(newinput) != 0 & sum(newinput) != 1000){
    x = nextGeneration(newinput)
    newinput = x
    mutantfreq = c(mutantfreq, (sum(newinput)/1000))
  }
  return(mutantfreq)
}

###plots genetic drift graph
plot(0, type='n',xlab="number of generations", ylab="frequency of mutant allele", ylim=c(0,1),xlim=c(1,500), main="Genetic Drift")
for(i in 1:10000){
  temp = fixSelection(newinput)
  lines(temp)
}


#########Number 2#############
popsize = seq(from = 10, to = 100, by = 10)
selectionstrength = seq(from = 0.005, to = 0.05, by = 0.005)


selection = function(popsize, selectionstrength){
  input = c(1)
  for(i in 1:(popsize-1)){
    input = c(input, 0)
  }
  
  mutantfreq = c() 
  for(i in 1:100000){
    biasedChancesVector = numeric(length(input))
    for(i in 1:length(input)){ 
      if(input[i] == 0){
        biasedChancesVector[i] = (1 - selectionstrength)
      }
      else{
        biasedChancesVector[i] = 1
      }
    }
    nextGen = sample(input, size = length(input), replace = TRUE, prob = biasedChancesVector)
    input = nextGen    ###uses temp variable x to allow the output of nextGeneration to become input of next iteration of nextGeneration
    mutantfreq = c(mutantfreq, sum(nextGen)/length(nextGen))    
    if(sum(nextGen) == 0 | sum(nextGen) == popsize){      ###if statement to check for fixation
      return(mutantfreq)
      break
    }
  }
}

avgGensToMutantFixation = function(popsize, selectionstrength){
  gens = c()
  for(x in 1:10000){
    vec = selection(popsize, selectionstrength)
    if(vec[length(vec)] == 1){
      gens = c(gens, length(vec))
    }
  }
  return(sum(gens)/length(gens))
}
### below plots average number of generations to reach fixation against selection strength for each population size in a selected mutation simulation


plot(0, type='n', xlab='avg gens to mutant fixation', ylab='selection strength', ylim = c(0,0.05), xlim=c(1,250), main="Selected Mutation")
vec = c()
for(i in popsize){
  xtemp = c()
  ytemp = c()
  for(j in selectionstrength){
    vec = avgGensToMutantFixation(i, j)
    xtemp = c(xtemp, vec)
    ytemp = c(ytemp, j)
  }
  lines(xtemp, ytemp)
}


####below is the plot for no selection bias
neutralVec = c()
yvals = c()
for(i in popsize){
    vec = avgGensToMutantFixation(i, 0)
    neutralVec = c(neutralVec, vec)
    yvals = c(yvals, 0)
}
plot(neutralVec, yvals, xlab='avg gens to mutant fixation', ylab='selection strength', main='Neutral Mutation', ylim=c(-0.05,0.05))
