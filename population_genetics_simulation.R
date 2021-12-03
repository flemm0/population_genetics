###
#creates initial population of 1,000 with 1 mutant
initialpop = c(1)
for(i in 1:999){
  initialpop = c(initialpop, 0)
}
newinput = initialpop  #assign name to initial population

nextGeneration = function(inputpop){ #function to return the next generation based on the input population
  biasedChancesVector = numeric(length(inputpop))
  for(i in 1:length(inputpop)){ ### creating vector containing biased chances of different selection rates for mutant vs wt
    if(inputpop[i] == 0){
      biasedChancesVector[i] = 0.97 #probability of wild-type reproducing
    }
    else{
      biasedChancesVector[i] = 1 #probability of mutant allele reproducing
    }
  }
  
  nextgen = sample(inputpop, size = 1000, replace = TRUE, prob = biasedChancesVector) #use sample function to generate the next population
  return(nextgen)
}

fixation = function(newinput){ #iterates nextGeneration until fixation
  mutationfreq = c()
  while(sum(newinput) != 0 & sum(newinput) != length(newinput)){ #while loop checks for fixation
    x = nextGeneration(newinput)
    newinput = x   
    mutationfreq = c(mutationfreq, (sum(newinput)/length(newinput))
  }
  
  return(mutationfreq) #returns the frequency of the mutant allele
}

probOfFixation=function(newinput){  #function to calculate the probability that the mutant allele will fix in 10,000 iterations
  mutantfixationscore = 0
  mutantextinctionscore = 0
  for(i in 1:10000){
    mutantfrequency = fixation(newinput)
    if(mutantfrequency[length(mutantfrequency)] == 1){
      mutantfixationscore = mutantfixationscore + 1
    }
    else{
      mutantextinctionscore = mutantextinctionscore + 1
    }
  }
  return(mutantfixationscore/10000)
  
}

fixSelection = function(newinput){    
  mutationfreq = c(sum(newinput)/1000)    ###initialize vector to store frequency of mutant in population of each generation 
  for(i in 1:100000){
    x = nextGeneration(newinput)
    newinput = x    ###uses temp variable x to allow the output of nextGeneration to become input of next iteration of nextGeneration
    mutationfreq = c(mutationfreq, (sum(newinput)/1000))     ###adds frequency of mutant to vector 
    if(sum(x) == 0 | sum(x) > 999){      ###if statement to check for fixation
      return(mutationfreq)
      break
    }
  }
}  
#fixSelection(newinput)


### below plots selection graph
plot(0, type='n',xlab="number of generations", ylab="frequency of mutant allele", ylim=c(0,1),xlim=c(1,500), main="Selection")
for(i in 1:10000){
  x = fixSelection(newinput)
  lines(x)
}

nextGenNoBias = function(newinput){ #simulates changes in population if the mutant and wild-type allele had the same fitness
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


###part 2
popsize = seq(from = 10, to = 100, by = 10) #create vector to store different sized populations
selectionstrength = seq(from = 0.005, to = 0.05, by = 0.005) #create vector to store different fitnesses of wild-type allele


selection = function(popsize, selectionstrength){ #create populations of individuals based on predetermined population sizes
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

###below plots average number of generations to reach fixation against selection strength for each population size in a selected mutation simulation
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


####below plots changes in population for when no selection bias for mutant or wild-type alleles
neutralVec = c()
yvals = c()
for(i in popsize){
  vec = avgGensToMutantFixation(i, 0)
  neutralVec = c(neutralVec, vec)
  yvals = c(yvals, 0)
}
plot(neutralVec, yvals, xlab='avg gens to mutant fixation', ylab='selection strength', main='Neutral Mutation', ylim=c(-0.05,0.05))