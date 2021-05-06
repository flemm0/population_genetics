# population_genetics
Let's say there is a gene in the population that is fixed, meaning that there is only one allele present in the population.  Suddenly there is a new mutation so that there are now two alleles present.  
Moreover, let's say this new allele is actually beneficial:  an organism with this new mutant allele is slightly more likely to survive and reproduce than an organism with the more common (wild-type) allele.
In this project I simulate the evolution of the frequencies of these two alleles at this genetic locus to ultimately answer the question:  what is the chance that the new, mutant allele will take over and become fixed in the population?  
To do this, I set up a population of alleles, assigned them fitnesses, had them reproduce to make a new generation depending on their fitnesses, and iterate until one of two alleles is extinct and the other one is fixed.  Then I repeated this many times to determine what fraction of times the mutant is fixed and what fraction of the time it is the one that goes extinct.
I designated the wild-type allele by 0 and the mutant allele by 1. This is just a convention to keep them distinct. We will also assume that no more mutation happens.  

I plotted what happened to my 10,000 runs of genetic drift and selection to visualize what happens over time.
At the end, I estimated the number of generations needed for the mutation to go to fixation for a neutral and a selected mutation and varied population sizes and intensity of selection to answer the question: 
Is there a relationship between fixation, population size, and strength of selection?
