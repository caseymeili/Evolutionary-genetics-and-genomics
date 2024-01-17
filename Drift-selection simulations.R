# Evolutionary Genetics and Genomics (HGEN 6092)
# R Session 1 - Drift-selection simulations
# Casey Meili

# Some exercises were adapted from: https://evolutionarygenetics.github.io/

# Set population size
n <- 10
# Set allele frequency
p <- 0.5

# Sample a set of alleles for the next generation
gen1 <- rbinom(2*n,size=1,prob=p)
#print out the sampled alleles, where 1 is A1 and 0 is A2.
print(gen1)

# Calculate the allele frequency in this generation
p_gen1 <- sum(gen1)/(2*n)

# Print the allele frequency
p_gen1 

# Function to simulate genetic drift
drift_sim <- function(n, p, ngen){
  #p is the initial allele frequency
  #n is the sample size (# of diploid individuals)
  #ngen is the number of generations to simulate
  #ptime will keep a vector of the allele frequencies, adding the updated frequency after each  generation of sampling
  
  ptime <- p
  for(i in 1:ngen){
    newgen <- rbinom(2*n,size=1,prob=p) #sample a new set of 2n alleles
    p <- sum(newgen)/(2*n) #calculate the allele frequency in the next generation
    ptime<-append(ptime,p) #append the allele frequency for this generation to the previous generations
  }
  return(ptime)
}

# Run drift simulation for 50 generations
drift_sim(n = 10, p = 0.5, ngen = 100)

# Run 100 times with n=10, n=100 or n=1000
sim_n10 <- replicate(100,drift_sim(n=10,p=0.5,ngen=1000))
sim_n100 <- replicate(100,drift_sim(n=100,p=0.5,ngen=1000))
sim_n1000 <- replicate(100,drift_sim(n=1000,p=0.5,ngen=1000))

par(mfrow=c(3,1))
matplot(sim_n10,type="l",main="n=10",xlab="Generation",ylab="A1 frequency (p)")
matplot(sim_n100,type="l",main="n=100",xlab="Generation",ylab="A1 frequency (p)")
matplot(sim_n1000,type="l",main="n=1000",xlab="Generation",ylab="A1 frequency (p)")

# Show last row of matrix
table(sim_n10[1001,])
table(sim_n100[1001,])
table(sim_n1000[1001,])


