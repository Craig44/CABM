# Initialisating the agents with B0
# A script that I use to test different methods for
# setting the number of recruits during intialisation before
# we scale up to population level values


## The problem

# we have a specific B0 in mind that is represented by N-agents 
# how do we keep that constraint true? This gets very difficult
# if we have spatially varying M

# if we don't have spatial M we could get an M


### Mortality bias on how you add selectivity to an instantaneous rate
M = 0.2
Survival = 1 - exp(-M)
N = 100000;
age_specific_m = 0.33

chance = runif(N, 0,1)

sum(chance <= Survival) / N
Survival

sum(chance <= Survival * age_specific_m) / N
sum(chance <=  1 - exp(-M * age_specific_m)) / N


