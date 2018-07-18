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
age_specific_m = 1

chance = runif(N, 0,1)

sum(chance <= Survival) / N
Survival

sum(chance <= Survival * age_specific_m) / N
sum(chance <=  1 - exp(-M * age_specific_m)) / N

start_year = 1990;
init_year = 50
first_year = start_year - init_year
length(first_year:(start_year - 1))

## initial distribution
M = 0.2
N = 10000
dist = rexp(N,0.2)
hist(dist, breaks = 30)

actual_dist = cut(dist, breaks = c(0:29,100))
age_dist = as.numeric(table(actual_dist))

seed = sum(dist < 1)

## track a cohort and create the same plot
nums = vector()
nums[1] = seed
for (i in 2:30)  {
  nums[i] = nums[i - 1] *  sum(chance <=  1 - exp(-M * age_specific_m)) / N
}

nums

## length and weight calculation
L_inf = 60
k = 0.2
t0 = 0

a = 





