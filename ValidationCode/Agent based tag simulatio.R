# Routine to test agent based tagging process
# June 13th 2016 J McKenzie
# @author J McKenzie

#############
# FUNCTIONS #
#############               

#RMSD

RMSD<- function(x,y){
  
  
  CV<-round(sqrt( sum((x-y)^2)/length(y))/x,2)
  bias<-round((mean(y)/x)-1,2)
  return(cbind('true'=x,'est_mean'=round(mean(y),2),'bias'=bias,'RMCV'=CV,'obs'= length(y)))
}

# Resequence
reseq<- function(x){
  seq_ref<- 1
  for(i in 1:length(x[,1])){
    x[i,3]<-seq_ref
    seq_ref<-seq_ref+x[i,2]
    x[i,4]<-seq_ref-1
  }
  x[,1]<-c(1:length(x[,1]))
  return(x)
}

# point sequence
## What does this function do?
#
# x = ndx relating to individuals (1-pop) who are randomly tagged from the entire population
# y = data frame that maps agents to individuals from original population
# assign tags to agents based on individuals (agents can consist of multiple tagged individuals)
pointer<- function(x,y){
  point<-NULL
  for(i in x){
    point<-c(point,grep(1,(y[,3] <= i) * (y[,4] >= i)) )
  }
  point<-point[order(point)]
  point<-aggregate(point,list(point),length)
  colnames(point)<-c('point','freq') 
  return(point)
}







################
# Main Routine #
################

# Population size
pop = 1000

# maximum agent scaler value
mscalr = 1


# Number of tags
tags = 50

# Number of catch (abundance)
catch = 50 
## We won't have this information we have biomass, but we won't know numbers because this is stochastic and
## because there is a selectivity involved it, makes it tricky

boot<- 50
#######################

RunDet<-list()


RunDet[c('pop','mscalr','tags','catch','boot','runtime','stats')]<-c(pop,mscalr,tags,catch,boot,0,0)

pop_est<-NULL
cyrils_pop_est<-NULL

startTime<-Sys.time()

for(bt in 1:boot){  
  print(bt)
  # trap silly mscalar values
  if(mscalr >= pop ) mscalr <- pop-1
  if(mscalr <= 0 ) mscalr <- 1
  
  #### Create agent matrix
  agent_vec<- rep(mscalr,as.integer(pop/mscalr))
  if(pop %% mscalr > 0) 
  	agent_vec <- c(agent_vec,pop %% mscalr)
  
  agent_mat<-as.data.frame(cbind('agent'=1:length(agent_vec), 'number'=agent_vec, 'S_seq'=(c(0,cumsum(agent_vec)[1:(length(agent_vec)-1)]))+1,
                             'E_seq'=((c(0,cumsum(agent_vec)[1:(length(agent_vec)-1)]))+1)+(agent_vec-1), 'tag'=0))


  ##########################################
  
  # Apply Tagging
  ## agent ndx with tags
  tagseq<-sample(c(1:pop), tags, replace = FALSE, prob = NULL)
  point<-pointer(tagseq,agent_mat)
  
  ## pull out true tagged individuals from the agent dataframe, and rejig the sequence
  for(i in 1:length(point[,1])){
      for(j in 1:point[i,2]){
        if(agent_mat[point[i,1],2]==1){
          agent_mat[point[i,1],5]= 1
        }else{
            zz<-agent_mat[point[i,1],]
            zz[,3]<-zz[,4]
            zz[,5]<-1
            zz[,2]<-1
            agent_mat[point[i,1],4]<- agent_mat[point[i,1],4]-1
            agent_mat[point[i,1],2]<- agent_mat[point[i,1],2]-1
            agent_mat<-rbind(agent_mat,zz)
        }
  
      }
      
  }
  
  agent_mat<-agent_mat[order(agent_mat[,4]),]
  ## reset agent ndx to be new agent ndx 
  agent_mat[,1]<-seq(1,length(agent_mat[,1]),1)
  
  ## check everything was applied correctly
  pop<- sum(agent_mat[,2])
  tagsum<- sum(agent_mat[,5]) 
  #############################
  # Take catch routines 
  #############################
  # catch tag as individual  catch untagged as agent
  # recover tags
  catseq<-sample(c(1:pop), catch, replace = FALSE, prob = NULL)
  point_tag<-pointer(catseq,agent_mat)
  
  catch_tag<-NULL
  del_vec<-NULL
  for(i in 1:length(point_tag[,1])){
      if(agent_mat[point_tag[i,1],5]==1){
        catch_tag<-rbind(catch_tag,agent_mat[point_tag[i,1],])
        del_vec<-c(del_vec,-1*point_tag[i,1])
      }
  }
  norectags<- sum(catch_tag[,5])
  ## remove tags from original partition
  if(!is.null(del_vec)) 
    agent_mat<-agent_mat[del_vec,]

  #print( catch_tag, row.names = F)
  
  # split out uncaught tagged fish
  Tags<- agent_mat[grep(1,agent_mat[,5]),]
  untagged<-agent_mat[grep(0,agent_mat[,5]),]
  untagged<-reseq(untagged)

  
  ############
  # untagged catch removal
  cat_untag<- catch - norectags
  while(cat_untag>0) {
    catseq<-sample(sum(untagged[,2]), 1, replace = FALSE, prob = NULL)
    point <- pointer(catseq,untagged)
    cat <- untagged[point[1,1],2]
    if(cat_untag >= cat)
    {
      cat_untag<-cat_untag-cat
      untagged<-untagged[-point[1,1],]
      untagged<-reseq(untagged)
    } else{

      untagged[point[1,1],2]<-untagged[point[1,1],2]-cat_untag
      untagged[point[1,1],4]<-untagged[point[1,1],4]-cat_untag
      untagged <- reseq(untagged)
      cat_untag <- 0
    }  
  }
  #####################
  ## Cyrils alternative
  #####################
  # The Algorthm
  # - probability of sampling an Agent, depends on how many individuals an Agent represents.
  # - So it's a weighted resampling algorithm
  # - This a little hard, in the ABM C++ mechanism, firstly because multiple agents will have 
  # different amount of individuals, if there are different 'stocks', so what we do is calculate
  # the average, 
  agent_mat$weight = agent_mat$number / sum(agent_mat$number)
  cat_untag<- catch
  agents_sampled = c();
  tag_count = 0;
  max_attempts = nrow(agent_mat) * 10
  attempt = 1;
  agents_to_sample = sample(size = nrow(agent_mat), x = 1:nrow(agent_mat), prob = agent_mat$weight)
  for(i in 1:length(agents_to_sample)) {
    ## check if agent vulnerable or already been sampled
    cat_untag = cat_untag - agent_mat[agents_to_sample[i],"number"]
    if(agent_mat[agents_to_sample[i],"tag"] == 1)
      tag_count = tag_count + 1;
    if (cat_untag <= 0)
      break;
  }

  
  agent_mat<-reseq(rbind(Tags,untagged))
  # new_pop<-sum(agent_mat[,2])
  #stats
  tagrec<-sum(catch_tag[,5])
  pop_est<- rbind(pop_est,cbind('pop_est'= ((tags+1)*(catch+1)/(tagrec+1))-1, 'tagsrec'=tagrec)) # Petersen estimator with Chapman bias correction
  cyrils_pop_est <- rbind(cyrils_pop_est, cbind('pop_est'= ((tags+1)*(catch+1)/(tag_count+1))-1, 'tagsrec'=tag_count))

} #endboot

endTime<-Sys.time()
runTime<-endTime-startTime

RMSD(pop,pop_est[,1])
RMSD(pop,cyrils_pop_est[,1])
# stats


RunDet[['runtime']]<-runTime
RunDet[['stats']]<-RMSD(pop,pop_est[,1])


RunDet





## issues:
## we deal with biomass not abundance, so can't do these ahead calculations as is laid out in this code.
## don't know vulnerable population if selectivity is involved as this is a stochastic process
## 



# I think we generate a sample() like in R, instead of current uniform selection the probability
# of being selected for catch is less than an agent. Problem this will be dynamic, due to sample with 
# out replacement
# an experiment about generating random numbers
set.seed(123);
weights = floor(rlnorm(2000,log(5), 0.4))
vals = 1:max(weights)
Freq.y=tabulate(weights,nbins=max(weights))/length(weights)
## we can attach a weight to each agent before executing the calculation, which when summed over all 
## agents would be = 1
prob = weights / sum(weights)
## resample ndx 
N_rep = 50
N_sample = 1000
R_sample = my_sample = alt_sample = matrix(NA,nrow = N_rep, ncol = N_sample);

for(i in 1:N_rep) {
	full_draw = sample(1:length(weights), replace = F, size = length(weights), prob = prob)
	my_sample[i,] = weights[full_draw[1:N_sample]]
	R_sample[i,] = sample(weights, replace = F, size = N_sample, prob = prob)
	alt_sample[i,] = sample(x = vals, size = N_sample, prob = Freq.y, replace = T)
}

Freq.R=tabulate(as.vector(R_sample),nbins=max(weights))/(N_rep* N_sample)
Freq.my=tabulate(as.vector(my_sample),nbins=max(weights))/(N_rep* N_sample)
Freq.alt=tabulate(as.vector(alt_sample),nbins=max(weights))/(N_rep* N_sample)

barplot(Freq.y,main="", names = vals)
for(i in 1:N_rep) {
  check = tabulate(my_sample[i,], nbins =max(weights))
  check = check / sum(check)
}

barplot(rbind(Freq.y,Freq.alt,Freq.R, Freq.my),beside=T,main="")


## Also find out at what point should we even care about this. when does this really become an issue
## vulnerable in a fishing event
boots = 1000;
set.seed(123)
tags_to_look_at = c(100, 1000, 2000, 5000)
scalars = c(100,50,10,5,1)
tag_counts = results = array(0, dim = c(boots, length(scalars), length(tags_to_look_at)))
population = 1000000
## should possible look at these two
catch = 5000 # abundance - individuals, scanning 

startTime<-Sys.time()
scalar_ndx = 1;
tag_ndx = 1;
for(scalar in scalars) {
	print(scalar)
  tag_ndx = 1;
	for(n_tags in tags_to_look_at) {
		print(paste0("tag = ", n_tags))
	  n_agents = as.integer(population / scalar)
		tags = rep(0,n_agents)
		weights = rep(scalar, n_agents)
		# apply tagging, all agents have equal prob of being tagged
		tag_agent_ndx = sample(x = n_agents, size = n_tags, replace = F)
		## pull out tagged individuals take away the scalar
		weights[tag_agent_ndx] = weights[tag_agent_ndx] - 1
		untagged_pop = data.frame(weight = weights, tagged = 0)
		tagged_pop = data.frame(weight = rep(1,n_tags), tagged = 1)
		
		## now do recapture event first treating this is as equal probability
	  prob_untagged = sum(untagged_pop$weight) / population
	  prob_tagged = sum(tagged_pop$weight) / population
	  #full_prob = c(prob_tagged, prob_untagged)
	  ## Adapt probability for sampling of entire agents at a time
	  prob_untagged = prob_untagged / scalar
	  prob_tagged = prob_tagged / (prob_untagged + prob_tagged)
	  prob_untagged = 1 - prob_tagged
	  pop_est = NULL
	  for(bt in 1:boots) {
	    temp_untagged_pop = untagged_pop
	    temp_tagged_pop = tagged_pop
	    temp_catch = catch
	    ## now do recapture event first treating this is as equal probability
	    tag_count = 0
	    agent_ndx_sampled = c()
	    while(temp_catch > 0) {
	      agent_tag = rbinom(n = c(1), size = 1, prob = prob_tagged)
	      ## now find agent and sample
	      if (agent_tag == 0) {
	        agent_ndx = sample(1:nrow(temp_untagged_pop), size = 1)
	        counter = 0;
	        while ((temp_untagged_pop$weight[agent_ndx] == 0) | (counter > 5)) {
	          agent_ndx = sample(1:nrow(temp_untagged_pop), size = 1)
	          counter = counter + 1
	        }
          if (counter > 5)
            break
	        agent_ndx_sampled = c(agent_ndx_sampled,agent_ndx)
	        temp_catch = temp_catch - temp_untagged_pop$weight[agent_ndx];
	        temp_untagged_pop$weight[agent_ndx] = 0.0
        } else  {
	        agent_ndx = sample(1:nrow(temp_tagged_pop), size = 1)
	        counter = 0;
	        while ((temp_tagged_pop$weight[agent_ndx] == 0) | (counter > 5)) {
	          agent_ndx = sample(1:nrow(temp_tagged_pop), size = 1)
	          counter = counter + 1
	        }
          if (counter > 5)
            break
	        agent_ndx_sampled = c(agent_ndx_sampled,agent_ndx)
	        temp_catch = temp_catch - temp_tagged_pop$weight[agent_ndx];
	        temp_tagged_pop$weight[agent_ndx] = 0.0
	        tag_count = tag_count + 1;
	      }
	    }
	    #tag_count 
	    #catch * sum(tagged_pop$weight) / population
	    pop_est = ((n_tags+1)*(catch+1)/(tag_count+1))-1
	    results[bt, scalar_ndx, tag_ndx] = pop_est
	    tag_counts[bt, scalar_ndx, tag_ndx] = tag_count
		}
		tag_ndx = tag_ndx + 1;
	}
	scalar_ndx = scalar_ndx + 1;
}
endTime<- Sys.time() - startTime
# 8.2 mins
bias = estimated = cv = true = matrix(0,nrow = length(scalars), ncol = length(tags_to_look_at))
scalar_ndx = 1;
tag_ndx = 1;
for(scalar_ndx in 1:length(scalars)) {
  n_agents = population / scalars[scalar_ndx]
	for(tag_ndx in 1:length(tags_to_look_at)) {
		rmsd = RMSD(population, results[,scalar_ndx, tag_ndx])
		bias[scalar_ndx,tag_ndx] = rmsd[1,"bias"]
		estimated[scalar_ndx,tag_ndx] = rmsd[1,"est_mean"]
		cv[scalar_ndx,tag_ndx] = rmsd[1,"RMCV"]
		true[scalar_ndx,tag_ndx] = rmsd[1,"true"]
	}
}


dimnames(true) = dimnames(bias) = dimnames(cv) = dimnames(estimated) = list(scalars, tags_to_look_at)
bias
estimated
true
setwd("C:/Work/Software/IBM/ValidationCode")
write.table(bias, file = "Bias.txt", quote = F)
write.table(true, file = "true.txt", quote = F)
write.table(estimated,file = "estimated.txt", quote = F)

## so makes a big difference if small tag population to overall population.
weighted_Random_Sample <- function(
    .data,
    .weights,
    .n
    ){

    key <- runif(length(.data)) ^ (1 / .weights)
    return(.data[order(key, decreasing=TRUE)][1:.n])
}


## repeat with new algorithm
boots = 500;
set.seed(123)
tags_to_look_at = c(100, 1000, 2000)
catches = c(1000, 5000)
scalars = c(10,50,100)
results = tag_counts = array(0, dim = c(boots,length(catches), length(scalars), length(tags_to_look_at)))
n_tags = 1000
## should possible look at these two

population = 100000
scalar = 1
catch = 2000
startTime<-Sys.time()
scalar_ndx = 1;
tag_ndx = 1;
catch_ndx = 1;
#for(catch in catches) {
	scalar_ndx = 1;
	#for(scalar in scalars) {
		print(scalar)
		tag_ndx = 1;
		#for(n_tags in tags_to_look_at) {
			print(paste0("tag = ", n_tags))
		  n_agents = as.integer(population / scalar)
			population = n_agents * scalar
			tags = rep(0,n_agents)
			weights = rep(scalar, n_agents)
			# apply tagging, all agents have equal prob of being tagged
			tag_agent_ndx = sample(x = n_agents, size = n_tags, replace = F)
			## pull out tagged individuals take away the scalar
			weights[tag_agent_ndx] = weights[tag_agent_ndx] - 1
			## initially append taged agents at end and then reshuffle
			n_new_agents = n_agents + n_tags
			weights = c(weights, rep(1,n_tags))
			tags = c(tags, rep(1,n_tags))
			reshuffle_ndx = sample(1:n_new_agents)
			weights = weights[reshuffle_ndx]
			prob_untagged = sum(weights[tags == 0]) / sum(weights)
			prob_tagged = sum(weights[tags == 1]) / sum(weights)
			#full_prob = c(prob_tagged, prob_untagged)
			prob = weights / sum(weights)
			tags = tags[reshuffle_ndx]
			pop_est = NULL
			for(bt in 1:boots) {
				temp_weights = weights
				temp_catch = catch
				## now do recapture event first treating this is as equal probability
				tag_count = 0
				while(temp_catch > 0) {
					individual_ndx = round(runif(1, 1, sum(temp_weights, na.rm = T)))
					this_agent_ndx = which.min(abs(cumsum(temp_weights) - individual_ndx))
					if (temp_weights[this_agent_ndx] == 0)
					  next
					if (tags[this_agent_ndx] == 1)
					  tag_count = tag_count + 1;
					temp_catch = temp_catch - temp_weights[this_agent_ndx]
					temp_weights[this_agent_ndx] = 0
				}	
				print(tag_count)
				pop_est = ((n_tags+1)*(catch+1)/(tag_count+1))-1
				results[bt,catch_ndx, scalar_ndx, tag_ndx] = pop_est
				tag_counts[bt,catch_ndx, scalar_ndx, tag_ndx] = tag_count
			}
			#tag_ndx = tag_ndx + 1;
		#}
		#scalar_ndx = scalar_ndx + 1;
	
	#catch_ndx = catch_ndx + 1;
#}
endTime<- Sys.time() - startTime
endTime

## run with some specific settings
# -scalar = 1, catch = 2000, pop = 1e+5, n_tags = 1000, boots = 500 summary(tags) = 
summary(tag_counts[,catch_ndx, scalar_ndx, tag_ndx])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8.00   17.00   20.00   20.04   23.00   33.00 

## we should get something similar, with scalar = 10
# -scalar = 10, catch = 2000, pop = 1e+5, n_tags = 1000, boots = 500 summary(tags) = 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00    9.00   11.00   11.54   14.00   23.00 

# 8.2 mins
bias = estimated = cv = true = matrix(NA,nrow = length(scalars), ncol = length(tags_to_look_at))
scalar_ndx = 1;
tag_ndx = 1;
## for the first catch ndx
for(scalar_ndx in 1:length(scalars)) {
	population = n_agents * scalars[scalar_ndx]
	for(tag_ndx in 1:length(tags_to_look_at)) {
		rmsd = RMSD(population, results[,1,scalar_ndx, tag_ndx])
		bias[scalar_ndx,tag_ndx] = rmsd[1,"bias"]
		estimated[scalar_ndx,tag_ndx] = rmsd[1,"est_mean"]
		cv[scalar_ndx,tag_ndx] = rmsd[1,"RMCV"]
		true[scalar_ndx,tag_ndx] = rmsd[1,"true"]
	}
}
dimnames(true) = dimnames(bias) = dimnames(cv) = dimnames(estimated) = list(scalars, tags_to_look_at)
bias






## repeat with new algorithm
boots = 500;
set.seed(125)
tags_to_look_at = c(100, 1000, 2000)
catches = c(1000, 5000)
scalars = c(10,50,100)
results = tag_counts = array(0, dim = c(boots,length(catches), length(scalars), length(tags_to_look_at)))
n_tags = 1000
## should possible look at these two
population = 100000
scalar = 10
catch = 2000
startTime<-Sys.time()
scalar_ndx = 1;
tag_ndx = 1;
catch_ndx = 1;
#for(catch in catches) {
scalar_ndx = 1;
#for(scalar in scalars) {
print(scalar)
tag_ndx = 1;
#for(n_tags in tags_to_look_at) {
print(paste0("tag = ", n_tags))
n_agents = as.integer(population / scalar)
population = n_agents * scalar
tags = rep(0,n_agents)
weights = rep(scalar, n_agents)
# apply tagging, all agents have equal prob of being tagged
tag_agent_ndx = sample(x = n_agents, size = n_tags, replace = F)
## pull out tagged individuals take away the scalar
weights[tag_agent_ndx] = weights[tag_agent_ndx] - 1
## initially append taged agents at end and then reshuffle
n_new_agents = n_agents + n_tags
weights = c(weights, rep(1,n_tags))
tags = c(tags, rep(1,n_tags))
reshuffle_ndx = sample(1:n_new_agents)
weights = weights[reshuffle_ndx]
tags = tags[reshuffle_ndx]
prob_untagged = sum(weights[tags == 0]) / sum(weights)
prob_tagged = sum(weights[tags == 1]) / sum(weights)
## account for sampling agents not inidviduals
prob_untagged = prob_untagged / (scalar)
prob_untagged = prob_untagged / (prob_untagged + prob_tagged)
prob_tagged = 1 - prob_untagged
#full_prob = c(prob_tagged, prob_untagged)
pop_est = NULL
for(bt in 1:boots) {
  temp_weights = weights
  temp_catch = catch
  ## now do recapture event first treating this is as equal probability
  tag_count = 0
  while(temp_catch > 0) {
    agent_tag = rbinom(n = c(1), size = 1, prob = prob_tagged)
    ## now find agent and sample
    for(i in 1:length(tags)) { ## could be issues if this vector isn't random.
      if (tags[i] == agent_tag & temp_weights[i] != 0) {
        temp_catch = temp_catch - temp_weights[i]
        temp_weights[i] = 0
        if (tags[i] == 1)
          tag_count = tag_count + 1;
        break;
      }
    }
  }
  pop_est = ((n_tags+1)*(catch+1)/(tag_count+1))-1
  results[bt,catch_ndx, scalar_ndx, tag_ndx] = pop_est
  tag_counts[bt,catch_ndx, scalar_ndx, tag_ndx] = tag_count
}
#tag_ndx = tag_ndx + 1;
#}
#scalar_ndx = scalar_ndx + 1;

#catch_ndx = catch_ndx + 1;
#}
endTime<- Sys.time() - startTime
endTime

summary(tag_counts[,catch_ndx, scalar_ndx, tag_ndx])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8.00   17.00   20.00   20.31   23.00   36.00 

