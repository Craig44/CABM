##########################################
# A little experiment on how to move individuals around box's (Box transfer)
##########################################
N = 1000000 ## in cell 1
prob = c(0.333333,0.033333,0.6333333)
area = c(1,2,3)
## prob = probabiliity of moving to another cell
area_freq = area_move2 = area_move = vector(length = 3, mode = "numeric");
for (i in 1:N) {
  ## Current implementation
  possible_area = sample(area,1)
  area_freq[possible_area] = area_freq[possible_area] + 1;
  if( runif(1) <= prob[possible_area]) {
    area_move[possible_area] = area_move[possible_area] + 1;
  }
  ## proposed implementation
  ndx = which(as.numeric(rmultinom(1,size = 1, prob)) > 0)
  area_move2[ndx] = area_move2[ndx] + 1;
}
area_move2 / sum(area_move2)
area_move / sum(area_move)
N - sum(area_move2)
N - sum(area_move)

# standard approach for sampling from a multinomial like this is to compare a random standard 
# uniform value to the cumulative sums of the probabilities and return the first index for which 
# the cumulative sum is greater than the random uniform
