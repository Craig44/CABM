#September 8th 2020 J McKenzie


#############
# FUNCTIONS #
#############               

#RMSD



RMSD<- function(x,y){
  
  
  CV<-round(sqrt( sum((x-y)^2)/length(y))/x,2)
  bias<-round((mean(y)/x)-1,2)
  return(cbind('true'=x,'est_mean'=round(mean(y),2),'bias'=bias,'RMCV'=CV,'obs'= length(y)))
}



LFPlot<- function(den1,den2){
  pmax<- max(c(den1[,3],den2[,3]))
  pmin<- min(c(den1[,3],den2[,3]))
  plot(0,0,type = 'n',ylab='proportional freq',xlab='length cm', 
       xlim=c(min(den1$length),max(den1$length)),ylim =c(pmin,pmax*1.01))
  
  lines(den1$length,den1$pfreq, col="black",lwd=2)
  points(den2$length,den2$pfreq, col="blue")
  
  mtext(side=3,adj=0,line=0.3,paste("Sample proportional frequency"),cex=par()$cex*1.1)
}






################
# Main Routine #
################

#Creat lf density
meanlog = 4

sdlog = 0.5

freq_dist<-round(rlnorm(1000000, meanlog, sdlog),0)
den<-aggregate(freq_dist,list(freq_dist),length)
names(den)<-c('length','freq')
den$pfreq<- den$freq/sum(den$freq)
den<-subset(den, den$pfreq >0.005)
den$pfreq<- den$freq/sum(den$freq)
freq_dist<- subset(freq_dist,freq_dist >= min(den$length) & freq_dist <= max(den$length))






pop = 100000
agents = 1000





# Number of tags
tags = 500

# Number of catch
catch = 2500

boot<- 100
#######################
# agent scaler value
mscalr = pop/agents



RunDet<-list()


RunDet[c('pop','mscalr','tags','catch','boot','runtime','stats')]<-c(pop,mscalr,tags,catch,boot,0,0)






pop_est<-NULL
startTime<-Sys.time()
bt =1
for(bt in 1:boot){  

  #### Create agent matrix
  agent_vec<- rep(mscalr,pop/mscalr)

  
  agent_mat<-as.data.frame(cbind('agent'=1:length(agent_vec), 'scale'=agent_vec, 'tag'=0,'length'=0))
  for(i in 1:length(agent_mat[,1])){
    agent_mat[i,4] <- sample(freq_dist, 1, replace = FALSE, prob = NULL)
  }


  
  densc <- aggregate(agent_mat[,2],list(agent_mat[,4]),sum)
  names(densc)<-c('length','freq')
  densc$pfreq<-densc$freq/sum(densc$freq)
  
  #LFPlot(den,densc)
  
  
  ##########################################
  
  # Tagging sequence
  
  for( i in 1:tags){
    notag = TRUE
    while(notag){
      point<-sample(agent_mat[,1], 1, replace = FALSE, prob = NULL)
      #if(agent_mat[point,3] < as.integer(agent_mat[point,2])){
      if(agent_mat[point,3] < round(agent_mat[point,2],0)){        
        
        agent_mat[point,3] <-agent_mat[point,3] +1
        notag = FALSE
      }
    }
  }    
    
  dentg<-aggregate(agent_mat[,3],list(agent_mat[,4]),sum)
  names(dentg)<-c('length','freq')
  dentg$pfreq<-dentg$freq/sum(dentg$freq)
  
  #LFPlot(den,dentg)
  hold<-agent_mat

  
  #############################
  # Take catch routines 
  #############################
  
  cat_tag<-catch
  rectag<-NULL
  surplus_catch<-0
  while(cat_tag>0) {
    
    point<-sample(agent_mat[,1], 1, replace = FALSE, prob = NULL)
    cat<-agent_mat[point,2]
    if(cat_tag>=cat)
    {
      if(cat_tag-cat < 1e-10){ # deals with floating point inequality
         cat_tag <-0
      }else{
      
        cat_tag<-cat_tag-cat
      }   
           
      
      
      #recover tags
      if(agent_mat[point,3]>0){
          rectag<-rbind(rectag,agent_mat[point,3:4])
      }
      
      
      agent_mat<-agent_mat[-point,]
      if(length(agent_mat[,1])>0){
         agent_mat[,1]<-seq(1,length(agent_mat[,1]),1)
      }   
    }else{
      
      
      #recover tags
      if(agent_mat[point,3]>0){
        rectag<-rbind(rectag,agent_mat[point,3:4])
      }
      surplus_catch<- cat - cat_tag
      
      cat_tag<-0
      agent_mat<-agent_mat[-point,]
      if(length(agent_mat[,1]>0)){ 
        agent_mat[,1]<-seq(1,length(agent_mat[,1]),1)
      }
    }  
    cat_tag
  }
  
  
  # new_pop<-sum(agent_mat[,2])
  #stats
  tagrec<-sum(rectag[,1])
  pop_est<- rbind(pop_est,cbind('pop_est'= ((tags+1)*(catch+1+surplus_catch)/(tagrec+1))-1, 'tagsrec'=tagrec)) # Petersen estimator with Chapman bias correction

  pop_est<-pop_est[order(pop_est[,1]),]
} #endboot


endTime<-Sys.time()
runTime<-endTime-startTime


# stats


RunDet[['runtime']]<-runTime
RunDet[['stats']]<-RMSD(pop,pop_est[,1])


RunDet









