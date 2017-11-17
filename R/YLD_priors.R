population_priors<-function(x,scale){ # using hyper priors and norm distribution for uncertainty of population
  SD_priors<-rgamma(1, shape=x,scale = scale) # vague hyper-priors, non-informative
  return(rnorm(1,x,SD_priors))
}

death_priors<-function(x){ # using poisson distribution for incident data
  return(rpois(1,x))
}
#require(triangle)
Durantion_priors<-function(x,index){
  #x=Duration_interval
  #index=8
  return(round(triangle::rtriangle(n=1,a=x[index,1],b=x[index,2]),digits = 2))
}

DisabilityWeight_priors<-function(x,index){
  #x=Duration_interval
  #index=8
  return(round(triangle::rtriangle(n=1,a=x[index,1],b=x[index,2]),digits = 2))
}

