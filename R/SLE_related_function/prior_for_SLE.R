population_priors<-function(x,scale){ # using hyper priors and norm distribution for uncertainty of population
  SD_priors<-rgamma(1, shape=x,scale = scale) # vague hyper-priors, non-informative
  return(rnorm(1,x,SD_priors))
}

death_priors<-function(x){ # using poisson distribution for death count data
  return(rpois(1,x))
}

AAAdeath<-function(x,index){
  return(round(runif(1,x[index,1],x[index,2]),digits = 1))
}