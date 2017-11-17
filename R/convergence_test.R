#' testing_SLE_convergence() function
#'
#' This function perform convergence test for SLE, and output in a dataset full of sampling data
#' The argument data can only take in a SLE_object produced by SLE()
#'
#' @param iter indicates number of iteration
#' @param data  must take in a SLE_object produced by SLE()
#' @param ...
#' @return A list of Convergence_Test_object
#' @keywords data, iter
#' @export
testing_SLE_convergence<-function(iter=100,
                                  data#,#= SLE_object, # must take in a SLE_object produced by SLE() 
                                  #filepath=paste0(getwd(),'/R')
								  ){ # remember to save SE output as SLE_object, for sake of later use
  #nTrials
  options(scipen = 999)
  if(data$uncertainty_output==FALSE){
    stop('\n
         SLE_object does not specify the analysis of uncertainty\n
         Please look for the argument "uncertainty_output" in SLE() and set it to TRUE')
  }
  # filepath,'/SLE_related_function',
  #if(file.exists(paste0('./prior_for_SLE.R'))==FALSE){
  #  stop('\n
  #       Missing prior_for_SLE.R file, cannot estimate priors,\n
  #       Please check the integrity of source files')
  #} #filepath,'/SLE_related_function',
  #source(paste0('/prior_for_SLE.R'))
  age_labels<-data$Original_data$Age_groups
  population<-data$Original_data$population
  death_counts<-data$Original_data$death
  standard_age<-data$Original_life_table$Age_group
  prior_population<-data$priors_info
  
  up_bound<-unlist(lapply(2:length(age_labels),function(i) standard_age[1:length(age_labels)][i]-0.1))
  uniform_prior_range<-data.frame(lr=standard_age[1:length(age_labels)],
                                  up=c(up_bound,100))

  LE_samples<-list()
  for(i in 1:length(age_labels)){
    sampling_outcome<-matrix(c(replicate(iter,expr = population_priors(x=population[i],scale = prior_population)),
                               replicate(iter,expr = death_priors(death_counts[i])),
                               replicate(iter,expr = AAAdeath(uniform_prior_range,i))),nrow=iter,byrow=FALSE)
    dimnames(sampling_outcome)[[2]]<-c('Population','counts','average_age_at_death')
    LE_samples[[i]]<-sampling_outcome
  }
  class(LE_samples) <- "Convergence_Test_object"
  attr(LE_samples,"Call") <- sys.call()
  return(LE_samples)
}
#testing<-c()
#for(i in 1:nTrials){
#  testing<-c(testing,mean(testing_convergence(i)[[1]][,1]))
#}

