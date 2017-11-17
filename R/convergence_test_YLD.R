#' testing_YLD_convergence() function
#'
#' This function perform convergence test for YLD, and output in a dataset full of sampling data
#' The argument data can only take in a YLD_object produced by YLD_incident(), not by YLD_prevalent()  
#'
#' @param iter indicates number of iteration
#' @param data  must take in a YLD_object produced by YLD_incident(), not by YLD_prevalent()
#' @param ...
#' @return A list of Convergence_Test_object
#' @keywords data, iter
#' @export
testing_YLD_convergence<-function(iter=100, 
                                  data#,#= YLD_range, # must take in a YLD_object produced by YLD_incident(), not YLD_prevalent
                                  #filepath=paste0(getwd(),'/R')
								  ){ # remember to save YLD output as YLD_object, for sake of later use
  #nTrials
  options(scipen = 999)
  if(data$YLD_uncertainty$Calculate_range==FALSE){
    stop('\n
         YLD_object does not specify the analysis of uncertainty\n
         Please look for the argument "uncertainty_output" in SLE() and set it to TRUE')
  }
  # filepath,'/YLD_related_function',
  #if(file.exists(paste0('./YLD_priors.R'))==FALSE){
  #  stop('\n
  #       Missing YLD_priors.R file, cannot estimate priors,\n
  #       Please check the integrity of source files')
  #} # filepath,'/YLD_related_function',
  #source(paste0('./YLD_priors.R'))
  age_labels<-data$Original_data$Age
  population<-data$Original_data$population
  incidence_forcal<-data$Original_data$incidence
  Duration_interval<-data$input_interval$Duration_interval
  DisabilityWeight_interval<-data$input_interval$DisabilityWeight_interval
  prior_population<-data$priors_info
  
  
  #+++++
    
  YLD_priors_samples<-list()
  for(i in 1:length(age_labels)){
    sampling_outcome<-matrix(c(replicate(iter,expr = population_priors(x=population[i],scale = prior_population)),
                               replicate(iter,expr = death_priors(incidence_forcal[i])),
                               replicate(iter,expr = Durantion_priors(Duration_interval,i)),
                               replicate(iter,expr = DisabilityWeight_priors(DisabilityWeight_interval,i))),nrow=iter,byrow=FALSE)
    dimnames(sampling_outcome)[[2]]<-c('Population','Incidence','Duration','Disability_weight')
    YLD_priors_samples[[i]]<-sampling_outcome
  }
  #+++++
  
  #LE_samples<-list()
  #for(i in 1:length(age_labels)){
  #  sampling_outcome<-matrix(c(replicate(iter,expr = population_priors(x=population[i],scale = prior_population)),
  #                             replicate(iter,expr = death_priors(incident_counts[i])),
  #                             replicate(iter,expr = AAAdeath(uniform_prior_range,i))),nrow=iter,byrow=FALSE)
  #  dimnames(sampling_outcome)[[2]]<-c('Population','counts','average_age_at_death')
  #  LE_samples[[i]]<-sampling_outcome
  #}
  class(YLD_priors_samples) <- "Convergence_Test_object"
  attr(YLD_priors_samples,"Call") <- sys.call()
  return(YLD_priors_samples)
}


#' testing_prevYLD_convergence() function
#'
#' This function perform convergence test for YLD, and output in a dataset full of sampling data
#' The argument data can only take in a YLD_object produced by YLD_prevalent(), not by YLD_incident()
#'
#' @param iter indicates number of iteration
#' @param data  must take in a YLD_object produced by YLD_prevalent(), not by YLD_incident()
#' @param ...
#' @return A list of Convergence_Test_object
#' @keywords data, iter
#' @export
testing_prevYLD_convergence<-function(iter=100, 
                                  data#,#= YLD_range, # must take in a YLD_object produced by YLD_prevalent
                                  #filepath=paste0(getwd(),'/R')
								  ){ # remember to save YLD output as YLD_object, for sake of later use
  #nTrials
  options(scipen = 999)
  if(data$YLD_uncertainty$Calculate_range==FALSE){
    stop('\n
         YLD_object does not specify the analysis of uncertainty\n
         Please look for the argument "uncertainty_output" in SLE() and set it to TRUE')
  }
   # filepath,'/YLD_related_function',
  #if(file.exists(paste0('./YLD_priors.R'))==FALSE){
  #  stop('\n
  #       Missing YLD_priors.R file, cannot estimate priors,\n
  #       Please check the integrity of source files')
  #} # filepath,'/YLD_related_function',
  #source(paste0('./YLD_priors.R'))
  age_labels<-data$Original_data$Age
  population<-data$Original_data$population
  prevalence_forcal<-data$Original_data$prevalence
  DisabilityWeight_interval<-data$input_interval$DisabilityWeight_interval
  prior_population<-data$priors_info
  
  
  #+++++
  
  YLD_priors_samples<-list()
  for(i in 1:length(age_labels)){
    sampling_outcome<-matrix(c(replicate(iter,expr = population_priors(x=population[i],scale = prior_population)),
                               replicate(iter,expr = death_priors(prevalence_forcal[i])),
                               replicate(iter,expr = DisabilityWeight_priors(DisabilityWeight_interval,i))),nrow=iter,byrow=FALSE)
    dimnames(sampling_outcome)[[2]]<-c('Population','Prevalence','Disability_weight')
    YLD_priors_samples[[i]]<-sampling_outcome
  }
  #+++++
  
  #LE_samples<-list()
  #for(i in 1:length(age_labels)){
  #  sampling_outcome<-matrix(c(replicate(iter,expr = population_priors(x=population[i],scale = prior_population)),
  #                             replicate(iter,expr = death_priors(incident_counts[i])),
  #                             replicate(iter,expr = AAAdeath(uniform_prior_range,i))),nrow=iter,byrow=FALSE)
  #  dimnames(sampling_outcome)[[2]]<-c('Population','counts','average_age_at_death')
  #  LE_samples[[i]]<-sampling_outcome
  #}
  class(YLD_priors_samples) <- "Convergence_Test_object"
  attr(YLD_priors_samples,"Call") <- sys.call()
  return(YLD_priors_samples)
  }