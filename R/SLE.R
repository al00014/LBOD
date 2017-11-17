
SLE<-function(age_labels,population,death_counts,age_average_at_death,
              standard_LE,standard_age,
              prior_population=0.001,
              nTrials=2000,
              uncertainty_range=TRUE,
              uncertainty_alpha=0.05,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
              max_age=100#,
              #filepath=paste0(getwd(),'/R')
){
  options(scipen = 999)
  #standard_age[1:19]
  #length(standard_age)
  #length(age_average_at_death)
  age_labels<-as.character(age_labels)
  
#  if(sum(!c('xlsx','ggplot2','parallel','tidyr') %in% rownames(installed.packages()) )!= 0) {
#    stop('\n
#         Please installed the dependencies first:
#         xlsx, ggplot2','parallel','tidyr')
#  }
  #require(parallel)
  if (length(death_counts)>=length(standard_LE)){
    stop('\n
         There are less age groups in standard life table than those in the epidemiological data. \n
         Currently not acceptable to these data form.\n
         Check the age groups of both input data, \n
         and try to make age groups in epidemiological data less than those in life table')
  }
  if (class(population)=='character'){
    stop('Please input numeric values for population argument')
  }
  if (class(death_counts)=='character'){
    stop('Please input numeric values for death_counts argument')
  }
  if (class(age_average_at_death)=='character'){
    stop('Please input numeric values for age_average_at_death argument')
  }
  if (class(standard_LE)=='character'){
    stop('Please input numeric values for standard_LE argument')
  }
  if (class(standard_age)=='character'){
    stop('Please input numeric values for standard_age argument')
  }
  data_plugin<-data.frame(age_labels,population,death_counts,age_average_at_death)
  life_table_plugin<-data.frame(standard_age,standard_LE)
  names(life_table_plugin)<-c('Age','LE')
  names(data_plugin)<-c('Age','population','death','age_av_death')
  
  #Stand_life_expectancy<-data.frame()
  #for(j in 1:length(data_plugin$age_av_death)){
  #  LE<-life_table_plugin$LE[j]+
  #    (life_table_plugin$LE[j+1]-life_table_plugin$LE[j])*
  #    (data_plugin$age_av_death[j]-
  #       life_table_plugin$Age[j])/(life_table_plugin$Age[j+1]-
  #                                    life_table_plugin$Age[j])
  #  
  #  LE<-data.frame(LE)
  #  Stand_life_expectancy<-rbind(Stand_life_expectancy,LE)
  #}
  
  if(uncertainty_range==TRUE){
    
    #filepath,'/SLE_related_function',
    #if(file.exists(paste0('./prior_for_SLE.R'))==FALSE){
    #  stop('\n
    #     Missing prior_for_SLE.R file, cannot estimate priors,\n
    #     Please check the integrity of source files')
    #}
    #source(paste0('./prior_for_SLE.R'))
    # filepath,'/SLE_related_function',
    
    #population_priors<-function(x,scale=prior_population){ # using hyper priors and norm distribution for uncertainty of population
    #  SD_priors<-rgamma(1, shape=x,scale = scale) # vague hyper-priors, non-informative
    #  return(rnorm(1,x,SD_priors))
    #}
    #replicate(10,expr = population_priors(population[2]))
    
    #death_priors<-function(x){ # using poisson distribution for death count data
    #  return(rpois(1,x))
    #}
    #AAAdeath<-function(x,index){
    #  return(round(runif(1,x[index,1],x[index,2]),digits = 1))
    #}
    
    up_bound<-unlist(parallel::mclapply(2:length(age_labels),function(i) standard_age[1:length(age_labels)][i]-0.1))
    
    uniform_prior_range<-data.frame(lr=standard_age[1:length(age_labels)],
                                    up=c(up_bound,max_age))

    LE_samples<-list()
    for(i in 1:length(age_labels)){
      sampling_outcome<-matrix(c(replicate(nTrials,expr = population_priors(x=population[i],scale = prior_population)),
                                 replicate(nTrials,expr = death_priors(death_counts[i])),
                                 replicate(nTrials,expr = AAAdeath(uniform_prior_range,i))),nrow=nTrials,byrow=FALSE)
      dimnames(sampling_outcome)[[2]]<-c('Population','counts','average_age_at_death')
      LE_samples[[i]]<-sampling_outcome
    }
    
    
    
    uncertainty_results<-parallel::mclapply(1:length(LE_samples), function(i){
      tempdata<-apply(LE_samples[[i]],2,quantile,probs=c(uncertainty_alpha/2,1-uncertainty_alpha/2))
      
      tempdata<-as.matrix(tempdata)
      
    })
    
    Population_uncertainty<-data.frame(age_labels,population,
                                       round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,1]})),digits = 0))
    names(Population_uncertainty)<-c('Age_groups','Recorded_population','Lr','Up')
    
    Counts_uncertainty<-data.frame(age_labels,death_counts,
                                   round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,2]})),digits = 0))
    names(Counts_uncertainty)<-c('Age_groups','Recorded_death','Lr','Up')
    
    Aver_uncertainty<-data.frame(age_labels,age_average_at_death,
                                 round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,3]})),digits = 1))
    names(Aver_uncertainty)<-c('Age_groups','Recorded_average_age_at_death','Lr','Up')
    
    ### getting stand life expectancy for calculation of YLL
    Stand_life_expectancy<-parallel::mclapply(1:length(data_plugin$age_av_death),FUN = function(j){
      LE<-life_table_plugin$LE[j]+
        (life_table_plugin$LE[j+1]-life_table_plugin$LE[j])*
        (data_plugin$age_av_death[j]-
           life_table_plugin$Age[j])/(life_table_plugin$Age[j+1]-
                                        life_table_plugin$Age[j])
    })
    
    Stand_life_expectancy<-unlist(Stand_life_expectancy)
    Stand_life_expectancy<-data.frame(age_labels,Stand_life_expectancy)
    names(Stand_life_expectancy)<-c('Age','LE')
    
    ## SLE lower bound
    
    
    
    data_plugin_lr<-data.frame(Age=age_labels,
                               population=Population_uncertainty$Lr,
                               death=Counts_uncertainty$Lr,
                               age_av_death=Aver_uncertainty$Lr)
    
    
    Stand_life_expectancy_lr<-parallel::mclapply(1:length(data_plugin_lr$age_av_death),FUN = function(j){
      LE<-life_table_plugin$LE[j]+
        (life_table_plugin$LE[j+1]-life_table_plugin$LE[j])*
        (data_plugin_lr$age_av_death[j]-
           life_table_plugin$Age[j])/(life_table_plugin$Age[j+1]-
                                        life_table_plugin$Age[j])
    })
    
    Stand_life_expectancy_lr<-unlist(Stand_life_expectancy_lr)
    
    Stand_life_expectancy_lr<-data.frame(age_labels,Stand_life_expectancy_lr)
    names(Stand_life_expectancy_lr)<-c('Age','LE')
    
    ## SLE upper bound
    
    
    data_plugin_up<-data.frame(Age=age_labels,
                               population=Population_uncertainty$Up,
                               death=Counts_uncertainty$Up,
                               age_av_death=Aver_uncertainty$Up)
    
    
    Stand_life_expectancy_up<-parallel::mclapply(1:length(data_plugin_up$age_av_death),FUN = function(j){
      LE<-life_table_plugin$LE[j]+
        (life_table_plugin$LE[j+1]-life_table_plugin$LE[j])*
        (data_plugin_up$age_av_death[j]-
           life_table_plugin$Age[j])/(life_table_plugin$Age[j+1]-
                                        life_table_plugin$Age[j])
    })
    
    Stand_life_expectancy_up<-unlist(Stand_life_expectancy_up)
    
    Stand_life_expectancy_up<-data.frame(age_labels,Stand_life_expectancy_up)
    names(Stand_life_expectancy_up)<-c('Age','LE')
    
    #SLE_upper<-Stand_life_expectancy_up
    #SLE_lower<-Stand_life_expectancy_lr
    
    # this should be the correct way to do the function
    SLE_lower<-Stand_life_expectancy_up
    SLE_upper<-Stand_life_expectancy_lr
    
    SLE_lower<-unlist(parallel::mclapply(1:nrow(Stand_life_expectancy),FUN=function(i){
      ifelse(is.na(Stand_life_expectancy$LE[i])==TRUE,NA,SLE_lower$LE[i])
    }))
    SLE_upper<-unlist(parallel::mclapply(1:nrow(Stand_life_expectancy),FUN=function(i){
      ifelse(is.na(Stand_life_expectancy$LE[i])==TRUE,NA,SLE_upper$LE[i])
    }))
    
    Aver_uncertainty$Lr<-unlist(parallel::mclapply(1:nrow(Aver_uncertainty),FUN=function(i){
      ifelse(is.na(Aver_uncertainty$Recorded_average_age_at_death[i])==TRUE,NA,Aver_uncertainty$Lr[i])
    }))
    
    Aver_uncertainty$Up<-unlist(parallel::mclapply(1:nrow(Aver_uncertainty),FUN=function(i){
      ifelse(is.na(Aver_uncertainty$Recorded_average_age_at_death[i])==TRUE,NA,Aver_uncertainty$Up[i])
    }))
    
    ####
    SLE_object<-list(SLE=Stand_life_expectancy,
                     SLE_lr=data.frame(Age=age_labels,
                                       LE=SLE_lower),
                     SLE_up=data.frame(Age=age_labels,
                                       LE=SLE_upper),
                     uncertainty_output=uncertainty_range, 
                     alpha=uncertainty_alpha,
                     uncertainty_results=uncertainty_results, # uncertainty_results list the Uncertainty interval for 
                                                              # each variable (population, death, average at death) for each specific age group
                                                              # interval range refers to arguments uncertainty_alpha
                     LE_MC_samples=LE_samples,
                     Population_uncertainty=Population_uncertainty,
                     Death_uncertainty=Counts_uncertainty,
                     Average_at_death_uncertainty=Aver_uncertainty,
                     Original_data=data.frame(Age_groups=age_labels,
                                              population=population,
                                              death=death_counts,
                                              age_av_death=age_average_at_death),
                     Original_life_table=data.frame(Age_group=standard_age,
                                                    Life_expectancy=standard_LE),
                     priors_info=prior_population)
    
  } else if(uncertainty_range==FALSE){
    
    ### getting stand life expectancy for calculation of YLL
    
    Stand_life_expectancy<-parallel::mclapply(1:length(data_plugin$age_av_death),FUN = function(j){
      LE<-life_table_plugin$LE[j]+
        (life_table_plugin$LE[j+1]-life_table_plugin$LE[j])*
        (data_plugin$age_av_death[j]-
           life_table_plugin$Age[j])/(life_table_plugin$Age[j+1]-
                                        life_table_plugin$Age[j])
    })
    
    Stand_life_expectancy<-unlist(Stand_life_expectancy)
    Stand_life_expectancy<-data.frame(age_labels,Stand_life_expectancy)
    names(Stand_life_expectancy)<-c('Age','LE')
    
    SLE_object<-list(SLE=Stand_life_expectancy,uncertainty_output=uncertainty_range,
                     Original_data=data.frame(Age_groups=age_labels,
                                              population=population,
                                              death=death_counts,
                                              age_av_death=age_average_at_death),
                     Original_life_table=data.frame(Age_group=standard_age,
                                                    Life_expectancy=standard_LE),
                     priors_info=prior_population)
  }
  
  
  
  
  class(SLE_object) <- "Standard_Life_Expectancy_object"
  attr(SLE_object,"Call") <- sys.call()
  
  
  return(SLE_object)
  
  
  }