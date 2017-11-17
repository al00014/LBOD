YLD_prevalent<-function(age_labels,
                        population,
                        prevalence,
                        DisabilityWeight,
                        DisabilityWeight_interval,
                        RateUNIT=1000,
                        YLD_perUnit=1000,
                        Input_Prevalent_Rate=TRUE,
                        nTrials=2000,
                        prior_population=0.001,
                        uncertainty_range=TRUE, ## in calculating uncertainty for YLD, CI for age_weight should be accounted for
                        uncertainty_alpha=0.05,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                        max_age=90#,
                        #filepath=paste0(getwd(),'/R')
   
){
  options(scipen = 999)
  #if(sum(!c('parallel','triangle') %in% rownames(installed.packages()) )!= 0) {
  #  stop('\n
  #        Please installed the dependencies first:\n
  #       "parallel","triangle" ')
  #}
  
  
  if(Input_Prevalent_Rate==TRUE){
    prevalence_forcal<-population*prevalence/RateUNIT
  } else if(Input_Prevalent_Rate==FALSE){
    prevalence_forcal<-prevalence
  }
  
  data_sort<-data.frame(Age=age_labels,
                        population=population,
                        prevalence=prevalence_forcal,
                        DisabilityWeight=DisabilityWeight)
  
  if(uncertainty_range==TRUE){
    
    #if(length(age_ForTurningonset)!=length(age_labels)){
    #  stop('Make sure the input values of age_ForTurningonset is consistent with age_labels,\n
    #       for the uncertainty analysis of YLD requires the specific use of age_ForTurningonset')
    #}
    
    
    ### for uncertainty of YLD, mainly consider the priors for incidence, duration, disability weight
    ### this can be supported by Australian BOD study, report link http://www.aihw.gov.au/WorkArea/DownloadAsset.aspx?id=60129547710 at page 142
    
    ### so no need to consider uncertainty for age_at_onset, 
    ### this variable happens to be the mean of each age group (age group is defined by interval, so mean of the interval) 
    
    #age_ForTurningonset<-as.numeric(age_ForTurningonset)
    #up_bound<-unlist(parallel::mclapply(2:length(age_labels),function(i) age_ForTurningonset[1:length(age_labels)][i]-0.1))
    
    #uniform_prior_range_age_at_onset<-data.frame(lr=age_ForTurningonset[1:length(age_labels)],
    #                                up=c(up_bound,max_age))
    #round(rowMeans(uniform_prior_range),digits = 1)
    
    # to comply with uncertainty analysis of YLL, uncertainty for population variable is also considered
    
    source(paste0('./YLD_priors.R'))
    # filepath,'/YLD_related_function',
    YLD_priors_samples<-list()
    for(i in 1:length(age_labels)){
      sampling_outcome<-matrix(c(replicate(nTrials,expr = population_priors(x=population[i],scale = prior_population)),
                                 replicate(nTrials,expr = death_priors(prevalence_forcal[i])),
                                 replicate(nTrials,expr = DisabilityWeight_priors(DisabilityWeight_interval,i))),nrow=nTrials,byrow=FALSE)
      dimnames(sampling_outcome)[[2]]<-c('Population','Prevalence','Disability_weight')
      YLD_priors_samples[[i]]<-sampling_outcome
    }
    
    uncertainty_results<-parallel::mclapply(1:length(YLD_priors_samples), function(i){
      tempdata<-apply(YLD_priors_samples[[i]],2,quantile,probs=c(uncertainty_alpha/2,1-uncertainty_alpha/2))
      
      tempdata<-as.matrix(tempdata)
      
    })
    
    ## population
    Population_uncertainty<-data.frame(age_labels,population,
                                       round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,1]})),digits = 0))
    names(Population_uncertainty)<-c('Age_groups','Recorded_population','Lr','Up')
    
    ## incidence
    Counts_uncertainty<-data.frame(age_labels,prevalence_forcal,
                                   round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,2]})),digits = 0))
    names(Counts_uncertainty)<-c('Age_groups','Recorded_prevalence','Lr','Up')
    
    
    ## DisabilityWeight
    DW_uncertainty<-data.frame(age_labels,DisabilityWeight,
                               round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,3]})),digits = 1))
    names(DW_uncertainty)<-c('Age_groups','Recorded_DW','Lr','Up')
    
    YLD_point<-Counts_uncertainty$Recorded_prevalence*DW_uncertainty$Recorded_DW
    YLD_Lr<-Counts_uncertainty$Lr*DW_uncertainty$Lr
    YLD_Up<-Counts_uncertainty$Up*DW_uncertainty$Up
    
    ### point est
    population_withTotal<-c(population,sum(population))
    prevalence_forcal_total<-c(prevalence_forcal,sum(prevalence_forcal))
    DisabilityWeight_withTotal<-c(DisabilityWeight,sum(DisabilityWeight*prevalence_forcal)/sum(prevalence_forcal))
    
    YLD_withTotal<-c(YLD_point,sum(YLD_point))
    YLD_per<-YLD_withTotal/population_withTotal*YLD_perUnit
    
    ### lr
    population_withTotal_lr<-c(Population_uncertainty$Lr,sum(Population_uncertainty$Lr))
    prevalence_forcal_total_lr<-c(Counts_uncertainty$Lr,sum(Counts_uncertainty$Lr))
    DisabilityWeight_withTotal_lr<-c(DW_uncertainty$Lr,sum(DW_uncertainty$Lr*Counts_uncertainty$Lr)/sum(Counts_uncertainty$Lr))
    
    YLD_withTotal_lr<-c(YLD_Lr,sum(YLD_Lr))
    YLD_per_lr<-YLD_withTotal_lr/population_withTotal_lr*YLD_perUnit
    
    ### Up
    population_withTotal_up<-c(Population_uncertainty$Up,sum(Population_uncertainty$Up))
    prevalence_forcal_total_up<-c(Counts_uncertainty$Up,sum(Counts_uncertainty$Up))
    DisabilityWeight_withTotal_up<-c(DW_uncertainty$Up,sum(DW_uncertainty$Up*Counts_uncertainty$Up)/sum(Counts_uncertainty$Up))
    
    YLD_withTotal_up<-c(YLD_Up,sum(YLD_Up))
    YLD_per_up<-YLD_withTotal_up/population_withTotal_up*YLD_perUnit
    
    
    YLD_DataFrame<-data.frame(Age=c(as.character(age_labels),'Total'),
                              population=population_withTotal,
                              prevalent_case=prevalence_forcal_total,
                              DisabilityWeight=DisabilityWeight_withTotal,
                              YLD=YLD_withTotal,
                              YLD_per=YLD_per)
    
    YLD_DataFrame_lr<-data.frame(Age=c(as.character(age_labels),'Total'),
                                 population=population_withTotal_lr,
                                 prevalent_case=prevalence_forcal_total_lr,
                                 DisabilityWeight=DisabilityWeight_withTotal_lr,
                                 YLD=YLD_withTotal_lr,
                                 YLD_per=YLD_per_lr)
    YLD_DataFrame_up<-data.frame(Age=c(as.character(age_labels),'Total'),
                                 population=population_withTotal_up,
                                 prevalent_case=prevalence_forcal_total_up,
                                 DisabilityWeight=DisabilityWeight_withTotal_up,
                                 YLD=YLD_withTotal_up,
                                 YLD_per=YLD_per_up)
    
    names(YLD_DataFrame)[ncol(YLD_DataFrame)]<-paste0('YLD','_per',YLD_perUnit)
    names(YLD_DataFrame_lr)[ncol(YLD_DataFrame_lr)]<-paste0('YLD','_lr','_per',YLD_perUnit)
    names(YLD_DataFrame_up)[ncol(YLD_DataFrame_up)]<-paste0('YLD','_up','_per',YLD_perUnit)
    
    
    Original_data<-data.frame(Age=c(as.character(age_labels),'Total'),
                              population=population_withTotal,
                              prevalence=prevalence_forcal_total,
                              prevalence_rate=prevalence_forcal_total/population_withTotal*RateUNIT,
                              DisabilityWeight=DisabilityWeight_withTotal)
    names(Original_data)[4]<-paste0( "prevalence_rate" ,'_per',RateUNIT)
    
    YLD_object<-list(Age_group=c(as.character(age_labels),'Total'),
                     YLD=YLD_withTotal,
                     YLD_lr=YLD_withTotal_lr,
                     YLD_up=YLD_withTotal_up,
                     YLD_per=YLD_per,
                     YLD_per_lr=YLD_per_lr,
                     YLD_per_up=YLD_per_up,
                     perUnit=YLD_perUnit,
                     YLD_uncertainty=list(Calculate_range=uncertainty_range,
                                          Alpha=uncertainty_alpha),
                     priors_info=prior_population,
                     parameters_MC_samples=YLD_priors_samples,
                     uncertainty_results=uncertainty_results, # uncertainty_results list the Uncertainty interval for 
                     # each variable (population, death, average at death) for each specific age group
                     # interval range refers to arguments uncertainty_alpha
                     Original_data_withtotal=Original_data,
                     Original_data=Original_data[-nrow(Original_data),],
                     YLD_DataFrame=YLD_DataFrame,
                     YLD_DataFrame_lr=YLD_DataFrame_lr,
                     YLD_DataFrame_up=YLD_DataFrame_up,
                     UncertaintyFile=list(population=Population_uncertainty,
                                          prevalence_case=Counts_uncertainty,
                                          disability_weight=DW_uncertainty),
                     input_interval=list(DisabilityWeight_interval=DisabilityWeight_interval))
  }else if(uncertainty_range==FALSE){
    YLD_point<-data_sort$prevalence*data_sort$DisabilityWeight
    
    population_withTotal<-c(population,sum(population))
    prevalence_forcal_total<-c(prevalence_forcal,sum(prevalence_forcal))
    DisabilityWeight_withTotal<-c(DisabilityWeight,sum(DisabilityWeight*prevalence_forcal)/sum(prevalence_forcal))
    
    YLD_withTotal<-c(YLD_point,sum(YLD_point))
    YLD_per<-YLD_withTotal/population_withTotal*YLD_perUnit
    
    
    YLD_DataFrame<-data.frame(Age=c(as.character(age_labels),'Total'),
                              population=population_withTotal,
                              prevalence=prevalence_forcal_total,
                              DisabilityWeight=DisabilityWeight_withTotal,
                              YLD=YLD_withTotal,
                              YLD_per=YLD_per)
    
    names(YLD_DataFrame)[ncol(YLD_DataFrame)]<-paste0( "YLD" ,'_per',YLD_perUnit)
    
    Original_data<-data.frame(Age=c(as.character(age_labels),'Total'),
                              population=population_withTotal,
                              prevalence=prevalence_forcal_total,
                              prevalence_rate=prevalence_forcal_total/population_withTotal*RateUNIT,
                              DisabilityWeight=DisabilityWeight_withTotal)
    names(Original_data)[4]<-paste0( "prevalence_rate" ,'_per',RateUNIT)
    
    YLD_object<-list(Age_group=c(as.character(age_labels),'Total'),
                     YLD=YLD_withTotal,
                     YLD_per=YLD_per,
                     #YLD_uncertainty=list(uncertainty_range,uncertainty_alpha),
                     Original_data_withtotal=Original_data,
                     Original_data=Original_data[-nrow(Original_data),],
                     YLD_DataFrame=YLD_DataFrame)
    
  }
  class(YLD_object) <- "Prevalence_YLD_object"
  attr(YLD_object,"Call") <- sys.call()
  return(YLD_object)
}