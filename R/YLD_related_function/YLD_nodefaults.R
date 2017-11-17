YLD_incident<-function(age_labels,
                       population,
                       incidence,
                       age_at_onset,
                       duration,
                       DisabilityWeight,
                       Rate,#=0.03,
                       Beta,#=0.04,
                       Const,#=0.1658,
                       Agewt,#=0,
                       RateUNIT,#=1000,
                       YLD_perUnit,#=1000,
                       Input_Incident_Rate,#=TRUE,
                       nTrials,#=2000,
                       prior_population,#=0.001,
                       uncertainty_range,#=TRUE, ## in calculating uncertainty for YLD, CI for age_weight should be accounted for
                       uncertainty_alpha,#=0.05,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                       #age_ForTurningonset,
                       Duration_interval,
                       DisabilityWeight_interval,
                       max_age,#=90,
                       filepath#='A:/r references/R functions/Disease_Burden_projection_project'
){
  population<-population
  incidence<-incidence
  age_labels<-age_labels
  duration<-duration
  DisabilityWeight<-DisabilityWeight
  
  #warning(str(population))
  #warning(str(incidence))
  #warning(str(age_labels))
  #warning(str(duration))
  #warning(str(DisabilityWeight))
  
  
  options(scipen = 999)
  if(sum(!c('parallel','triangle') %in% rownames(installed.packages()) )!= 0) {
    stop('\n
          Please installed the dependencies first:\n
         "parallel","triangle" ')
  }
  
  
  if(Input_Incident_Rate==TRUE){
    incidence_forcal<-population*incidence/RateUNIT
  } else if(Input_Incident_Rate==FALSE){
    incidence_forcal<-incidence
  }
  BplusR= -(Beta+Rate)
  data_sort<-data.frame(Age=as.character(age_labels),
                          population=population,
                          incidence=incidence_forcal,
                          age_at_onset=age_at_onset,
                          duration=duration,
                          DisabilityWeight=DisabilityWeight)
  #print(str(data_plugin))
  require(parallel)
  require(triangle)
  source(paste0(filepath,'/YLD_related_function','/YLD_subfunction.R'))
  
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
    #up_bound<-unlist(mclapply(2:length(age_labels),function(i) age_ForTurningonset[1:length(age_labels)][i]-0.1))
    
    #uniform_prior_range_age_at_onset<-data.frame(lr=age_ForTurningonset[1:length(age_labels)],
    #                                up=c(up_bound,max_age))
    #round(rowMeans(uniform_prior_range),digits = 1)
    
    # to comply with uncertainty analysis of YLL, uncertainty for population variable is also considered
    
    source(paste0(filepath,'/YLD_related_function','/YLD_priors.R'))
    
    YLD_priors_samples<-list()
    for(i in 1:length(age_labels)){
      sampling_outcome<-matrix(c(replicate(nTrials,expr = population_priors(x=population[i],scale = prior_population)),
                                 replicate(nTrials,expr = death_priors(incidence_forcal[i])),
                                 replicate(nTrials,expr = Durantion_priors(Duration_interval,i)),
                                 replicate(nTrials,expr = DisabilityWeight_priors(DisabilityWeight_interval,i))),nrow=nTrials,byrow=FALSE)
      dimnames(sampling_outcome)[[2]]<-c('Population','Incidence','Duration','Disability_weight')
      YLD_priors_samples[[i]]<-sampling_outcome
    }
    
    
    
    uncertainty_results<-mclapply(1:length(YLD_priors_samples), function(i){
      tempdata<-apply(YLD_priors_samples[[i]],2,quantile,probs=c(uncertainty_alpha/2,1-uncertainty_alpha/2))
      
      tempdata<-as.matrix(tempdata)
      
    })
    
    ## population
    Population_uncertainty<-data.frame(age_labels,population,
                                       round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,1]})),digits = 0))
    names(Population_uncertainty)<-c('Age_groups','Recorded_population','Lr','Up')
    
    ## incidence
    Counts_uncertainty<-data.frame(age_labels,incidence_forcal,
                                   round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,2]})),digits = 0))
    names(Counts_uncertainty)<-c('Age_groups','Recorded_incidence','Lr','Up')
    
    ## duration
    Dur_uncertainty<-data.frame(age_labels,duration,
                                round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,3]})),digits = 1))
    names(Dur_uncertainty)<-c('Age_groups','Recorded_Duration','Lr','Up')
    
    ## DisabilityWeight
    DW_uncertainty<-data.frame(age_labels,DisabilityWeight,
                               round(t(sapply(1:length(uncertainty_results),FUN = function(i){uncertainty_results[[i]][,4]})),digits = 1))
    names(DW_uncertainty)<-c('Age_groups','Recorded_DW','Lr','Up')
    
    YLD_point_est<-YLD_sub_function(data_plugin = data.frame(Age=age_labels,
                                                             population=Population_uncertainty$Recorded_population,
                                                             incidence=Counts_uncertainty$Recorded_incidence,
                                                             age_at_onset=age_at_onset,
                                                             duration=Dur_uncertainty$Recorded_Duration,
                                                             DisabilityWeight=DW_uncertainty$Recorded_DW),
                                    Rate=Rate,
                                    Beta=Beta,
                                    Const=Const,
                                    Agewt=Agewt,
                                    BplusR=BplusR)
    ## lr
    YLD_lr<-YLD_sub_function(data_plugin = data.frame(Age=age_labels,
                                                      population=Population_uncertainty$Lr,
                                                      incidence=Counts_uncertainty$Lr,
                                                      age_at_onset=age_at_onset,
                                                      duration=Dur_uncertainty$Lr,
                                                      DisabilityWeight=DW_uncertainty$Lr),
                             Rate=Rate,
                             Beta=Beta,
                             Const=Const,
                             Agewt=Agewt,
                             BplusR=BplusR)
    
    ## Up
    YLD_up<-YLD_sub_function(data_plugin = data.frame(Age=age_labels,
                                                      population=Population_uncertainty$Up,
                                                      incidence=Counts_uncertainty$Up,
                                                      age_at_onset=age_at_onset,
                                                      duration=Dur_uncertainty$Up,
                                                      DisabilityWeight=DW_uncertainty$Up),
                             Rate=Rate,
                             Beta=Beta,
                             Const=Const,
                             Agewt=Agewt,
                             BplusR=BplusR)
    
    ### point est
    population_withTotal<-c(population,sum(population))
    incidence_forcal_total<-c(incidence_forcal,sum(incidence_forcal))
    age_at_onset_wihtTotal<-c(age_at_onset,sum(age_at_onset*incidence_forcal)/sum(incidence_forcal))
    duration_withTotal<-c(duration,sum(duration*incidence_forcal)/sum(incidence_forcal))
    DisabilityWeight_withTotal<-c(DisabilityWeight,sum(DisabilityWeight*incidence_forcal)/sum(incidence_forcal))
    
    YLD_withTotal<-c(YLD_point_est,sum(YLD_point_est))
    YLD_per<-YLD_withTotal/population_withTotal*YLD_perUnit
    
    
    ### lr
    population_withTotal_lr<-c(Population_uncertainty$Lr,sum(Population_uncertainty$Lr))
    incidence_forcal_total_lr<-c(Counts_uncertainty$Lr,sum(Counts_uncertainty$Lr))
    age_at_onset_wihtTotal_lr<-c(age_at_onset,sum(age_at_onset*Counts_uncertainty$Lr)/sum(Counts_uncertainty$Lr))
    duration_withTotal_lr<-c(Dur_uncertainty$Lr,sum(Dur_uncertainty$Lr*Counts_uncertainty$Lr)/sum(Counts_uncertainty$Lr))
    DisabilityWeight_withTotal_lr<-c(DW_uncertainty$Lr,sum(DW_uncertainty$Lr*Counts_uncertainty$Lr)/sum(Counts_uncertainty$Lr))
    
    YLD_withTotal_lr<-c(YLD_lr,sum(YLD_lr))
    YLD_per_lr<-YLD_withTotal_lr/population_withTotal_lr*YLD_perUnit
    
    ### Up
    population_withTotal_up<-c(Population_uncertainty$Up,sum(Population_uncertainty$Up))
    incidence_forcal_total_up<-c(Counts_uncertainty$Up,sum(Counts_uncertainty$Up))
    age_at_onset_wihtTotal_up<-c(age_at_onset,sum(age_at_onset*Counts_uncertainty$Up)/sum(Counts_uncertainty$Up))
    duration_withTotal_up<-c(Dur_uncertainty$Up,sum(Dur_uncertainty$Up*Counts_uncertainty$Up)/sum(Counts_uncertainty$Up))
    DisabilityWeight_withTotal_up<-c(DW_uncertainty$Up,sum(DW_uncertainty$Up*Counts_uncertainty$Up)/sum(Counts_uncertainty$Up))
    
    YLD_withTotal_up<-c(YLD_up,sum(YLD_up))
    YLD_per_up<-YLD_withTotal_up/population_withTotal_up*YLD_perUnit
    
    YLD_DataFrame<-data.frame(Age=c(as.character(age_labels),'Total'),
                              population=population_withTotal,
                              incidence_case=incidence_forcal_total,
                              age_at_onset=age_at_onset_wihtTotal,
                              duration=duration_withTotal,
                              DisabilityWeight=DisabilityWeight_withTotal,
                              YLD=YLD_withTotal,
                              YLD_per=YLD_per)
    
    YLD_DataFrame_lr<-data.frame(Age=c(as.character(age_labels),'Total'),
                                 population=population_withTotal_lr,
                                 incidence_case=incidence_forcal_total_lr,
                                 age_at_onset=age_at_onset_wihtTotal_lr,
                                 duration=duration_withTotal_lr,
                                 DisabilityWeight=DisabilityWeight_withTotal_lr,
                                 YLD=YLD_withTotal_lr,
                                 YLD_per=YLD_per_lr)
    YLD_DataFrame_up<-data.frame(Age=c(as.character(age_labels),'Total'),
                                 population=population_withTotal_up,
                                 incidence_case=incidence_forcal_total_up,
                                 age_at_onset=age_at_onset_wihtTotal_up,
                                 duration=duration_withTotal_up,
                                 DisabilityWeight=DisabilityWeight_withTotal_up,
                                 YLD=YLD_withTotal_up,
                                 YLD_per=YLD_per_up)
    names(YLD_DataFrame)[ncol(YLD_DataFrame)]<-paste0('YLD','_per',YLD_perUnit)
    names(YLD_DataFrame_lr)[ncol(YLD_DataFrame_lr)]<-paste0('YLD','_lr','_per',YLD_perUnit)
    names(YLD_DataFrame_up)[ncol(YLD_DataFrame_up)]<-paste0('YLD','_up','_per',YLD_perUnit)
    
    
    Original_data<-data.frame(Age=c(as.character(age_labels),'Total'),
                              population=population_withTotal,
                              incidence=incidence_forcal_total,
                              incidence_rate=incidence_forcal_total/population_withTotal*RateUNIT,
                              age_at_onset=age_at_onset_wihtTotal,
                              duration=duration_withTotal,
                              DisabilityWeight=DisabilityWeight_withTotal)
    names(Original_data)[4]<-paste0( "incidence_rate" ,'_per',RateUNIT)
    
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
                     Original_data=Original_data[-nrow(Original_data),],
                     YLD_DataFrame=YLD_DataFrame,
                     YLD_DataFrame_lr=YLD_DataFrame_lr,
                     YLD_DataFrame_up=YLD_DataFrame_up,
                     UncertaintyFile=list(population=Population_uncertainty,
                                          incidence_case=Counts_uncertainty,
                                          duration=Dur_uncertainty,
                                          disability_weight=DW_uncertainty),
                     input_interval=list(Duration_interval=Duration_interval,
                                         DisabilityWeight_interval=DisabilityWeight_interval))
    
    } else if(uncertainty_range==FALSE){
      
      YLD<-YLD_sub_function(data_plugin = data_sort,
                            Rate=Rate,
                            Beta=Beta,
                            Const=Const,
                            Agewt=Agewt,
                            BplusR=BplusR)
      
      population_withTotal<-c(population,sum(population))
      incidence_forcal_total<-c(incidence_forcal,sum(incidence_forcal))
      age_at_onset_wihtTotal<-c(age_at_onset,sum(age_at_onset*incidence_forcal)/sum(incidence_forcal))
      duration_withTotal<-c(duration,sum(duration*incidence_forcal)/sum(incidence_forcal))
      DisabilityWeight_withTotal<-c(DisabilityWeight,sum(DisabilityWeight*incidence_forcal)/sum(incidence_forcal))
      
      YLD_withTotal<-c(YLD,sum(YLD))
      YLD_per<-YLD_withTotal/population_withTotal*YLD_perUnit
      
      YLD_DataFrame<-data.frame(Age=c(as.character(age_labels),'Total'),
                                population=population_withTotal,
                                incidence=incidence_forcal_total,
                                age_at_onset=age_at_onset_wihtTotal,
                                duration=duration_withTotal,
                                DisabilityWeight=DisabilityWeight_withTotal,
                                YLD=YLD_withTotal,
                                YLD_per=YLD_per)
      
      names(YLD_DataFrame)[ncol(YLD_DataFrame)]<-paste0( "YLD" ,'_per',YLD_perUnit)
      Original_data<-data.frame(Age=c(as.character(age_labels),'Total'),
                                population=population_withTotal,
                                incidence=incidence_forcal_total,
                                incidence_rate=incidence_forcal_total/population_withTotal*RateUNIT,
                                age_at_onset=age_at_onset_wihtTotal,
                                duration=duration_withTotal,
                                DisabilityWeight=DisabilityWeight_withTotal)
      names(Original_data)[4]<-paste0( "incidence_rate" ,'_per',RateUNIT)
      
      YLD_object<-list(Age_group=c(as.character(age_labels),'Total'),
                       YLD=YLD_withTotal,
                       YLD_per=YLD_per,
                       #YLD_uncertainty=list(uncertainty_range,uncertainty_alpha),
                       Original_data=Original_data,
                       YLD_DataFrame=YLD_DataFrame)
      
    }
  
  class(YLD_object) <- "Incident_YLL_object"
  attr(YLD_object,"Call") <- sys.call()
  return(YLD_object)
  
  }