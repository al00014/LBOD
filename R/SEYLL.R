SEYLL<-function(SLE_object,
                personunits=1000,
                Rate=0.03,
                Beta=0.04,
                Const=0.1658,
                Agewt=0#,
                #filepath=paste0('./BOD/R')
				){
  # Discount rate (r)	Standard discount rate is 0.03	
  # Beta (b)		Standard age weights use beta=0.04
  # Constant (C)		Standard age weights use C=0.1658
  # Agewt(K)		K=0 (no age weights) to 1 (full age weights)
  options(scipen = 999)
  #if(sum(!c('parallel') %in% rownames(installed.packages()) )!= 0) {
  #  stop('\n
  #        Please installed the dependencies first:\n
  #       "parallel" ')
  #}
  # require(parallel)
  # filepath,'/SEYLL_related_function',
  if(file.exists(paste0('./YLL_subfunction.R'))==FALSE){
    stop('\n
         Missing YLL_subfunction.R file, cannot estimate YLL,\n
         Please check the integrity of source files')
  }
  source(paste0('./YLL_subfunction.R'))
  # filepath,'/SEYLL_related_function',
  BplusR= -(Beta+Rate)
  Stand_life_expectancy<-SLE_object$SLE
  data_plugin<-data.frame(Age=SLE_object$Original_data$Age_groups,
                          population=SLE_object$Original_data$population,
                          death=SLE_object$Original_data$death,
                          age_av_death=SLE_object$Original_data$age_av_death)
  
  if(SLE_object$uncertainty_output==TRUE){
    
    YLL_results_partly<-YLL_sub_function(data_plugin = data_plugin,
                                         Stand_life_expectancy = Stand_life_expectancy,
                                         Rate=Rate,
                                         Beta=Beta,
                                         Const=Const,
                                         Agewt=Agewt,
                                         BplusR=BplusR)
    population_withTotal<-c(data_plugin$population,
                            sum(data_plugin$population,na.rm = TRUE))
    
    death_withTotal<-c(data_plugin$death,
                       sum(data_plugin$death,na.rm = TRUE))
    
    Age_av_withTotal<-c(data_plugin$age_av_death,
                        sum(data_plugin$death*data_plugin$age_av_death,na.rm = TRUE)/death_withTotal[length(death_withTotal)])
    YLL_withTotal<-c(YLL_results_partly,sum(YLL_results_partly,na.rm = TRUE))
    
    YLL_per_withTotal<-YLL_withTotal/population_withTotal*personunits
    
    
    YLL_DataFrame=data.frame(Age_group=as.factor(c(as.character(SLE_object$Original_data$Age_groups),c('Total'))),
                             population=population_withTotal,
                             death=death_withTotal,
                             age_av_death=Age_av_withTotal,
                             YLL=YLL_withTotal,
                             YLL_perUnit=YLL_per_withTotal)
    names(YLL_DataFrame)[ncol(YLL_DataFrame)]<-paste0('YLL_per',personunits)
    
    ### lower
    data_plugin_lr<-data.frame(Age=SLE_object$Original_data$Age_groups,
                               population=SLE_object$Population_uncertainty$Lr,
                               death=SLE_object$Death_uncertainty$Lr,
                               age_av_death=SLE_object$Average_at_death_uncertainty$Lr)
    YLL_DataFrame_lr=data.frame(Age_group=as.factor(c(as.character(SLE_object$Original_data$Age_groups),c('Total'))),
                                population_lr=c(SLE_object$Population_uncertainty$Lr,
                                                sum(SLE_object$Population_uncertainty$Lr,na.rm = TRUE)),
                                death_lr=c(SLE_object$Death_uncertainty$Lr,
                                           sum(SLE_object$Death_uncertainty$Lr,na.rm = TRUE)),
                                
                                age_av_death_lr=c(SLE_object$Average_at_death_uncertainty$Lr,
                                                  sum(SLE_object$Death_uncertainty$Lr*SLE_object$Average_at_death_uncertainty$Lr,na.rm = TRUE)/sum(SLE_object$Death_uncertainty$Lr,na.rm = TRUE)),
                                YLL_lr=c(YLL_sub_function(data_plugin = data_plugin_lr,
                                                          Stand_life_expectancy = SLE_object$SLE_lr,
                                                          Rate=Rate,
                                                          Beta=Beta,
                                                          Const=Const,
                                                          Agewt=Agewt,
                                                          BplusR=BplusR),sum(YLL_sub_function(data_plugin = data_plugin_lr,
                                                                                              Stand_life_expectancy = SLE_object$SLE_lr,
                                                                                              Rate=Rate,
                                                                                              Beta=Beta,
                                                                                              Const=Const,
                                                                                              Agewt=Agewt,
                                                                                              BplusR=BplusR),na.rm = TRUE)))
    YLL_DataFrame_lr$YLL_lr_per<-YLL_DataFrame_lr$YLL_lr/YLL_DataFrame_lr$population_lr*personunits
    names(YLL_DataFrame_lr)[ncol(YLL_DataFrame_lr)]<-paste0('YLL','_lr','_per',personunits)
    Stand_life_expectancy
    ### uper
    data_plugin_up<-data.frame(Age=SLE_object$Original_data$Age_groups,
                               population=SLE_object$Population_uncertainty$Up,
                               death=SLE_object$Death_uncertainty$Up,
                               age_av_death=SLE_object$Average_at_death_uncertainty$Up)
    YLL_DataFrame_up<-data.frame(Age_group=as.factor(c(as.character(SLE_object$Original_data$Age_groups),c('Total'))),
                                 population_up=c(SLE_object$Population_uncertainty$Up,
                                                 sum(SLE_object$Population_uncertainty$Up,na.rm = TRUE)),
                                 death_up=c(SLE_object$Death_uncertainty$Up,
                                            sum(SLE_object$Death_uncertainty$Up,na.rm = TRUE)),
                                 
                                 age_av_death_up=c(SLE_object$Average_at_death_uncertainty$Up,
                                                   sum(SLE_object$Death_uncertainty$Up*SLE_object$Average_at_death_uncertainty$Up,na.rm = TRUE)/sum(SLE_object$Death_uncertainty$Up,na.rm = TRUE)),
                                 YLL_up=c(YLL_sub_function(data_plugin = data_plugin_up,
                                                           Stand_life_expectancy = SLE_object$SLE_up,
                                                           Rate=Rate,
                                                           Beta=Beta,
                                                           Const=Const,
                                                           Agewt=Agewt,
                                                           BplusR=BplusR),sum(YLL_sub_function(data_plugin = data_plugin_up,
                                                                                               Stand_life_expectancy = SLE_object$SLE_up,
                                                                                               Rate=Rate,
                                                                                               Beta=Beta,
                                                                                               Const=Const,
                                                                                               Agewt=Agewt,
                                                                                               BplusR=BplusR),na.rm = TRUE)))
    YLL_DataFrame_up$YLL_up_per<-YLL_DataFrame_up$YLL_up/YLL_DataFrame_up$population_up*personunits
    names(YLL_DataFrame_up)[ncol(YLL_DataFrame_up)]<-paste0('YLL','_up','_per',personunits)
    #Original_data
    YLL_object<-list(Age_group=as.factor(c(as.character(SLE_object$Original_data$Age_groups),c('Total'))),
                     YLL=YLL_withTotal,
                     YLL_lr=YLL_DataFrame_lr$YLL_lr,
                     YLL_up=YLL_DataFrame_up$YLL_up,
                     YLL_perUnit=YLL_per_withTotal,
                     YLL_perUnit_lr=YLL_DataFrame_lr$YLL_lr_per1000,
                     YLL_perUnit_up=YLL_DataFrame_up$YLL_up_per1000,
                     perUnit=personunits,
                     YLL_uncertainty=list(Calculate_range=SLE_object$uncertainty_output,
                                      Alpha=SLE_object$alpha),
                     Original_data=SLE_object$Original_data, # for the original data set from YLL_object, no total sum is retained as it is designed conveniently for convergence_test function
                                  # rbind(SLE_object$Original_data,
                                  #       c(factor('Total'),colSums(SLE_object$Original_data[,2:3]),
                                  #                                  sum(SLE_object$Original_data$death*SLE_object$Original_data$age_av_death)/sum(SLE_object$Original_data$death))),
                     Original_life_table=SLE_object$Original_life_table,
                     YLL_DataFrame=YLL_DataFrame,
                     YLL_DataFrame_lr=YLL_DataFrame_lr,
                     YLL_DataFrame_up=YLL_DataFrame_up,
                     Original_uncertaintyFile=list(population=SLE_object$Population_uncertainty,
                                                   death_counts=SLE_object$Death_uncertainty,
                                                   age_av_death=SLE_object$Average_at_death_uncertainty))
    
  } else if(SLE_object$uncertainty_output==FALSE){
    
    YLL_results_partly<-YLL_sub_function(data_plugin = data_plugin,
                                         Stand_life_expectancy = Stand_life_expectancy,
                                         Rate=Rate,
                                         Beta=Beta,
                                         Const=Const,
                                         Agewt=Agewt,
                                         BplusR=BplusR)
    population_withTotal<-c(data_plugin$population,
                            sum(data_plugin$population,na.rm = TRUE))
    
    death_withTotal<-c(data_plugin$death,
                       sum(data_plugin$death,na.rm = TRUE))
    
    Age_av_withTotal<-c(data_plugin$age_av_death,
                        sum(data_plugin$death*data_plugin$age_av_death,na.rm = TRUE)/death_withTotal[length(death_withTotal)])
    YLL_withTotal<-c(YLL_results_partly,sum(YLL_results_partly,na.rm = TRUE))
    
    YLL_per_withTotal<-YLL_withTotal/population_withTotal*personunits
    
    YLL_DataFrame=data.frame(Age_group=as.factor(c(as.character(SLE_object$Original_data$Age_groups),c('Total'))),
                             population=population_withTotal,
                             death=death_withTotal,
                             age_av_death=Age_av_withTotal,
                             YLL=YLL_withTotal,
                             YLL_perUnit=YLL_per_withTotal)
    names(YLL_DataFrame)[ncol(YLL_DataFrame)]<-paste0('YLL_per',personunits)
    YLL_object<-list(Age_group=as.factor(c(as.character(SLE_object$Original_data$Age_groups),c('Total'))),
                     YLL=YLL_withTotal,
                     YLL_perUnit=YLL_per_withTotal,
                     perUnit=personunits,
                     YLL_uncertainty=SLE_object$uncertainty_output,
                     Original_data=SLE_object$Original_data,
                     Original_life_table=SLE_object$Original_life_table,
                     YLL_DataFrame=YLL_DataFrame)
  }
  
  
  class(YLL_object) <- "SEYLL_object"
  attr(YLL_object,"Call") <- sys.call()
  return(YLL_object)
  }