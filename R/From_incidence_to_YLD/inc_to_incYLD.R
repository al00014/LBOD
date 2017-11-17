inc_to_incYLD<-function(incident_age_labels,
                        incident_data,
                        input_age_at_onset,
                        input_population, # population should be annual, and should have the same annual data as incident data
                        input_duration,  
                        input_DisabilityWeight,
                        input_Duration_interval,
                        input_DisabilityWeight_interval,
                        year_range,
                        Output_YLD_noUI=TRUE,# if only output YLD, no uncertainty is analyzed, no further detail information will be output,#
                                             #  it should be a quicker way to get a teaser of the output YLD
                        RateUNIT=1000,
                        Input_Incident_Rate=FALSE,
                        uncertainty_alpha=0.05,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                        prior_population=0.001,
                        nTrials=2000,
                        YLD_perUnit=1000,
                        Rate=0.03,
                        Beta=0.04,
                        Const=0.1658,
                        Agewt=0,
                        max_age=90,
                        filepath=paste0(getwd(),'/R'),
                        verbose=TRUE,...){
  #incidence_data=incidence_data
  #population=population
  options(scipen = 999)
  #if(sum(!c('parallel') %in% rownames(installed.packages()) )!= 0) {
  #  stop('\n
  #       Please installed the dependencies first:\n
  #       "parallel" ')
  #}
  
  if(verbose==TRUE){
    print('Incidence age groups are as follows,')
    print(incident_age_labels)
    print('-----------------------------')
    
  }
  if(!(file.exists(paste0(filepath,'/YLD_related_function','/YLD.R')))){
    stop('\n
         Cannot locate YLD.R files,
         Please check the integrity of source files.')
  }
  
  source(paste0(filepath,'/YLD_related_function','/YLD.R')) # YLD_nodefaults.R
  incident_cols<-length(dimnames(incident_data)[[2]])
  
  
            if(class(input_age_at_onset)=='numeric'){
              message('\n
                      Values of input_age_at_onset will be considered uniform across all the years (year by columns),\n
                      as there is only one vector of data for input_age_at_onset')
              
              age_at_onset_final<-matrix(rep(input_age_at_onset,times=as.numeric(ncol(incident_data))),ncol=ncol(incident_data))
              
            } else if(class(input_age_at_onset)=='matrix' |class(input_age_at_onset)=='data.frame'){
              print(paste0('There are a total of ',ncol(input_age_at_onset),
                           ' columns in the input_age_at_onset argument'))
              print('-------------------')
              message('It is a matrix or data.frame, each column of input_age_at_onset will be used separately in calculation of YLD')
              
              if(ncol(incident_data)!=ncol(input_age_at_onset)){
                stop('\n 
                     The number of columns of input_age_at_onset are not the same as those in incidence data,\n
                     Check the input_age_at_onset argument, and make sure its columns are the same as incidence data,\n
                     If several years of input_age_at_onset are considered uniform, make replication of that/those year(s) of columns\n
                     to ensure conformity')
              }
              age_at_onset_final<-input_age_at_onset
            }
            
            if(ncol(input_population)>=ncol(incident_data)){
              population_data<-input_population[,1:ncol(incident_data)]
            } else if(ncol(input_population)<ncol(incident_data)){
              stop('\n
                   Not enough years of data for input_population, cannot calculate YLD,\n
                   Please ensure that the input_population data and the incident data \n
                   has the same annual data')
            }
  
  
  #age_at_onset_final<-age_at_onset ###
  
  #incident_age_labels<-factor(incident_age_labels)
  

  #require(parallel)
  
  
  #population_data<-input_population[,1:incident_cols] ###
  
  if(Output_YLD_noUI==TRUE){
    YLD_byyear<-matrix(unlist(parallel::mclapply(1:incident_cols,FUN=function(k){
      #print(k)
      tem<-YLD_incident(age_labels=incident_age_labels,
                        population=population_data[,k],
                        incidence=incident_data[,k],
                        age_at_onset=age_at_onset_final[,k],
                        duration=input_duration,
                        DisabilityWeight=input_DisabilityWeight,
                        Duration_interval=input_Duration_interval,
                        DisabilityWeight_interval=input_DisabilityWeight_interval,
                        prior_population=prior_population,
                        nTrials=nTrials,
                        uncertainty_range=FALSE,
                        Input_Incident_Rate=Input_Incident_Rate,
                        uncertainty_alpha=uncertainty_alpha,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                        max_age=max_age,
                        Rate=Rate,
                        Beta=Beta,
                        Const=Const,
                        Agewt=Agewt,
                        RateUNIT=RateUNIT,
                        YLD_perUnit=YLD_perUnit,
                        filepath=filepath
      )
      return(tem$YLD)
    }
    )),ncol=incident_cols)
    
    
    YLD_byyear<-data.frame(YLD_byyear)
    names(YLD_byyear)<-paste0('Y',year_range)
    row.names(YLD_byyear)<-c(incident_age_labels,'Total')
    
    YLD_transformed<-list(YLD_byyear=YLD_byyear)
    
  } else {
    
    YLD_object_by_year<-list()
    for(i in 1:incident_cols){
      tem<-YLD_incident(age_labels=incident_age_labels,
                               population=population_data[,i],
                               incidence=incident_data[,i],
                               age_at_onset=age_at_onset_final[,i],
                               duration=input_duration,
                               DisabilityWeight=input_DisabilityWeight,
                               Duration_interval=input_Duration_interval,
                               DisabilityWeight_interval=input_DisabilityWeight_interval,
                               prior_population=prior_population,
                               nTrials=nTrials,
                               uncertainty_range=TRUE,
                               Input_Incident_Rate=Input_Incident_Rate,
                               uncertainty_alpha=uncertainty_alpha,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                               max_age=max_age,
                               Rate=Rate,
                               Beta=Beta,
                               Const=Const,
                               Agewt=Agewt,
                               RateUNIT=RateUNIT,
                               YLD_perUnit=YLD_perUnit,
                               filepath=filepath
      )
      
      YLD_object_by_year[[i]]<-tem
    }
    
    names(YLD_object_by_year)<-paste0('Y',year_range)
    
    YLD_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(YLD_object_by_year),FUN=function(i){
      YLD_object_by_year[[i]]$YLD
    })),ncol=incident_cols))
    
    YLD_lr_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(YLD_object_by_year),FUN=function(i){
      YLD_object_by_year[[i]]$YLD_lr
    })),ncol=incident_cols))
    
    YLD_up_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(YLD_object_by_year),FUN=function(i){
      YLD_object_by_year[[i]]$YLD_up
    })),ncol=incident_cols))
    
    
    ### personunits
    YLDper_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(YLD_object_by_year),FUN=function(i){
      YLD_object_by_year[[i]]$YLD_per
    })),ncol=incident_cols))
    
    YLDper_lr_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(YLD_object_by_year),FUN=function(i){
      YLD_object_by_year[[i]]$YLD_per_lr
    })),ncol=incident_cols))
    
    YLDper_up_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(YLD_object_by_year),FUN=function(i){
      YLD_object_by_year[[i]]$YLD_per_up
    })),ncol=incident_cols))
    
    
    
    ### defined names and cols
    names(YLD_byyear)<-paste0('Y',year_range)
    row.names(YLD_byyear)<-c(incident_age_labels,'Total')
    
    names(YLD_lr_byyear)<-paste0('Y',year_range)
    row.names(YLD_lr_byyear)<-c(incident_age_labels,'Total')
    
    names(YLD_up_byyear)<-paste0('Y',year_range)
    row.names(YLD_up_byyear)<-c(incident_age_labels,'Total')
    
    names(YLDper_byyear)<-paste0('Y',year_range)
    row.names(YLDper_byyear)<-c(incident_age_labels,'Total')
    
    names(YLDper_lr_byyear)<-paste0('Y',year_range)
    row.names(YLDper_lr_byyear)<-c(incident_age_labels,'Total')
    
    names(YLDper_up_byyear)<-paste0('Y',year_range)
    row.names(YLDper_up_byyear)<-c(incident_age_labels,'Total')
    
	print('========================================') 
	message(paste0('YLD rates were calculated in the unit of ',YLD_perUnit,'.'))
	print('========================================') 
	
    YLD_transformed<-list(YLD_byyear=YLD_byyear,
              YLD_lr_byyear=YLD_lr_byyear,
              YLD_up_byyear=YLD_up_byyear,
              YLDper_byyear=YLDper_byyear,
              YLDper_lr_byyear=YLDper_lr_byyear,
              YLDper_up_byyear=YLDper_up_byyear,
              YLD_object_by_year=YLD_object_by_year)
    names(YLD_transformed)[4:6]<-paste0(rep('YLDper',times=3),#YLD_perUnit,
										c('_byyear','_lr_byyear','_up_byyear'))
    
    
    
  }
  class(YLD_transformed) <- "YLD_object"
  attr(YLD_transformed,"Call") <- sys.call()
  return(YLD_transformed)
  
}
