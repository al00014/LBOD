prev_to_prevYLD<-function(prevalent_age_labels,
                        prevalent_data,
                        input_population, # population should be annual, and should have the same annual data as incident data
                        input_DisabilityWeight,
                        input_DisabilityWeight_interval,
                        year_range,
                        Output_YLD_noUI=TRUE,# if only output YLD, no uncertainty is analyzed, no further detail information will be output,#
                        #  it should be a quicker way to get a teaser of the output YLD
                        RateUNIT=1000,
                        Input_Prevalent_Rate=FALSE,
                        uncertainty_alpha=0.05,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                        prior_population=0.001,
                        nTrials=2000,
                        YLD_perUnit=1000,
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
    print('Prevalent age groups are as follows,')
    print(prevalent_age_labels)
    print('-----------------------------')
    
  }
  if(!(file.exists(paste0(filepath,'/YLD_related_function','/YLD_prevalent.R')))){
    stop('\n
         Cannot locate YLD_prevalent.R files,
         Please check the integrity of source files.')
  }
  
  if(ncol(input_population)>=ncol(prevalent_data)){
    population_data<-input_population[,1:ncol(prevalent_data)]
  } else if(ncol(input_population)<ncol(prevalent_data)){
    stop('\n
         Not enough years of data for input_population, cannot calculate YLD,\n
         Please ensure that the input_population data and the prevalent data \n
         has the same annual data')
  }
  
  
  source(paste0(filepath,'/YLD_related_function','/YLD_prevalent.R')) # YLD_nodefaults.R
  incident_cols<-length(dimnames(prevalent_data)[[2]])
  
  

  #age_at_onset_final<-age_at_onset ###
  
  #prevalent_age_labels<-factor(prevalent_age_labels)
  
  
  #require(parallel)
  
  
  #population_data<-input_population[,1:incident_cols] ###
  
  if(Output_YLD_noUI==TRUE){
    YLD_byyear<-matrix(unlist(parallel::mclapply(1:incident_cols,FUN=function(k){
      #print(k)
      tem<-YLD_prevalent(age_labels=prevalent_age_labels,
                        population=population_data[,k],
                        prevalence=prevalent_data[,k],
                        DisabilityWeight=input_DisabilityWeight,
                        DisabilityWeight_interval=input_DisabilityWeight_interval,
                        prior_population=prior_population,
                        nTrials=nTrials,
                        uncertainty_range=FALSE,
                        Input_Prevalent_Rate=Input_Prevalent_Rate,
                        uncertainty_alpha=uncertainty_alpha,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                        max_age=max_age,
                        RateUNIT=RateUNIT,
                        YLD_perUnit=YLD_perUnit,
                        filepath=filepath
      )
      return(tem$YLD)
    }
    )),ncol=incident_cols)
    
    
    YLD_byyear<-data.frame(YLD_byyear)
    names(YLD_byyear)<-paste0('Y',year_range)
    row.names(YLD_byyear)<-c(prevalent_age_labels,'Total')
    
    YLD_transformed<-list(YLD_byyear=YLD_byyear)
    
  } else {
    
    YLD_object_by_year<-list()
    for(i in 1:incident_cols){
      tem<-YLD_prevalent(age_labels=prevalent_age_labels,
                        population=population_data[,i],
                        prevalence=prevalent_data[,i],
                        DisabilityWeight=input_DisabilityWeight,
                        DisabilityWeight_interval=input_DisabilityWeight_interval,
                        prior_population=prior_population,
                        nTrials=nTrials,
                        uncertainty_range=TRUE,
                        Input_Prevalent_Rate=Input_Prevalent_Rate,
                        uncertainty_alpha=uncertainty_alpha,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                        max_age=max_age,
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
    row.names(YLD_byyear)<-c(prevalent_age_labels,'Total')
    
    names(YLD_lr_byyear)<-paste0('Y',year_range)
    row.names(YLD_lr_byyear)<-c(prevalent_age_labels,'Total')
    
    names(YLD_up_byyear)<-paste0('Y',year_range)
    row.names(YLD_up_byyear)<-c(prevalent_age_labels,'Total')
    
    names(YLDper_byyear)<-paste0('Y',year_range)
    row.names(YLDper_byyear)<-c(prevalent_age_labels,'Total')
    
    names(YLDper_lr_byyear)<-paste0('Y',year_range)
    row.names(YLDper_lr_byyear)<-c(prevalent_age_labels,'Total')
    
    names(YLDper_up_byyear)<-paste0('Y',year_range)
    row.names(YLDper_up_byyear)<-c(prevalent_age_labels,'Total')
    
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
