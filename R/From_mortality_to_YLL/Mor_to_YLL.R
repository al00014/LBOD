Mor_to_YLL<-function(mortality_age_labels,
                                       mortality_data,
                                       year_range,
                                       age_average_at_death,
                                       population, # population should be annual, and should have the same annual data as mortality data
                                       standard_LE,  
                                       standard_LT_age, # standard life table recommended by the WHO, currently refers to the Coale_and_Demeny_level_west_25_26 life-table of Japen
                                       uncertainty_alpha=0.05,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                                       prior_population=0.001,
                                       nTrials=2000,
                                       Output_YLL_noUI=TRUE,# if only output YLL, no uncertainty is analyzed, no further detail information will be output,#
                                                            #  it should be a quicker way to get a teaser of the output YLL
                                       personunits=1000,
                                       Rate=0.03,
                                       Beta=0.04,
                                       Const=0.1658,
                                       Agewt=0,
                                       max_age=100,
                                       filepath=paste0('./BOD/R'),
                                       verbose=TRUE){
  options(scipen = 999)
  if(verbose==TRUE){
    print('Mortality age groups are as follows,')
    print(mortality_age_labels)
    print('-----------------------------')
    print('Standard life table age groups are as follows:')
    print(standard_LT_age)
    print('-----------------------------')
    message('\n
            Before further calculation, please check the age groups of standard life table, and make changes where appropriate, \n
            Normally, take average of Life expectancy of several age groups should work, \n
            and discard a few age groups, make the age interval wider in other words,')
  }
  if(class(age_average_at_death)=='numeric'){
  
	print('-----------------------------')
    message('\n
            Values of age_average_at_death will be considered uniform across all the years (year by columns),\n
            as there is only one vector of data for age_average_at_death')
    age_average_at_death_final<-matrix(rep(age_average_at_death,times=ncol(mortality_data)),ncol=ncol(mortality_data))
  } else if(class(age_average_at_death)=='matrix' |class(age_average_at_death)=='data.frame'){
	print('-----------------------------')
    print(paste0('There are a total of ',ncol(age_average_at_death),' columns in the age_average_at_death argument'))
    message('It is a matrix or data.frame, each column of age_average_at_death will be used separately in calculation of YLL')
    if(ncol(mortality_data)!=ncol(age_average_at_death)){
      stop('\n 
           The number of columns of age_average_at_death are not the same as those in mortality data,\n
           Check the age_average_at_death argument, and make sure its columns are the same as mortality data,\n
           If several years of age_average_at_death are considered uniform, make replication of that/those year(s) of columns\n
           to ensure conformity')
    }
    age_average_at_death_final<-age_average_at_death
    }
  
  if(!(file.exists(paste0(filepath,'/SLE_related_function','/SLE.R')))){
    stop('\n
         Cannot locate SLE.R files,
         Please check the integrity of source files.')
  }
  if(!(file.exists(paste0(filepath,'/SEYLL_related_function','/SEYLL.R')))){
    stop('\n
         Cannot locate SEYLL.R files,
         Please check the integrity of source files.')
  }
  if(ncol(population)>=ncol(mortality_data)){
    population_data<-population[,1:ncol(mortality_data)]
  } else if(ncol(population)<ncol(mortality_data)){
    stop('\n
         Not enough years of data for population, cannot calculate YLL,\n
         Please ensure that the population data and the mortality data \n
         has the same annual data')
  }
  
  source(paste0(filepath,'/SLE_related_function','/SLE.R'))
  source(paste0(filepath,'/SEYLL_related_function','/SEYLL.R'))
  #require(parallel)
  
  if(Output_YLL_noUI==TRUE){
    YLL_byyear<-matrix(unlist(parallel::mclapply(1:ncol(mortality_data),FUN=function(i){
      SLE_object<-SLE(age_labels=mortality_age_labels,
                      population=population_data[,i],
                      death_counts=mortality_data[,i],
                      age_average_at_death=age_average_at_death_final[,i],
                      standard_LE=standard_LE,
                      standard_age=standard_LT_age,
                      prior_population=prior_population,
                      nTrials=nTrials,
                      uncertainty_range=FALSE,
                      
                      
                      uncertainty_alpha=uncertainty_alpha,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                      max_age=max_age,
                      filepath=filepath
      )
      return(SEYLL(SLE_object = SLE_object,
                   personunits=personunits,
				   filepath=filepath,
                   Rate=Rate,
                   Beta=Beta,
                   Const=Const,
                   Agewt=Agewt)$YLL)
    })),ncol=ncol(mortality_data))
    YLL_byyear<-data.frame(YLL_byyear)
    names(YLL_byyear)<-paste0('Y',year_range)
    row.names(YLL_byyear)<-c(mortality_age_labels,'Total')
    
    YLL<-list(YLL_byyear=YLL_byyear)
  } else {
    SEYLL_object_by_year<-list()
    for(i in 1:ncol(mortality_data)){
      SLE_object<-SLE(age_labels=mortality_age_labels,
                      population=population_data[,i],
                      death_counts=mortality_data[,i],
                      age_average_at_death=age_average_at_death_final[,i],
                      standard_LE=standard_LE,
                      standard_age=standard_LT_age,
                      prior_population=prior_population,
                      nTrials=nTrials,
                      uncertainty_range=TRUE,
                      uncertainty_alpha=uncertainty_alpha,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                      max_age=max_age,
                      filepath=filepath
      )
      SEYLL_object_by<-SEYLL(SLE_object = SLE_object,
                             personunits=personunits,
							 filepath=filepath,
                             Rate=Rate,
                             Beta=Beta,
                             Const=Const,
                             Agewt=Agewt)
      SEYLL_object_by_year[[i]]<-SEYLL_object_by
    }
    
    names(SEYLL_object_by_year)<-paste0('Y',year_range)
    
    YLL_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(SEYLL_object_by_year),FUN=function(i){
      SEYLL_object_by_year[[i]]$YLL
    })),ncol=ncol(mortality_data)))
    
    YLL_lr_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(SEYLL_object_by_year),FUN=function(i){
      SEYLL_object_by_year[[i]]$YLL_lr
    })),ncol=ncol(mortality_data)))
    
    YLL_up_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(SEYLL_object_by_year),FUN=function(i){
      SEYLL_object_by_year[[i]]$YLL_up
    })),ncol=ncol(mortality_data)))
    
    
    ### personunits
    YLLper_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(SEYLL_object_by_year),FUN=function(i){
      SEYLL_object_by_year[[i]]$YLL_perUnit
    })),ncol=ncol(mortality_data)))
    
    YLLper_lr_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(SEYLL_object_by_year),FUN=function(i){
      SEYLL_object_by_year[[i]]$YLL_perUnit_lr
    })),ncol=ncol(mortality_data)))
    
    YLLper_up_byyear<-data.frame(matrix(unlist(parallel::mclapply(1:length(SEYLL_object_by_year),FUN=function(i){
      SEYLL_object_by_year[[i]]$YLL_perUnit_up
    })),ncol=ncol(mortality_data)))
    
    
    
    ### defined names and cols
    names(YLL_byyear)<-paste0('Y',year_range)
    row.names(YLL_byyear)<-c(mortality_age_labels,'Total')
    
    names(YLL_lr_byyear)<-paste0('Y',year_range)
    row.names(YLL_lr_byyear)<-c(mortality_age_labels,'Total')
    
    names(YLL_up_byyear)<-paste0('Y',year_range)
    row.names(YLL_up_byyear)<-c(mortality_age_labels,'Total')
    
    names(YLLper_byyear)<-paste0('Y',year_range)
    row.names(YLLper_byyear)<-c(mortality_age_labels,'Total')
    
    names(YLLper_lr_byyear)<-paste0('Y',year_range)
    row.names(YLLper_lr_byyear)<-c(mortality_age_labels,'Total')
    
    names(YLLper_up_byyear)<-paste0('Y',year_range)
    row.names(YLLper_up_byyear)<-c(mortality_age_labels,'Total')
    
	print('========================================') 
	message(paste0('YLD rates were calculated in the unit of ',personunits,'.'))
	print('========================================') 
	
    YLL<-list(YLL_byyear=YLL_byyear,
              YLL_lr_byyear=YLL_lr_byyear,
              YLL_up_byyear=YLL_up_byyear,
              YLLper_byyear=YLLper_byyear,
              YLLper_lr_byyear=YLLper_lr_byyear,
              YLLper_up_byyear=YLLper_up_byyear,
              SEYLL_object=SEYLL_object_by_year)
    names(YLL)[4:6]<-paste0(rep('YLLper',times=3),#personunits,
							c('_byyear','_lr_byyear','_up_byyear'))
  }
  
  
  class(YLL) <- "YLL_object"
  attr(YLL,"Call") <- sys.call()
  return(YLL)
  
  }