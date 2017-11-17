#setwd("./BOD")

#' get_burden function
#'
#' This function estimates YLL, incident-YLD and DALY and their rates
#' based on input of mortality and incident cases.
#' Rates are calculated based on input of argument population_std. 
#' YLL also requires input of age_average_at_death and standard lifetable.
#' YLD also requires input of input_age_at_onset, input_duration, input_Duration_interval,
#' input_DisabilityWeight and input_DisabilityWeight_interval.
#' Current default lifetable is Coale and Demeny's West level 26.
#' Mortality and incident cases should be in tabular form, row by age and column by year.
#' The default age group should be 19 groups, starting from 0~1, 1~4, and then by 5 years up to 85+. 
#' 
#' @param disease a character indicator of which disease to extract, use when input_list=TRUE, the disease names can be found in description.
#' @param gender a character indicator for gender.
#' @param year_range a vector showing the range of years.
#' @param incident_data a list of diseases or a dataframe of disease. Each disease is represented in tabular data with row by age and column by year. It requires input of incidence cases.
#' @param mortality_data a list of diseases or a dataframe of disease. Each disease is represented in tabular data with row by age and column by year. It requires input of mortality cases.
#' @param population_std a dataframe of population, with row by age and column by year. 
#' @param input_list indicates whether mortality_data and incident_data is a list by itself, the example dataset is a list sorted by cancer names, so the (default)input_list=TRUE. 
#' @param ...
#' @return A list of tabular YLL, YLD, DALY and their rates.
#' @keywords incident_data, mortality_data
#' @export
#' @examples
#' get_burden(disease='lung',
#'                     year_range=2004:2011,
#'                     gender='female',
#'                     mortality_data=mortality_cases_bycause_female,
#'                     incident_data=incident_cases_bycause_female,
#'					           input_list=TRUE,
#'                     population_std=population_std,
#'                     age_average_at_death=age_average_at_death,
#'                     input_age_at_onset=age_at_onset,
#'                     standard_LE=WHO_LE_table[,2],
#'                     standard_LT_age=WHO_LE_table[,1],  
#'                     input_duration=input_duration[,1],
#'                     input_Duration_interval=input_duration[,2:3],
#'                     input_DisabilityWeight=input_DisabilityWeight,
#'                     input_DisabilityWeight_interval=input_DisabilityWeight_interval)
get_burden<-function(disease='lung',
                     gender='male',
                     year_range=2004:2011,
                     mortality_data,#=mortality_cases_bycause_male,
                     incident_data,#=incident_cases_bycause_male,
                     input_list=TRUE,
                     #population=population_GZ2004_2011$male,
                     population_std,#=China_pop_GZ2004_2011$male,
                     age_average_at_death,#=rowMeans(age_at_onset_df),
                     input_age_at_onset,#=rowMeans(age_at_onset_df),
                     standard_LE,#=WHO_LE_table[,2],
                     standard_LT_age,#=WHO_LE_table[,1],  
                     input_duration,#=c(rep(0,13),rep(10,2),rep(5,2),rep(3,2)),
                     input_DisabilityWeight,#=rep(0.54,times=length(age_labels)),
                     input_Duration_interval,#=data.frame(duration-duration*1/4,
                     #           duration+duration*1/4),
                     input_DisabilityWeight_interval,#=DisabilityWeight_interval,
                     age_label=c('0'  ,   '1-4' ,  '5-9'  , 
                                 '10-14', '15-19' ,'20-24' ,
                                 '25-29' ,'30-34', '35-39' ,
                                 '40-44' ,'45-49' ,'50-54' ,
                                 '55-59' ,'60-64', '65-69' ,
                                 '70-74' ,'75-79', '80-84', '85+'  ),
                     uncertainty_alpha=0.05,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                     prior_population=0.001,
                     nTrials=2000,
                     YLD_perUnit=1000,
                     personunits=1000,
                     RateUNIT=1000,
                     Rate=0.03,
                     Beta=0.04,
                     Const=0.1658,
                     Agewt=0,
                     max_age_YLL=100,
                     max_age_YLD=90#,# 'A:/r references/R functions/Disease_Burden_projection_project'
                     #filepath=paste0('./R') 
					 ){
  if(YLD_perUnit!=personunits){
    stop('\n
         To calculate DALY, YLL and YLD should be in the same unit,
         Just make YLD_perUnit==personunits')
  }
  output_rate_unit<-YLD_perUnit<-personunits
  print('---------------------------------------------------')
  message(paste0('Burden calculation for gender: ',gender))
  print('---------------------------------------------------')
  source(paste0('./Mor_to_YLL.R'))
  source(paste0('./inc_to_incYLD.R'))
  # filepath,'/From_incidence_to_YLD',
  if(input_list==TRUE){
    input_mortality<-mortality_data[[grep(disease,names(mortality_data))]]
    input_incident<-incident_data[[grep(disease,names(incident_data))]]
  } else {
    input_mortality<-mortality_data
    input_incident<-incident_data
  }
  
  
  ## get YLL from mortality, 
  ## input population is only used in standardization, which means that it only influences YLL rate
  YLL_std<-Mor_to_YLL(year_range=year_range,
                      mortality_age_labels=age_label,
                      mortality_data=input_mortality,
                      age_average_at_death=age_average_at_death,#[,1]#[,-1]#
                      standard_LE=standard_LE,
                      standard_LT_age=standard_LT_age, #
                      population=population_std,
                      max_age=max_age_YLL,
                      personunits=personunits,
                      prior_population=prior_population,
                      uncertainty_alpha=uncertainty_alpha,
                      nTrials=nTrials,
                      Rate=Rate,
                      Beta=Beta,
                      Const=Const,
                      Agewt=Agewt,
                      Output_YLL_noUI=FALSE,
                      #filepath=filepath,
                      verbose = FALSE)
  
  
  ## get YLD from incidence, 
  ## input population is only used in standardization, which means that it only influences YLD rate
  YLD_std<-inc_to_incYLD(incident_age_labels=age_label,
                         incident_data=input_incident,
                         input_age_at_onset=input_age_at_onset,
                         input_population=population_std,
                         input_duration=input_duration,
                         input_DisabilityWeight=input_DisabilityWeight,
                         input_Duration_interval=input_Duration_interval,
                         input_DisabilityWeight_interval=input_DisabilityWeight_interval,
                         Input_Incident_Rate=FALSE,
                         year_range=year_range,
                         verbose = FALSE,
                       #  filepath=filepath,
                         uncertainty_alpha=uncertainty_alpha,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                         prior_population=prior_population,
                         nTrials=nTrials,
                         YLD_perUnit=YLD_perUnit,
                         RateUNIT=RateUNIT,
                         Rate=Rate,
                         Beta=Beta,
                         Const=Const,
                         Agewt=Agewt,
                         max_age=max_age_YLD,
                         Output_YLD_noUI=FALSE)
  
  
  DALY<-YLL_std$YLL_byyear+YLD_std$YLD_byyear
  DALY_lr<-YLL_std$YLL_lr_byyear+YLD_std$YLD_lr_byyear
  DALY_up<-YLL_std$YLL_up_byyear+YLD_std$YLD_up_byyear
  
  DALY_rate<-YLL_std$YLLper_byyear+YLD_std$YLDper_byyear
  DALY_rate_lr<-YLL_std$YLLper_lr_byyear+YLD_std$YLDper_lr_byyear
  DALY_rate_up<-YLL_std$YLLper_up_byyear+YLD_std$YLDper_up_byyear
  
  burden_object<-list(gender=gender,
                      unit_for_standardization=output_rate_unit,
                      YLL=list(not_standardized=list(point_est=YLL_std$YLL_byyear,
                                                     lr=YLL_std$YLL_lr_byyear,
                                                     up=YLL_std$YLL_up_byyear),
                               standardized=list(point_est=YLL_std$YLLper_byyear,
                                                 lr=YLL_std$YLLper_lr_byyear,
                                                 up=YLL_std$YLLper_up_byyear)),
                      YLD=list(not_standardized=list(point_est=YLD_std$YLD_byyear,
                                                     lr=YLD_std$YLD_lr_byyear,
                                                     up=YLD_std$YLD_up_byyear),
                               standardized=list(point_est=YLD_std$YLDper_byyear,
                                                 lr=YLD_std$YLDper_lr_byyear,
                                                 up=YLD_std$YLDper_up_byyear)),
                      DALY=list(not_standardized=list(point_est=DALY,
                                                      lr=DALY_lr,
                                                      up=DALY_up),
                                standardized=list(point_est=DALY_rate,
                                                  lr=DALY_rate_lr,
                                                  up=DALY_rate_up)),
                      YLL_object=YLL_std,
                      YLD_object=YLD_std,
                      population_for_std=population_std)
  
  class(burden_object) <- "burden_object"
  attr(burden_object,"Call") <- sys.call()
  return(burden_object)
  
  }

#' get_burden_prev() function
#'
#' This function estimates YLL, prev-YLD and DALY and their rates
#' based on input of mortality and prevalent cases.
#' Rates are calculated based on input of argument population_std. 
#' YLL also requires input of age_average_at_death and standard lifetable.
#' YLD also requires input of input_DisabilityWeight and input_DisabilityWeight_interval.
#' Current default lifetable is Coale and Demeny's West level 26.
#' Mortality and prevalent cases should be in tabular form, row by age and column by year.
#' The default age group should be 19 groups, starting from 0~1, 1~4, and then by 5 years up to 85+. 
#' 
#' @param disease a character indicator of which disease to extract, use when input_list=TRUE, the disease names can be found in description.
#' @param gender a character indicator for gender.
#' @param year_range a vector showing the range of years.
#' @param prev_data a list of diseases or a dataframe of disease. Each disease is represented in tabular data with row by age and column by year. It requires input of prevalence cases.
#' @param mortality_data a list of diseases or a dataframe of disease. Each disease is represented in tabular data with row by age and column by year. It requires input of mortality cases.
#' @param population_std a dataframe of population, with row by age and column by year. 
#' @param input_list indicates whether mortality_data and incident_data is a list by itself, the example dataset is a list sorted by cancer names, so the (default)input_list=TRUE. 
#' @param ...
#' @return A list of tabular YLL, YLD, DALY and their rates.
#' @keywords prev_data, mortality_data
#' @export
#' @examples
#' get_burden_prev(disease='lung',year_range=2004:2011,
#'                     gender='female',
#'                     mortality_data=mortality_cases_bycause_female,
#'                     prev_data=incident_cases_bycause_female,
#'					   input_list=TRUE,
#'                     population_std=population_std,
#'                     age_average_at_death=age_average_at_death,
#'                     standard_LE=WHO_LE_table[,2],
#'                     standard_LT_age=WHO_LE_table[,1],  
#'                     input_DisabilityWeight=input_DisabilityWeight,
#'                     input_DisabilityWeight_interval=input_DisabilityWeight_interval)
get_burden_prev<-function(disease='lung',
                     gender='male',
                     year_range=2004:2011,
                     mortality_data,#=mortality_cases_bycause_male,
                     prev_data,#=#incident_cases_bycause_male,
					 input_list=TRUE,
                     population_std,#=China_pop_GZ2004_2011$male,
                     age_average_at_death,#=rowMeans(age_at_onset_df),
                     standard_LE,#=WHO_LE_table[,2],
                     standard_LT_age,#=WHO_LE_table[,1],  
                     input_DisabilityWeight,#=rep(0.54,times=length(age_labels)),
                     input_DisabilityWeight_interval,#=DisabilityWeight_interval,
                     age_label=c('0'  ,   '1-4' ,  '5-9'  , 
                                 '10-14', '15-19' ,'20-24' ,
                                 '25-29' ,'30-34', '35-39' ,
                                 '40-44' ,'45-49' ,'50-54' ,
                                 '55-59' ,'60-64', '65-69' ,
                                 '70-74' ,'75-79', '80-84', '85+'  ),
                     uncertainty_alpha=0.05,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                     prior_population=0.001,
                     nTrials=2000,
                     YLD_perUnit=1000,
                     personunits=1000,
                     RateUNIT=1000,
                     Rate=0.03,
                     Beta=0.04,
                     Const=0.1658,
                     Agewt=0,
                     max_age_YLL=100,
                     max_age_YLD=90#,
                     #filepath=paste0('./R') 
){
  if(YLD_perUnit!=personunits){
    stop('\n
         To calculate DALY, YLL and YLD should be in the same unit,
         Just make YLD_perUnit==personunits')
  }
  output_rate_unit<-YLD_perUnit<-personunits
  print('---------------------------------------------------')
  message(paste0('Burden calculation for gender: ',gender))
  print('---------------------------------------------------')
  source(paste0(filepath,'/From_mortality_to_YLL','/Mor_to_YLL.R'))
  #source(paste0(filepath,'/From_incidence_to_YLD','/inc_to_incYLD.R'))
  source(paste0(filepath,'/From_prevalence_to_YLD','/prev_to_prevYLD.R'))
  
  if(input_list==TRUE){
    input_mortality<-mortality_data[[grep(disease,names(mortality_data))]]
    input_incident<-prev_data[[grep(disease,names(prev_data))]]
  } else {
	input_mortality<-mortality_data
	input_incident<-prev_data
  }
  
  
  ## get YLL from mortality, 
  ## input population is only used in standardization, which means that it only influences YLL rate
  YLL_std<-Mor_to_YLL(year_range=year_range,
                      mortality_age_labels=age_label,
                      mortality_data=input_mortality,
                      age_average_at_death=age_average_at_death,#[,1]#[,-1]#
                      standard_LE=standard_LE,
                      standard_LT_age=standard_LT_age, #
                      population=population_std,
                      max_age=max_age_YLL,
                      personunits=personunits,
                      prior_population=prior_population,
                      uncertainty_alpha=uncertainty_alpha,
                      nTrials=nTrials,
                      Rate=Rate,
                      Beta=Beta,
                      Const=Const,
                      Agewt=Agewt,
                      Output_YLL_noUI=FALSE,
                      #filepath=filepath ,
                      verbose = FALSE)
  
  
  ## get YLD from incidence, 
  ## input population is only used in standardization, which means that it only influences YLD rate
  YLD_std<-prev_to_prevYLD(prevalent_age_labels=age_label,
                        prevalent_data=input_incident,
                        input_population=population_std, # population should be annual, and should have the same annual data as incident data
                        input_DisabilityWeight=input_DisabilityWeight,
                        input_DisabilityWeight_interval=input_DisabilityWeight_interval,
                        year_range=year_range,
                        Output_YLD_noUI=FALSE,# if only output YLD, no uncertainty is analyzed, no further detail information will be output,#
                        #  it should be a quicker way to get a teaser of the output YLD
                        RateUNIT=RateUNIT,
                        Input_Prevalent_Rate=FALSE,
                        uncertainty_alpha=uncertainty_alpha,  ## uncertainty range default to 0.95% CI, indicating an alpha of 0.05, it is on both tails
                        prior_population=prior_population,
                        nTrials=nTrials,
                        YLD_perUnit=YLD_perUnit,
                        max_age=max_age_YLD,
                        #filepath=filepath ,
                        verbose=FALSE)
  
  
  
  
  DALY<-YLL_std$YLL_byyear+YLD_std$YLD_byyear
  DALY_lr<-YLL_std$YLL_lr_byyear+YLD_std$YLD_lr_byyear
  DALY_up<-YLL_std$YLL_up_byyear+YLD_std$YLD_up_byyear
  
  DALY_rate<-YLL_std$YLLper_byyear+YLD_std$YLDper_byyear
  DALY_rate_lr<-YLL_std$YLLper_lr_byyear+YLD_std$YLDper_lr_byyear
  DALY_rate_up<-YLL_std$YLLper_up_byyear+YLD_std$YLDper_up_byyear
  
  burden_object<-list(gender=gender,
                      unit_for_standardization=output_rate_unit,
                      YLL=list(not_standardized=list(point_est=YLL_std$YLL_byyear,
                                                     lr=YLL_std$YLL_lr_byyear,
                                                     up=YLL_std$YLL_up_byyear),
                               standardized=list(point_est=YLL_std$YLLper_byyear,
                                                 lr=YLL_std$YLLper_lr_byyear,
                                                 up=YLL_std$YLLper_up_byyear)),
                      YLD=list(not_standardized=list(point_est=YLD_std$YLD_byyear,
                                                     lr=YLD_std$YLD_lr_byyear,
                                                     up=YLD_std$YLD_up_byyear),
                               standardized=list(point_est=YLD_std$YLDper_byyear,
                                                 lr=YLD_std$YLDper_lr_byyear,
                                                 up=YLD_std$YLDper_up_byyear)),
                      DALY=list(not_standardized=list(point_est=DALY,
                                                      lr=DALY_lr,
                                                      up=DALY_up),
                                standardized=list(point_est=DALY_rate,
                                                  lr=DALY_rate_lr,
                                                  up=DALY_rate_up)),
                      YLL_object=YLL_std,
                      YLD_object=YLD_std,
                      population_for_std=population_std)
  
  class(burden_object) <- "burden_object"
  attr(burden_object,"Call") <- sys.call()
  return(burden_object)
  
  }
  
#### function for recalculating the burden rate

#' recal_rate_forBothGen() function
#'
#' This function works with objects produced by get_burden() or get_burden_prev().
#' It recalculates burden rate for both gender, based on input of argument population.
#' 
#' @param burden_male a burden_object produced by get_burden() or get_burden_prev(). Only take in the male results.
#' @param burden_female a burden_object produced by get_burden() or get_burden_prev(). Only take in the female results.
#' @param output_rate_unit a value indicating the unit of rate.
#' @param population a list of population. Only two elements should be in the list, each represents a specific gender. Accepted values can be 'male' or 'female'.
#' @param ...
#' @return A list of tabular YLL, YLD, DALY and their rates, rates standardized by input of population.
#' @export
#' @examples
#' LC_burden_male<-get_burden(disease='lung',year_range=2004:2011,
#'                     gender='male',
#'                     mortality_data=mortality_cases_bycause_male,
#'                     incident_data=incident_cases_bycause_male,
#'					   input_list=TRUE,
#'                     population_std=population_std,
#'                     age_average_at_death=age_average_at_death,
#'                     input_age_at_onset=age_at_onset,
#'                     standard_LE=WHO_LE_table[,2],
#'                     standard_LT_age=WHO_LE_table[,1],  
#'                     input_duration=input_duration[,1],
#'                     input_Duration_interval=input_duration[,2:3],
#'                     input_DisabilityWeight=input_DisabilityWeight,
#'                     input_DisabilityWeight_interval=input_DisabilityWeight_interval)
#' LC_burden_female<-get_burden(disease='lung',year_range=2004:2011,
#'                     gender='female',
#'                     mortality_data=mortality_cases_bycause_female,
#'                     incident_data=incident_cases_bycause_female,
#'					   input_list=TRUE,
#'                     population_std=population_std,
#'                     age_average_at_death=age_average_at_death,
#'                     input_age_at_onset=age_at_onset,
#'                     standard_LE=WHO_LE_table[,2],
#'                     standard_LT_age=WHO_LE_table[,1],  
#'                     input_duration=input_duration[,1],
#'                     input_Duration_interval=input_duration[,2:3],
#'                     input_DisabilityWeight=input_DisabilityWeight,
#'                     input_DisabilityWeight_interval=input_DisabilityWeight_interval)
#'recal_rate_forBothGen(burden_male=LC_burden_male,# takes in a burden_object
#'                                burden_female=LC_burden_female,# takes in a burden_object
#'                                output_rate_unit=100000,
#'                                population=population_ex)
recal_rate_forBothGen<-function(burden_male,#=LC_burden,# takes in a burden_object
                                burden_female,#=LC_burden_fe,# takes in a burden_object
                                output_rate_unit=1000,
                                population#=China_pop_GZ2004_2011
								){
  #if(sum(!c('parallel') %in% rownames(installed.packages()) )!= 0) {
  #  warning('\n
  #       Required packages not found, installing now.
  #		 the dependent package(s) is(are):
  #       parallel')
  #	installed.packages('parallel')	 
  #}
  if(class(burden_male)!='burden_object'){
    
    stop('### must take-in a burden_object object produced by get_burden() function')
  }
  if(class(burden_female)!='burden_object'){
    
    stop('### must take-in a burden_object object produced by get_burden() function')
  }
  if(sum(names(population)!=c('male','female'))!=0){
    stop('This function requires input population from both genders.
         Please check the population argument')
  }
  rowNo<-length(dimnames(burden_male$DALY$not_standardized$point_est)[[1]])
  colNo<-length(dimnames(burden_male$DALY$not_standardized$point_est)[[2]])
  
  #require(parallel)
  burden_BothGen_final<-list()
  looper<-1
  for(i in 3:5){
    burden_BothGen_lst<-parallel::mclapply(1:3,FUN = function(j){
      burden_bothGen<-burden_male[[i]]$not_standardized[[j]]+burden_female[[i]]$not_standardized[[j]]
    })
    names(burden_BothGen_lst)<-paste0('not_standardized_',c('point_est','lr','up'))
    burden_BothGen_final[[looper]]<-burden_BothGen_lst
    looper<-looper+1
  }
  names(burden_BothGen_final)<-c('YLL','YLD','DALY')
  
  population_bothGen<-population$male+population$female
  burden_rate_BothGen_final<-list()
  for(i in 1:3){
    burden_BothGen_lst<-parallel::mclapply(1:3,FUN = function(j){
      burden_bothGen<-burden_BothGen_final[[i]][[j]][-rowNo,]/population_bothGen*output_rate_unit
      
      Total_rate<-colSums(burden_BothGen_final[[i]][[j]][-rowNo,])/colSums(population_bothGen)*output_rate_unit
      
      burden_bothGen_withTotal<-rbind(burden_bothGen,Total_rate)
      row.names(burden_bothGen_withTotal)[rowNo]<-'Total'
      return(burden_bothGen_withTotal)
    })
    names(burden_BothGen_lst)<-paste0('standardized_',c('point_est','lr','up'))
    burden_rate_BothGen_final[[i]]<-burden_BothGen_lst
    
  }
  names(burden_rate_BothGen_final)<-c('YLL','YLD','DALY')
  
  burden_object_BothGen<-list(YLL=list(not_standardized=list(point_est=burden_BothGen_final$YLL$not_standardized_point_est,
                                                             lr=burden_BothGen_final$YLL$not_standardized_lr,
                                                             up=burden_BothGen_final$YLL$not_standardized_up),
                                       standardized=list(point_est=burden_rate_BothGen_final$YLL$standardized_point_est,
                                                         lr=burden_rate_BothGen_final$YLL$standardized_lr,
                                                         up=burden_rate_BothGen_final$YLL$standardized_up)),
                              YLD=list(not_standardized=list(point_est=burden_BothGen_final$YLD$not_standardized_point_est,
                                                             lr=burden_BothGen_final$YLD$not_standardized_lr,
                                                             up=burden_BothGen_final$YLD$not_standardized_up),
                                       standardized=list(point_est=burden_rate_BothGen_final$YLD$standardized_point_est,
                                                         lr=burden_rate_BothGen_final$YLD$standardized_lr,
                                                         up=burden_rate_BothGen_final$YLD$standardized_up)),
                              DALY=list(not_standardized=list(point_est=burden_BothGen_final$DALY$not_standardized_point_est,
                                                              lr=burden_BothGen_final$DALY$not_standardized_lr,
                                                              up=burden_BothGen_final$DALY$not_standardized_up),
                                        standardized=list(point_est=burden_rate_BothGen_final$DALY$standardized_point_est,
                                                          lr=burden_rate_BothGen_final$DALY$standardized_lr,
                                                          up=burden_rate_BothGen_final$DALY$standardized_up)),
                              unit_for_standardization=output_rate_unit,
                              population_for_std=population
  )
  
  class(burden_object_BothGen) <- "burden_object_forBothGen"
  attr(burden_object_BothGen,"Call") <- sys.call()
  return(burden_object_BothGen)
  }