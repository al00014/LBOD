#' get_specific_agegroup_burden() function
#'
#' A helper function to sort the data into desired number of age groups.
#' The function can only take in burden_object produced by get_burden() function.
#' 
#' @param data input a burden_object object produced by get_burden() function
#' @param target_age_group_row it indicates which rows to be combined
#' @param target_age_group_lab it indentifies the combined rows with human-readable labels
#' @param ...
#' @return A list of tabular YLL, YLD, DALY and their rates, with re-grouped age.
#' @keywords data,target_age_group_row,target_age_group_lab
#' @export

get_specific_agegroup_burden<-function(data,#=LC_burden_std_GZ,
                                       target_age_group_row=list(c(6:10),
                                                                 c(11:14),
                                                                 c(15:19)),
                                       target_age_group_lab=c('20-44','45-64','65-89'),
                                       perunit=1000,
                                       digits=2){
  if(class(data)!='burden_object'){
    
    stop('### must take-in a burden_object object produced by get_burden() function')
  } else {
    YLL<-data$YLL$not_standardized
    YLD<-data$YLD$not_standardized
    DALY<-data$DALY$not_standardized
    population<-data$population_for_std

    ############################# get YLL and YLD in new age group    
    YLL_reAge_ls<-list()
    YLD_reAge_ls<-list()
    
    for(i in 1:length(YLL)){
      YLL_reAge<-data.frame()
      YLD_reAge<-data.frame()
      
      for(j in 1:length(target_age_group_row)){
        temp1<-colSums(YLL[[i]][target_age_group_row[[j]],])
        temp2<-colSums(YLD[[i]][target_age_group_row[[j]],])
        
        YLL_reAge<-rbind(YLL_reAge,temp1)
        YLD_reAge<-rbind(YLD_reAge,temp2)
        
      }
      YLL_reAge<-round(YLL_reAge,digits = digits)
      YLD_reAge<-round(YLD_reAge,digits = digits)
      
      row.names(YLL_reAge)<-row.names(YLD_reAge)<-target_age_group_lab
      colnames(YLL_reAge)<-colnames(YLD_reAge)<-colnames(data$YLL$not_standardized$point_est)
      YLL_reAge_ls[[i]]<-YLL_reAge
      
      YLD_reAge_ls[[i]]<-YLD_reAge
    }
    names(YLL_reAge_ls)<-c('point_est','lr','up')
    names(YLD_reAge_ls)<-c('point_est','lr','up')
    
    ############################ resort age group in population
    pop_reAge<-data.frame()
    for(i in 1:length(target_age_group_row)){
      pop_temp<-colSums(population[target_age_group_row[[i]],])
      pop_reAge<-rbind(pop_reAge,pop_temp)
    }

    row.names(pop_reAge)<-target_age_group_lab
    colnames(pop_reAge)<-colnames(data$YLL$not_standardized$point_est)

    ############################# get DALY in new age group    
    DALY_reAge_ls<-list()
    for(i in 1:length(YLL)){
      DALY_reAge_ls[[i]]<-YLL_reAge_ls[[i]]+YLD_reAge_ls[[i]]
      
    }
    names(DALY_reAge_ls)<-c('point_est','lr','up')
    
    ############################ get total
    for(i in 1:length(YLL)){
      ### YLL
      YLL_reAge_ls[[i]]<-rbind(YLL_reAge_ls[[i]],colSums(YLL_reAge_ls[[i]])) 
      row.names(YLL_reAge_ls[[i]])<-c(target_age_group_lab,'Total')
      
      ### YLD
      YLD_reAge_ls[[i]]<-rbind(YLD_reAge_ls[[i]],colSums(YLD_reAge_ls[[i]])) 
      row.names(YLD_reAge_ls[[i]])<-c(target_age_group_lab,'Total')
      
      ### DALY
      DALY_reAge_ls[[i]]<-rbind(DALY_reAge_ls[[i]],colSums(DALY_reAge_ls[[i]])) 
      row.names(DALY_reAge_ls[[i]])<-c(target_age_group_lab,'Total')
    }
    pop_reAge<-rbind(pop_reAge,colSums(pop_reAge))
    row.names(pop_reAge)<-c(target_age_group_lab,'Total')
    
    
    ############################ get rates for YLL, YLD and DALY
    YLLrate_reAge_ls<-list()
    YLDrate_reAge_ls<-list()
    DALYrate_reAge_ls<-list()
    for(i in 1:length(YLL)){
      YLLrate_reAge_ls[[i]]<-YLL_reAge_ls[[i]]/pop_reAge*perunit
      YLDrate_reAge_ls[[i]]<-YLD_reAge_ls[[i]]/pop_reAge*perunit
      DALYrate_reAge_ls[[i]]<-DALY_reAge_ls[[i]]/pop_reAge*perunit
    }
    
    names(YLLrate_reAge_ls)<-names(YLDrate_reAge_ls)<-names(DALYrate_reAge_ls)<-c('point_est','lr','up')
    
  }
  print('========================================')
  message(paste0('The dataset were regrouped by the following age groups: (years)'))
  for(i in 1:length(target_age_group_lab)){
    message(target_age_group_lab[i])
  }
  message(paste0('Total were also calculated for the re-age dataset.'))
 # message(paste(target_age_group_lab))
  print('========================================') 
  message(paste0('Rate was calculated in the unit of ',perunit))
  
  
  return_obj<-list(not_rate=list(YLL=YLL_reAge_ls,
                                 YLD=YLD_reAge_ls,
                                 DALY=DALY_reAge_ls),
                   rate=list(YLLrate=YLLrate_reAge_ls,
                             YLDrate=YLDrate_reAge_ls,
                             DALYrate=DALYrate_reAge_ls),
                   perunit=perunit,
                   population=pop_reAge)
  #LC_burden_std_GZ$YLL$standardized
  #LC_burden_std_China$YLL$standardized
  class(return_obj) <- "reAge_datatable_object"
  attr(return_obj,"Call") <- sys.call()
  return(return_obj)
}
#' get_specific_year_burden() function
#'
#' A helper function to sort the data into desired number of periods of year. 
#' The function can only take in reAge_datatable_object produced by get_specific_agegroup_burden() function.
#' 
#' @param data input a reAge_datatable_object object produced by get_specific_agegroup_burden() function
#' @param target_year_col it indicates which columns to be combined
#' @param target_year_col_lab it indentifies the combined columns with human-readable labels
#' @param ...
#' @return A list of tabular YLL, YLD, DALY and their rates, with re-grouped age and year.
#' @keywords data,target_year_col,target_year_col_lab
#' @export
get_specific_year_burden<-function(data,#data=BC_reAge_byCH,
                                   target_year_col=list(c(1:3),
                                                             c(4:6),
                                                             c(7:8),
                                                        c(9:11),
                                                        c(12:14),
                                                        c(15:17),
                                                        c(18:20),
                                                        c(21:23),
                                                        c(24:27)),
                                   target_year_col_lab=c('2004_2006','2007_2009','2010_2011',
                                                     '2012_2014','2015_2017','2018_2020',
                                                     '2021_2023','2024_2026','2027_2030'),
                                   perunit=1000,
                                   digits=2){
  if(class(data)!='reAge_datatable_object'){
    
    stop('### must take-in a reAge_datatable_object object produced by get_specific_agegroup_burden() function')
  } else {
  output_rowname<-row.names(data$not_rate$YLL$point_est)
  output_colnames<-paste0('Y',target_year_col_lab)
  YLL<-data$not_rate$YLL
  YLD<-data$not_rate$YLD
  population<-data$population
  
  #################### resort year group in YLL and YLD
  YLL_reYear_ls<-list()
  YLD_reYear_ls<-list()
  for(i in 1:length(YLL)){
    YLL_reYear<-data.frame()
    YLD_reYear<-data.frame()
    
    for(j in 1:length(target_year_col)){
      tran_YLL<-t(YLL[[i]])
      temp1<-colSums(tran_YLL[target_year_col[[j]],])
      YLL_reYear<-rbind(YLL_reYear,temp1)
      
      tran_YLD<-t(YLD[[i]])
      temp2<-colSums(tran_YLD[target_year_col[[j]],])
      YLD_reYear<-rbind(YLD_reYear,temp2)
    }
    YLL_reYear<-t(YLL_reYear)
    YLD_reYear<-t(YLD_reYear)
    
    row.names(YLL_reYear)<-output_rowname
    colnames(YLL_reYear)<-output_colnames
    
    row.names(YLD_reYear)<-output_rowname
    colnames(YLD_reYear)<-output_colnames
    
    YLL_reYear_ls[[i]]<-YLL_reYear
    YLD_reYear_ls[[i]]<-YLD_reYear
    }
  
  names(YLL_reYear_ls)<-names(YLD_reYear_ls)<-c('point_est','lr','up')
  
  ############################ resort age group in population
  pop_reYear<-data.frame()
  for(i in 1:length(target_year_col)){
    population_t<-t(population)
    pop_temp<-colSums(population_t[target_year_col[[i]],])
    pop_reYear<-rbind(pop_reYear,pop_temp)
  }
  pop_reYear<-t(pop_reYear)
  row.names(pop_reYear)<-output_rowname
  colnames(pop_reYear)<-output_colnames
  
  ########################### re-calculate DALY in new year group
  
   
  DALY_reYear_ls<-list()
  for(i in 1:length(YLL)){
    DALY_reYear_ls[[i]]<-YLL_reYear_ls[[i]]+YLD_reYear_ls[[i]]
    
  }
  names(DALY_reYear_ls)<-c('point_est','lr','up')
  
  ########################### get rates for YLL, YLD and DALY, sorted by new year group
  ############################ get rates for YLL, YLD and DALY
  YLLrate_reYear_ls<-list()
  YLDrate_reYear_ls<-list()
  DALYrate_reYear_ls<-list()
  for(i in 1:length(YLL)){
    YLLrate_reYear_ls[[i]]<-round(YLL_reYear_ls[[i]]/pop_reYear*perunit,digits = digits)
    YLDrate_reYear_ls[[i]]<-round(YLD_reYear_ls[[i]]/pop_reYear*perunit,digits = digits)
    DALYrate_reYear_ls[[i]]<-round(DALY_reYear_ls[[i]]/pop_reYear*perunit,digits = digits)
  }
  
  names(YLLrate_reYear_ls)<-names(YLDrate_reYear_ls)<-names(DALYrate_reYear_ls)<-c('point_est','lr','up')
  
  print('========================================')
  message(paste0('The dataset were regrouped by the following year groups:'))
  for(i in 1:length(target_year_col_lab)){
    message(target_year_col_lab[i])
  }
  message(paste0('Total were also calculated for the re-year dataset.'))
  # message(paste(target_age_group_lab))
  print('========================================') 
  message(paste0('Rate was calculated in the unit of ',perunit))
  
  }
  return_obj<-list(not_rate=list(YLL=YLL_reYear_ls,
                                 YLD=YLD_reYear_ls,
                                 DALY=DALY_reYear_ls),
                   rate=list(YLLrate=YLLrate_reYear_ls,
                             YLDrate=YLDrate_reYear_ls,
                             DALYrate=DALYrate_reYear_ls),
                   perunit=perunit,
                   population=pop_reYear)
  
  class(return_obj) <- "reYear_datatable_object"
  attr(return_obj,"Call") <- sys.call()
  return(return_obj)
}



#' get_specific_agegroup_cases() function
#'
#' A helper function to sort the cases into desired number of age groups.
#' The function can only take in dataframe of cases (incidence, prevalence or mortality).
#' 
#' @param cases_data input a dataframe of cases (mortality, incidence or prevalence, all in tabular form).
#' @param population_data input a dataframe of population, this population is the local population, for example, a study conducted in Guangzhou, than the Guangzhou population would be used as the population_data.
#' @param population_weighted_data input a dataframe of standardized population. Input of this population depends on study setting, in China, the national population or the world population is often used as a standard population.
#' @param target_age_group_row it indicates which rows to be combined.
#' @param target_age_group_lab it indentifies the combined rows with human-readable labels.
#' @param desired_year indicates which periods of years to extract from the cases_data
#' @param ...
#' @return A list of tabular cases and case rates, with re-grouped age.
#' @keywords cases_data,population_data,population_weighted_data
#' @export

get_specific_agegroup_cases<-function(cases_data, # cases_data=ALLcause_mor_female_projection$datatable[,1:length(2004:2011)]
                                       population_data,#population_data=GZ_population_GAMpred_to2030_correct$female[,1:length(2004:2011)],
                                       population_weighted_data,#=China_population_GAMpred_to2030$female,
                                       target_age_group_row=list(c(6:10),
                                                                 c(11:14),
                                                                 c(15:19)),
                                       desired_year=2004:2011,
                                       age_adjusted=TRUE,
                                       target_age_group_lab=c('20-44','45-64','65-89'),
                                       perunit=1000,
                                       digits=2){
  
  cases<-list(cases_data$datatable[,paste0('Y',desired_year)],
              cases_data$datatable_lr[,paste0('Y',desired_year)],
              cases_data$datatable_up[,paste0('Y',desired_year)])
  population<-population_data[,paste0('Y',desired_year)]
  population_std<-population_weighted_data[,paste0('Y',desired_year)]
  output_colnames<-colnames(population_std)
  cases_reAge_ls<-list()
  for(i in 1:length(cases)){
    YLL_reYear<-data.frame()
     
    
    for(j in 1:length(target_age_group_row)){
      #tran_YLL<-t(YLL[[i]])
      temp1<-colSums(cases[[i]][target_age_group_row[[j]],])
      YLL_reYear<-rbind(YLL_reYear,temp1)
       
    }
    #YLL_reYear<-t(YLL_reYear)
     
    
    row.names(YLL_reYear)<-target_age_group_lab
    colnames(YLL_reYear)<-output_colnames
     
    YLL_reYear<-rbind(YLL_reYear,colSums(YLL_reYear))
    row.names(YLL_reYear)<-c(target_age_group_lab,'Total')
    
    cases_reAge_ls[[i]]<-YLL_reYear
    
  }
  
  names(cases_reAge_ls)<-c('point_est','lr','up')
  
  ############################ resort age group in population
  
  ### here is the population for denominator
  pop_reAge<-data.frame()
  for(i in 1:length(target_age_group_row)){
    pop_temp<-colSums(population[target_age_group_row[[i]],])
    pop_reAge<-rbind(pop_reAge,pop_temp)
  }
  
  row.names(pop_reAge)<-target_age_group_lab
  colnames(pop_reAge)<-output_colnames
  pop_reAge<-rbind(pop_reAge,colSums(pop_reAge))
  row.names(pop_reAge)<-c(target_age_group_lab,'Total')
  
  ############################# below is just the standard population
  pop_reAge_std<-data.frame()
  for(i in 1:length(target_age_group_row)){
    pop_temp<-colSums(population_std[target_age_group_row[[i]],])
    pop_reAge_std<-rbind(pop_reAge_std,pop_temp)
  }
  
  row.names(pop_reAge_std)<-target_age_group_lab
  colnames(pop_reAge_std)<-output_colnames
  pop_reAge_std<-rbind(pop_reAge_std,colSums(pop_reAge_std))
  row.names(pop_reAge_std)<-c(target_age_group_lab,'Total')
  
   
  
  Age_weights<-matrix(unlist(lapply(1:(nrow(pop_reAge_std)),FUN = function(i){
    c(t(pop_reAge_std[i,]))/c(t(colSums(population_std)))
  })),ncol=length(desired_year),byrow=TRUE)
  row.names(Age_weights)<-c(target_age_group_lab,'Total')
  colnames(Age_weights)<-output_colnames
  
  if(age_adjusted==TRUE){
    Casesrate_reAge_ls<-list()
    for(i in 1:length(cases)){
      
      Casesrate_reAge_ls[[i]]<-round(cases_reAge_ls[[i]]/pop_reAge*Age_weights*perunit,digits =digits )
    }
    
    names(Casesrate_reAge_ls)<-c('point_est','lr','up')
    
  } else{
    Casesrate_reAge_ls<-list()
    for(i in 1:length(cases)){
      
      Casesrate_reAge_ls[[i]]<-round(cases_reAge_ls[[i]]/pop_reAge*perunit,digits =digits )
    }
    
    names(Casesrate_reAge_ls)<-c('point_est','lr','up')
  }
  
  
  
  
  print('========================================')
  message(paste0('The dataset were regrouped by the following age groups: (years)'))
  for(i in 1:length(target_age_group_lab)){
    message(target_age_group_lab[i])
  }
  message(paste0('Total were also calculated for the re-age dataset.'))
  # message(paste(target_age_group_lab))
  print('========================================') 
  message(paste0('Rate was calculated in the unit of ',perunit))
  
  
  return_obj<-list(not_rate=cases_reAge_ls,
                   rate=Casesrate_reAge_ls,
                   perunit=perunit,
                   population=pop_reAge,
                   population_std=pop_reAge_std,
                   population_std_TotalbyYear=c(t(colSums(population_std))),
                   Age_weights=Age_weights)
  #LC_burden_std_GZ$YLL$standardized
  #LC_burden_std_China$YLL$standardized
  class(return_obj) <- "reAge_datatable_object"
  attr(return_obj,"Call") <- sys.call()
  return(return_obj)
}

#' get_specific_agegroup_cases() function
#'
#' A helper function to sort the cases into desired number of periods of year.
#' The function can only take in reAge_datatable_object produced by get_specific_agegroup_cases() function.
#' 
#' @param data input a reAge_datatable_object object produced by get_specific_agegroup_cases() function
#' @param target_year_col it indicates which columns to be combined
#' @param target_year_col_lab it indentifies the combined columns with human-readable labels
#' @param ...
#' @return A list of tabular cases and case rates, with re-grouped age and year.
#' @keywords data,target_year_col,target_year_col_lab
#' @export
get_specific_Yeargroup_cases<-function(data,#=AllCause_mort_rate, 
                                      target_year_col=list(c(1:3),
                                                           c(4:6),
                                                           c(7:8)),
                                      target_year_col_lab=c('2004_2006','2007_2009','2010_2011'),
                                      age_adjusted=TRUE,
                                      perunit=1000,
                                      digits=2){
  if(class(data)!='reAge_datatable_object'){
    
    stop('### must take-in a burden_object object produced by get_specific_agegroup_cases() function')
  } else {
    cases<-data$not_rate
    population<-data$population
    population_std<-data$population_std
    populatioN_total<-data$population_std_TotalbyYear
    output_colnames<-colnames(population_std)
    output_rowname<-row.names(population_std)
    
    
    cases_reAge_ls<-list()
    for(i in 1:length(cases)){
      YLL_reYear<-data.frame()
      
      
      for(j in 1:length(target_year_col)){
        tran_YLL<-t(cases[[i]])
        temp1<-colSums(tran_YLL[target_year_col[[j]],])
        YLL_reYear<-rbind(YLL_reYear,temp1)
        
      }
      YLL_reYear<-t(YLL_reYear)
      
      
      row.names(YLL_reYear)<-output_rowname
      colnames(YLL_reYear)<-paste0('Y',target_year_col_lab)
      
      #YLL_reYear<-rbind(YLL_reYear,colSums(YLL_reYear))
      #row.names(YLL_reYear)<-output_rowname
      
      cases_reAge_ls[[i]]<-YLL_reYear
      
    }
    
    names(cases_reAge_ls)<-c('point_est','lr','up')
    
    cases_reYear_ls<-cases_reAge_ls
    
    ############################ resort age group in population
    pop_reYear<-data.frame()
    for(i in 1:length(target_year_col)){
      population_t<-t(population)
      pop_temp<-colSums(population_t[target_year_col[[i]],])
      pop_reYear<-rbind(pop_reYear,pop_temp)
    }
    
    populatioN_total_reyear<-unlist(lapply(1:length(target_year_col),FUN = function(i) sum(populatioN_total[target_year_col[[i]]])))
    
    pop_reYear<-t(pop_reYear)
    row.names(pop_reYear)<-output_rowname
    colnames(pop_reYear)<-paste0('Y',target_year_col_lab)
    
    #### for standard population
    pop_reYear_std<-data.frame()
    for(i in 1:length(target_year_col)){
      population_t<-t(population_std)
      pop_temp<-colSums(population_t[target_year_col[[i]],])
      pop_reYear_std<-rbind(pop_reYear_std,pop_temp)
    }
    pop_reYear_std<-t(pop_reYear_std)
    row.names(pop_reYear_std)<-output_rowname
    colnames(pop_reYear_std)<-paste0('Y',target_year_col_lab)
    
    
    Age_weights<-matrix(unlist(lapply(1:(nrow(pop_reYear_std)),FUN = function(i){
      c(t(pop_reYear_std[i,]))/populatioN_total_reyear
    })),ncol=length(populatioN_total_reyear),byrow=TRUE)
    row.names(Age_weights)<-output_rowname
    colnames(Age_weights)<-paste0('Y',target_year_col_lab)
    
    if(age_adjusted==TRUE){
      Casesrate_reYear_ls<-list()
      for(i in 1:length(cases)){
        
        Casesrate_reYear_ls[[i]]<-round(cases_reYear_ls[[i]]/pop_reYear*Age_weights*perunit,digits =digits )
      }
      
      names(Casesrate_reYear_ls)<-c('point_est','lr','up')
      
    } else{
      Casesrate_reYear_ls<-list()
      for(i in 1:length(cases)){
        
        Casesrate_reYear_ls[[i]]<-round(cases_reYear_ls[[i]]/pop_reYear*perunit,digits =digits )
      }
      
      names(Casesrate_reYear_ls)<-c('point_est','lr','up')
    }
    
  }
  print('========================================')
  message(paste0('The dataset were regrouped by the following age groups: (years)'))
  for(i in 1:length(target_year_col_lab)){
    message(target_year_col_lab[i])
  }
  message(paste0('Total were also calculated for the re-age dataset.'))
  # message(paste(target_year_col_lab))
  print('========================================') 
  message(paste0('Rate was calculated in the unit of ',perunit))
  
  
  return_obj<-list(not_rate=cases_reYear_ls,
                   rate=Casesrate_reYear_ls,
                   perunit=perunit,
                   population=pop_reYear,
                   population_std=pop_reYear_std,
                   population_std_TotalbyYear=populatioN_total_reyear,
                   Age_weights=Age_weights)
  #LC_burden_std_GZ$YLL$standardized
  #LC_burden_std_China$YLL$standardized
  class(return_obj) <- "reYear_datatable_object"
  attr(return_obj,"Call") <- sys.call()
  return(return_obj)
}
