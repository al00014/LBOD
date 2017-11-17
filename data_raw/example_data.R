#if(sum(!c('parallel') %in% rownames(installed.packages()) )!= 0) {
#    warning('\n
#         Required packages not found, installing now.
#		 the dependent package(s) is(are):\n
#         parallel')
#	installed.packages('parallel')	 
#}
#  if(sum(!c(#'xlsx',
#			#'ggplot2',
#			'parallel',
#			#'tidyr',
#			'triangle') %in% rownames(installed.packages()) )!= 0) {
#    warning("\n
#         Required packages not found, installing now.
#		 the dependent package(s) is(are):
#         'parallel','triangle'
#	")
#    #installed.packages('xlsx')
#	#installed.packages('ggplot2')
#	installed.packages('parallel')
#	#installed.packages('tidyr')
#	installed.packages('triangle')
#  }
### function for getting burden for specific input of disease
#setwd("./BOD")
setwd("./BOD")

WHO_LE_table<-readRDS(paste0(getwd(),'/inst/extdata',
                             '/Coale_Demeny_West_level26_lifetable.rds'))
mortality_cases_bycause_male<-readRDS(paste0(getwd(),'/inst/extdata',
                                             '/GZ_tumor_2004to2011_mortality_male_byCancerNames.rds'))
incident_cases_bycause_female<-readRDS(paste0(getwd(),'/inst/extdata',
                                              '/GZ_tumor_2004to2011_incidence_female_byCancerNames.rds'))
incident_cases_bycause_male<-readRDS(paste0(getwd(),'/inst/extdata',
                                            '/GZ_tumor_2004to2011_incidence_male_byCancerNames.rds'))
mortality_cases_bycause_female<-readRDS(paste0(getwd(),'/inst/extdata',
                                               '/GZ_tumor_2004to2011_mortality_female_byCancerNames.rds'))		

population_std<-readRDS(paste0(getwd(),'/inst/extdata',
                               '/World_population_project_2004to2030_byGenders.rds'))$female[,paste0('Y',2004:2011)]
population_ex<-list(male=readRDS(paste0(getwd(),'/inst/extdata',
                                        '/World_population_project_2004to2030_byGenders.rds'))$male[,paste0('Y',2004:2011)],
                    female=readRDS(paste0(getwd(),'/inst/extdata',
                                          '/World_population_project_2004to2030_byGenders.rds'))$female[,paste0('Y',2004:2011)])
#age_labels<-c('0'  ,   '1-4' ,  '5-9'  , 
#              '10-14', '15-19' ,'20-24' ,
#              '25-29' ,'30-34', '35-39' ,
#              '40-44' ,'45-49' ,'50-54' ,
#              '55-59' ,'60-64', '65-69' ,
#              '70-74' ,'75-79', '80-84', '85+'  )
age_at_onset<-readRDS(paste0(getwd(),'/inst/extdata',
                             '/LC_ageonset_female.rds'))$LC_age_onset_mean_female
age_average_at_death<-rowMeans(data.frame(c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85),
                                          c(1,4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,84,100)))
input_duration<-readRDS(paste0(getwd(),'/inst/extdata',
                               '/LC_duration_female.rds'))
input_DisabilityWeight<-rep(0.54,times=19)

input_DisabilityWeight_interval<-data.frame(lr=rep(0.377,19),up=rep(0.687,19))		


devtools::use_data(WHO_LE_table)
devtools::use_data(mortality_cases_bycause_male)
devtools::use_data(incident_cases_bycause_female)
devtools::use_data(incident_cases_bycause_male)
devtools::use_data(mortality_cases_bycause_female)
devtools::use_data(population_std)
devtools::use_data(population_ex)
devtools::use_data(age_at_onset)
devtools::use_data(age_average_at_death)
devtools::use_data(input_duration)
devtools::use_data(input_DisabilityWeight)
devtools::use_data(input_DisabilityWeight_interval)