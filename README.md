# LBOD
A burden of disease packages, developed specifically for localizing BOD study

token: 1a04ae09ee6c23af26162f643094392697b684b6

useage below:
devtools::install_github('al00014/LBOD',auth_token='1a04ae09ee6c23af26162f643094392697b684b6')

## bug fix to the working example
get_burden(disease='lung',
                    year_range=2004:2011,
                    gender='female',
                    mortality_data=mortality_cases_bycause_female,
                    incident_data=incident_cases_bycause_female,
				           input_list=TRUE,
                    population_std=population_std,
                    age_average_at_death=age_average_at_death,
                    input_age_at_onset=age_at_onset,
                    standard_LE=WHO_LE_table[,2],
                    standard_LT_age=WHO_LE_table[,1],  
                    input_duration=c(0.63310,input_duration[,1]),
                    input_Duration_interval=rbind(data.frame(Duration_lr=0.59965,Duration_up=2.43565),input_duration[,2:3]),
                    input_DisabilityWeight=input_DisabilityWeight,
                    input_DisabilityWeight_interval=input_DisabilityWeight_interval)
