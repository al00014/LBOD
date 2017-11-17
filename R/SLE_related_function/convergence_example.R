function_mainDIR<-"A:/r references/R functions"
sub_dir<-"Disease_Burden_projection_project"


function_path<-file.path(function_mainDIR,sub_dir)




library(xlsx)
library(parallel)
WHO_LE_table<-read.xlsx(paste0(function_path,'/example_data','/Coale_and_Demeny_level_west_25_26.xlsx'),sheetIndex = 1)
head(WHO_LE_table)
WHO_mortality_data<-read.xlsx(paste0(function_path,'/example_data','/mortality.xlsx'),sheetIndex = 1)

source(paste0(function_path,'/SLE_related_function','/SLE.R'))

SLE(age_labels=WHO_mortality_data[,1],
population<-WHO_mortality_data[,2],
death_counts<-WHO_mortality_data[,3],
age_average_at_death<-WHO_mortality_data[,4],
standard_LE<-WHO_LE_table[,2],  
standard_age<-WHO_LE_table[,1],
uncertainty_range = FALSE)

SLE_object<-SLE(age_labels=WHO_mortality_data[,1],
    population<-WHO_mortality_data[,2],
    death_counts<-WHO_mortality_data[,3],
    age_average_at_death<-WHO_mortality_data[,4],
    standard_LE<-WHO_LE_table[,2],  
    standard_age<-WHO_LE_table[,1],
    uncertainty_alpha = 0.1,
    uncertainty_range = TRUE)


SLE_object$SLE
SLE_object$SLE_lr
SLE_object$SLE_ur
SLE_object$alpha
SLE_object$LE_MC_samples
SLE_object$uncertainty_output
SLE_object$uncertainty_results
SLE_object$Population_uncertainty

str(SLE_object)

source(paste0(function_path,'/SLE_related_function','/convergence_test.R'))
nTrials<-2000

### test mean convergence for the first age group for variable population
start_time<-Sys.time() ## start timing for function below
testingPmean <-unlist(lapply(1:nTrials,FUN = function(i){mean(testing_SLE_convergence(iter = i,data=SLE_object)[[1]][,1],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Pmean <-end_time-start_time

### test mean convergence for the first age group for variable death counts
start_time<-Sys.time() ## start timing for function below
testingDmean <-unlist(lapply(1:nTrials,FUN = function(i){mean(testing_SLE_convergence(iter = i,data=SLE_object)[[1]][,2],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Dmean <-end_time-start_time


### test mean convergence for the first age group for variable average age at death
start_time<-Sys.time() ## start timing for function below
testingAmean <-unlist(lapply(1:nTrials,FUN = function(i){mean(testing_SLE_convergence(iter = i,data=SLE_object)[[1]][,3],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Amean <-end_time-start_time

### test median convergence for the first age group for variable population
start_time<-Sys.time() ## start timing for function below
testingPmedian<-unlist(lapply(1:nTrials,FUN = function(i){median(testing_SLE_convergence(iter = i,data=SLE_object)[[1]][,1],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Pmedian<-end_time-start_time

### test median convergence for the first age group for variable death counts
start_time<-Sys.time() ## start timing for function below
testingDmedian<-unlist(lapply(1:nTrials,FUN = function(i){median(testing_SLE_convergence(iter = i,data=SLE_object)[[1]][,2],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Dmedian<-end_time-start_time


### test median convergence for the first age group for variable average age at death
start_time<-Sys.time() ## start timing for function below
testingAmedian<-unlist(lapply(1:nTrials,FUN = function(i){median(testing_SLE_convergence(iter = i,data=SLE_object)[[1]][,3],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Amedian<-end_time-start_time


plot(x=1:nTrials,testingPmean)
plot(x=1:nTrials,testingDmean)
plot(x=1:nTrials,testingAmean)

plot(x=1:nTrials,testingPmedian)
plot(x=1:nTrials,testingDmedian)
plot(x=1:nTrials,testingAmedian)



save.image(paste0(function_path,'/SLE_related_function','/testing_convergence_meanANDmedian.RData'))



########################################################
source(paste0(function_path,'/YLD_related_function','/YLD.R'))

library(ggplot2)
YLD_norange<-YLD_incident(uncertainty_range=FALSE,
                          age_labels=dat_incidence[,1],
                          population=dat_incidence[,2],
                          incidence=dat_incidence[,3],
                          age_at_onset=dat_incidence[,4],
                          duration=dat_incidence[,5],
                          DisabilityWeight=dat_incidence[,6],
                          YLD_perUnit=100000)
YLD_norange
YLD_norange$YLD
YLD_norange$YLD_per
YLD_norange$YLD_DataFrame
YLD_norange$Original_data

YLD_range<-YLD_incident(uncertainty_range=TRUE,
                          age_labels=dat_incidence[,1],
                          population=dat_incidence[,2],
                          incidence=dat_incidence[,3],
                          age_at_onset=dat_incidence[,4],
                          duration=dat_incidence[,5],
                          DisabilityWeight=dat_incidence[,6],
                          Duration_interval=Duration_interval,
                          DisabilityWeight_interval=DisabilityWeight_interval,
                          prior_population = 0.001,
                          YLD_perUnit=100000)
TempData<-data.frame(Age_group=factor(YLD_range$Age_group),
                     YLD=YLD_range$YLD,
                     YLD_lr=YLD_range$YLD_lr,
                     YLD_up=YLD_range$YLD_up)
TempData<-TempData[-nrow(TempData),]
ReshapeYLD<-reshape(TempData,
                    varying = c(2:4),#sep = "",
                    direction = 'long',
                    v.names=c('YLD'),
                    timevar='Value_groups',
                    idvar = 'Age_group')

#TempData$Age_group
#levels(ReshapeYLD$Age_group)
ReshapeYLD$Age_group<-factor(ReshapeYLD$Age_group,levels = c( "0-4",   "5-14" , "15-29" ,"30-44", "45-59", "60-69" ,"70-79" ,"80+" ))
ReshapeYLD$Value_groups<-factor(ReshapeYLD$Value_groups,levels = c(1,2,3),labels = c('Point','Lower','Upper'))
ggplot(data = ReshapeYLD)+geom_line(aes(x=Age_group,
                                        y=YLD,
                                        #fill=time,
                                        color=as.factor(Value_groups),
                                        group=Value_groups))+geom_point(aes(x=Age_group,
                                                                    y=YLD,
                                                                    #fill=time,
                                                                    color=as.factor(Value_groups),
                                                                    group=Value_groups))

##### convergence test for YLD
source(paste0(function_path,'/YLD_related_function','/convergence_test_YLD.R'))
nTrials<-2000
YLD_object<-YLD_range
### test mean convergence for the first age group for variable population
start_time<-Sys.time() ## start timing for function below
testingPmean <-unlist(lapply(1:nTrials,FUN = function(i){mean(testing_YLD_convergence(iter = i,data=YLD_object)[[8]][,1],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Pmean <-end_time-start_time
plot(x=1:nTrials,testingPmean)

### test mean convergence for the first age group for variable death counts
start_time<-Sys.time() ## start timing for function below
testingImean <-unlist(lapply(1:nTrials,FUN = function(i){mean(testing_YLD_convergence(iter = i,data=YLD_object)[[8]][,2],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Imean <-end_time-start_time


### test mean convergence for the first age group for variable average age at death
start_time<-Sys.time() ## start timing for function below
testingDmean <-unlist(lapply(1:nTrials,FUN = function(i){mean(testing_YLD_convergence(iter = i,data=YLD_object)[[8]][,3],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Dmean <-end_time-start_time

### test mean convergence for the first age group for variable average age at death
start_time<-Sys.time() ## start timing for function below
testingDWmean <-unlist(lapply(1:nTrials,FUN = function(i){mean(testing_YLD_convergence(iter = i,data=YLD_object)[[8]][,4],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_DWmean <-end_time-start_time

### test median convergence for the first age group for variable population
start_time<-Sys.time() ## start timing for function below
testingPmedian<-unlist(lapply(1:nTrials,FUN = function(i){median(testing_YLD_convergence(iter = i,data=YLD_object)[[8]][,1],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Pmedian<-end_time-start_time

### test median convergence for the first age group for variable average age at death
start_time<-Sys.time() ## start timing for function below
testingImedian<-unlist(lapply(1:nTrials,FUN = function(i){median(testing_YLD_convergence(iter = i,data=YLD_object)[[8]][,2],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Imedian<-end_time-start_time

### test median convergence for the first age group for variable death counts
start_time<-Sys.time() ## start timing for function below
testingDmedian<-unlist(lapply(1:nTrials,FUN = function(i){median(testing_YLD_convergence(iter = i,data=YLD_object)[[8]][,3],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_Dmedian<-end_time-start_time


### test median convergence for the first age group for variable death counts
start_time<-Sys.time() ## start timing for function below
testingDWmedian<-unlist(lapply(1:nTrials,FUN = function(i){median(testing_YLD_convergence(iter = i,data=YLD_object)[[8]][,4],na.rm=TRUE)}))
end_time<-Sys.time() # stop timing

time_operate_DWmedian<-end_time-start_time

plot(x=1:nTrials,testingPmean)
plot(x=1:nTrials,testingImean)
plot(x=1:nTrials,testingDmean)
plot(x=1:nTrials,testingDWmean)

plot(x=1:nTrials,testingPmedian)
plot(x=1:nTrials,testingImedian)
plot(x=1:nTrials,testingDmedian)
plot(x=1:nTrials,testingDWmedian)


save.image(paste0(function_path,'/SLE_related_function','/testing_convergence_meanANDmedian2.RData'))







source(paste0(function_path,'/From_incidence_to_YLD','/inc_to_incYLD.R'))
IntoYLD_noUI<-inc_to_incYLD(incident_age_labels<-row.names(incidence_1958to1997),
                       incident_data=incidence_1958to1997,#[,1]
                       input_age_at_onset=matrix(rep(rowMeans(age_ForTurningonset),times=ncol(incidence_1958to1997)),ncol=ncol(incidence_1958to1997))[,1],#[,1],#
                       input_population=population_data_1958to2022,
                       input_duration<-c(rep(0,12),10,10,5,5,3,3),
                       input_DisabilityWeight<-rep(0.5,times=nrow(incidence_1958to1997)),
                       input_Duration_interval<-data.frame(c(rep(0,12),5,5,2.5,2.5,1.5,1.5),c(rep(0,12),5,5,2.5,2.5,1.5,1.5)*3),
                       input_DisabilityWeight_interval<-data.frame(DisabilityWeight/2,DisabilityWeight+DisabilityWeight/2),
                       Input_Incident_Rate=FALSE,
                       year_range=seq(1958,1997,5),
                       
                       Output_YLD_noUI=TRUE)
IntoYLD_noUI$YLD_byyear




save.image(paste0(function_path,'/SLE_related_function','/testing_convergence_meanANDmedian3.RData'))










