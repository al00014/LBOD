YLD_sub_function<-function(data_plugin, ##  YLD_sub_function is only for the calculation of incident YLD
                           Rate,
                           Beta,
                           Const,
                           Agewt,
                           BplusR){
  options(scipen = 999)
  #if(sum(!c('parallel') %in% rownames(installed.packages()) )!= 0) {
  #  stop('\n
  #       Please installed the dependencies first:\n
  #       "parallel" ')
  #}
  names(data_plugin)<-c('Age','population','incidence','age_at_onset','duration','DisabilityWeight')
 
  #require(parallel)
  if (Rate==0){
    YLD<-unlist(parallel::mclapply(1:nrow(data_plugin),FUN=function(i){(data_plugin$incidence[i]*data_plugin$DisabilityWeight[i])*
        (Agewt*Const*
           ((exp(-Beta*data_plugin$age_at_onset[i]))/Beta^2)*
           ((exp(-Beta*data_plugin$duration[i]))*
              (-Beta*(data_plugin$duration[i]+data_plugin$age_at_onset[i])-1)-
              (-Beta*data_plugin$age_at_onset[i]-1))+
           ((1-Agewt)*data_plugin$duration[i]))}))
  }else if (Rate>0){
    YLD<-unlist(parallel::mclapply(1:nrow(data_plugin),FUN=function(i){data_plugin$incidence[i]*data_plugin$DisabilityWeight[i]*
        (Agewt*((Const*exp(Rate*data_plugin$age_at_onset[i]))/(BplusR^2))*
           ((exp(BplusR*(data_plugin$duration[i]+data_plugin$age_at_onset[i]))*
               (BplusR*(data_plugin$duration[i]+data_plugin$age_at_onset[i])-1))-
              (exp(BplusR*data_plugin$age_at_onset[i])*
                 (BplusR*data_plugin$age_at_onset[i]-1)))+
           ((1-Agewt)/Rate)*
           ((1-exp(-Rate*data_plugin$duration[i]))))}))
  }
  return(YLD)
}