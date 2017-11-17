YLL_sub_function<-function(data_plugin,
                           Stand_life_expectancy,
                           Rate,
                           Beta,
                           Const,
                           Agewt,
                           BplusR){
  #if(sum(!c('parallel') %in% rownames(installed.packages()) )!= 0) {
  #  stop('\n
  #        Please installed the dependencies first:\n
  #       "parallel" ')
  #}
  
  #require(parallel)
  if (Rate==0){
    YLL<-unlist(parallel::mclapply(1:nrow(Stand_life_expectancy), 
                         function(i){(data_plugin$death[i]*1)*(Agewt*Const*((exp(-Beta*data_plugin$age_av_death[i]))/Beta^2)*
                                                                 ((exp(-Beta*Stand_life_expectancy$LE[i]))*
                                                                    (-Beta*(Stand_life_expectancy$LE[i]+data_plugin$age_av_death[i])-1)-
                                                                    (-Beta*data_plugin$age_av_death[i]-1))+
                                                                 ((1-Agewt)*Stand_life_expectancy$LE[i]))}))
  } else if(Rate>0){
    YLL<-unlist(parallel::mclapply(1:nrow(Stand_life_expectancy),
                         function(i){data_plugin$death[i]*1*(Agewt*((Const*exp(Rate*data_plugin$age_av_death[i]))/(BplusR^2))*
                                                               ((exp(BplusR*(Stand_life_expectancy$LE[i]+data_plugin$age_av_death[i]))*
                                                                   (BplusR*(Stand_life_expectancy$LE[i]+data_plugin$age_av_death[i])-1))-
                                                                  (exp(BplusR*data_plugin$age_av_death[i])*
                                                                     (BplusR*data_plugin$age_av_death[i]-1)))+((1-Agewt)/Rate)*
                                                               ((1-exp(-Rate*Stand_life_expectancy$LE[i]))))}))
  }
  return(YLL)
}