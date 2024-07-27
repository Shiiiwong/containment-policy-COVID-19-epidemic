library(dplyr)
library(survival)
library(survminer)
library(ipw)
library(boot)
library(progress)
library(survRM2)
library(forestplot)
setwd()
source("help_func_demo.R")

#Loading country-specific Data
covdata<-read.csv()#Country-specific epidemic data from JHU
covdata$date<-as.Date(covdata$date)

OXCGRT=read.csv() #NPI panel data from OxCGRT
NPI=data.frame("country"=OXCGRT$country,"country_name"=OXCGRT$country_name,"date"=OXCGRT$date,"PolicyValue_C1"=OXCGRT$C1,"PolicyValue_C2"=OXCGRT$C2,"PolicyValue_C3"=OXCGRT$C3,"PolicyValue_C4"=OXCGRT$C4,"PolicyValue_C5"=OXCGRT$C5,"PolicyValue_C6"=OXCGRT$C6,"PolicyValue_C7"=OXCGRT$C7,"PolicyValue_C8"=OXCGRT$C8)
dur=read.csv() #Country-specific data of cumulative duration of each NPI

#time-invariant covariates
L<-read.csv()
#time-varying 
stringency_index=read.csv()
stringency_index=data.frame("country"=stringency_index$country,"date"=stringency_index$date,"stringency_index_C1"=stringency_index$C1,"stringency_index_C2"=stringency_index$C2,"stringency_index_C3"=stringency_index$C3,"stringency_index_C4"=stringency_index$C4,"stringency_index_C5"=stringency_index$C5,"stringency_index_C6"=stringency_index$C6,"stringency_index_C7"=stringency_index$C7,"stringency_index_C8"=stringency_index$C8)
total_vaccinations=read.csv()

#Analysis Data preparation
studybegin=as.Date("2020/1/1")
ana.data_5=anadata.fun(studyend="2022/12/31","cases",targetcase=100000)
ana.data_4=anadata.fun(studyend="2022/12/31","cases",targetcase=10000)
ana.data_3=anadata.fun(studyend="2022/12/31","cases",targetcase=1000)
###########################################################################
##############        Step I: Within NPI comparison         ##############
###########################################################################
#C1:School closure
incidence.res_C1=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C1","incidence")),
                       data.frame("endpoint"="10^4",step1.func(ana.data_4,"C1","incidence")),
                       data.frame("endpoint"="10^3",step1.func(ana.data_3,"C1","incidence")))                      #incidence
best_C1L=rbind(incidence.res_C1%>%filter(policy=="C1L1")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C1%>%filter(policy=="C1L2")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C1%>%filter(policy=="C1L3")%>%group_by(endpoint)%>%filter(incidence==min(incidence))) #Best maintaining strategy by stringency level
diff.res_C1L=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C1","diff")),
                   data.frame("endpoint"="10^4",step1.func(ana.data_4,"C1","diff")),
                   data.frame("endpoint"="10^3",step1.func(ana.data_3,"C1","diff")))                            #incidence difference (little longer due to bootstrap CI)

#C2:Workplace closure
incidence.res_C2=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C2","incidence")),
                       data.frame("endpoint"="10^4",step1.func(ana.data_4,"C2","incidence")),
                       data.frame("endpoint"="10^3",step1.func(ana.data_3,"C2","incidence")))
best_C2L=rbind(incidence.res_C2%>%filter(policy=="C2L1")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C2%>%filter(policy=="C2L2")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C2%>%filter(policy=="C2L3")%>%group_by(endpoint)%>%filter(incidence==min(incidence)))
diff.res_C2L=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C2","diff")),
                   data.frame("endpoint"="10^4",step1.func(ana.data_4,"C2","diff")),
                   data.frame("endpoint"="10^3",step1.func(ana.data_3,"C2","diff")))   
#C3:Cancel public events
incidence.res_C3=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C3","incidence")),
                       data.frame("endpoint"="10^4",step1.func(ana.data_4,"C3","incidence")),
                       data.frame("endpoint"="10^3",step1.func(ana.data_3,"C3","incidence")))
best_C3L=rbind(incidence.res_C3%>%filter(policy=="C3L1")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C3%>%filter(policy=="C3L2")%>%group_by(endpoint)%>%filter(incidence==min(incidence)))
diff.res_C3L=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C3","diff")),
                   data.frame("endpoint"="10^4",step1.func(ana.data_4,"C3","diff")),
                   data.frame("endpoint"="10^3",step1.func(ana.data_3,"C3","diff"))) 
#C4:Restrictions on gathering
incidence.res_C4=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C4","incidence")),
                       data.frame("endpoint"="10^4",step1.func(ana.data_4,"C4","incidence")),
                       data.frame("endpoint"="10^3",step1.func(ana.data_3,"C4","incidence")))
best_C4L=rbind(incidence.res_C4%>%filter(policy=="C4L1")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C4%>%filter(policy=="C4L2")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C4%>%filter(policy=="C4L3")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C4%>%filter(policy=="C4L4")%>%group_by(endpoint)%>%filter(incidence==min(incidence)))
diff.res_C4L=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C4","diff")),
                   data.frame("endpoint"="10^4",step1.func(ana.data_4,"C4","diff")),
                   data.frame("endpoint"="10^3",step1.func(ana.data_3,"C4","diff"))) 
#C5:Closing public transport
incidence.res_C5=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C5","incidence")),
                       data.frame("endpoint"="10^4",step1.func(ana.data_4,"C5","incidence")),
                       data.frame("endpoint"="10^3",step1.func(ana.data_3,"C5","incidence")))
best_C5L=rbind(incidence.res_C5%>%filter(policy=="C5L1")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C5%>%filter(policy=="C5L2")%>%group_by(endpoint)%>%filter(incidence==min(incidence)))
diff.res_C5L=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C5","diff")),
                   data.frame("endpoint"="10^4",step1.func(ana.data_4,"C5","diff")),
                   data.frame("endpoint"="10^3",step1.func(ana.data_3,"C5","diff"))) 
#C6:Stay at home requirements
incidence.res_C6=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C6","incidence")),
                       data.frame("endpoint"="10^4",step1.func(ana.data_4,"C6","incidence")),
                       data.frame("endpoint"="10^3",step1.func(ana.data_3,"C6","incidence")))
best_C6L=rbind(incidence.res_C6%>%filter(policy=="C6L1")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C6%>%filter(policy=="C6L2")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C6%>%filter(policy=="C6L3")%>%group_by(endpoint)%>%filter(incidence==min(incidence)))
diff.res_C6L=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C6","diff")),
                   data.frame("endpoint"="10^4",step1.func(ana.data_4,"C6","diff")),
                   data.frame("endpoint"="10^3",step1.func(ana.data_3,"C6","diff")))  
#C7:Internal movement restrictions
incidence.res_C7=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C7","incidence")),
                       data.frame("endpoint"="10^4",step1.func(ana.data_4,"C7","incidence")),
                       data.frame("endpoint"="10^3",step1.func(ana.data_3,"C7","incidence")))
best_C7L=rbind(incidence.res_C7%>%filter(policy=="C7L1")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C7%>%filter(policy=="C7L2")%>%group_by(endpoint)%>%filter(incidence==min(incidence)))
diff.res_C7L=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C7","diff")),
                   data.frame("endpoint"="10^4",step1.func(ana.data_4,"C7","diff")),
                   data.frame("endpoint"="10^3",step1.func(ana.data_3,"C7","diff"))) 
#C8:International movement restrictions
incidence.res_C8=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C8","incidence")),
                       data.frame("endpoint"="10^4",step1.func(ana.data_4,"C8","incidence")),
                       data.frame("endpoint"="10^3",step1.func(ana.data_3,"C8","incidence")))
best_C8L=rbind(incidence.res_C8%>%filter(policy=="C8L1")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C8%>%filter(policy=="C8L2")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C8%>%filter(policy=="C8L3")%>%group_by(endpoint)%>%filter(incidence==min(incidence)),
               incidence.res_C8%>%filter(policy=="C8L4")%>%group_by(endpoint)%>%filter(incidence==min(incidence)))
diff.res_C8L=rbind(data.frame("endpoint"="10^5",step1.func(ana.data_5,"C8","diff")),
                   data.frame("endpoint"="10^4",step1.func(ana.data_4,"C8","diff")),
                   data.frame("endpoint"="10^3",step1.func(ana.data_3,"C8","diff"))) 

##########################################################################
##############       Step II: Between NPI comparison       ##############
##########################################################################
#C1
best_C1=best_C1L%>%group_by(endpoint)%>%filter(incidence==min(incidence))#Best stringency level and maintaining strategy
diff.res_C1=rbind(data.frame("endpoint"="10^5",step2.func(ana.data_5,"C1","diff")),
                  data.frame("endpoint"="10^4",step2.func(ana.data_4,"C1","diff")),
                  data.frame("endpoint"="10^3",step2.func(ana.data_3,"C1","diff")))
#C2
best_C2=best_C2L%>%group_by(endpoint)%>%filter(incidence==min(incidence))
diff.res_C2=rbind(data.frame("endpoint"="10^5",step2.func(ana.data_5,"C2","diff")),
                  data.frame("endpoint"="10^4",step2.func(ana.data_4,"C2","diff")),
                  data.frame("endpoint"="10^3",step2.func(ana.data_3,"C2","diff")))
#C3
best_C3=best_C3L%>%group_by(endpoint)%>%filter(incidence==min(incidence))
diff.res_C3=rbind(data.frame("endpoint"="10^5",step2.func(ana.data_5,"C3","diff")),
                  data.frame("endpoint"="10^4",step2.func(ana.data_4,"C3","diff")),
                  data.frame("endpoint"="10^3",step2.func(ana.data_3,"C3","diff")))
#C4
best_C4=best_C4L%>%group_by(endpoint)%>%filter(incidence==min(incidence))
diff.res_C4=rbind(data.frame("endpoint"="10^5",step2.func(ana.data_5,"C4","diff")),
                  data.frame("endpoint"="10^4",step2.func(ana.data_4,"C4","diff")),
                  data.frame("endpoint"="10^3",step2.func(ana.data_3,"C4","diff")))
#C5
best_C5=best_C5L%>%group_by(endpoint)%>%filter(incidence==min(incidence))
diff.res_C5=rbind(data.frame("endpoint"="10^5",step2.func(ana.data_5,"C5","diff")),
                  data.frame("endpoint"="10^4",step2.func(ana.data_4,"C5","diff")),
                  data.frame("endpoint"="10^3",step2.func(ana.data_3,"C5","diff")))
#C6
best_C6=best_C6L%>%group_by(endpoint)%>%filter(incidence==min(incidence))
diff.res_C6=rbind(data.frame("endpoint"="10^5",step2.func(ana.data_5,"C6","diff")),
                  data.frame("endpoint"="10^4",step2.func(ana.data_4,"C6","diff")),
                  data.frame("endpoint"="10^3",step2.func(ana.data_3,"C6","diff")))
#C7
best_C7=best_C7L%>%group_by(endpoint)%>%filter(incidence==min(incidence))
diff.res_C7=rbind(data.frame("endpoint"="10^5",step2.func(ana.data_5,"C7","diff")),
                  data.frame("endpoint"="10^4",step2.func(ana.data_4,"C7","diff")),
                  data.frame("endpoint"="10^3",step2.func(ana.data_3,"C7","diff")))
#C8
best_C8=best_C8L%>%group_by(endpoint)%>%filter(incidence==min(incidence))
diff.res_C8=rbind(data.frame("endpoint"="10^5",step2.func(ana.data_5,"C8","diff")),
                  data.frame("endpoint"="10^4",step2.func(ana.data_4,"C8","diff")),
                  data.frame("endpoint"="10^3",step2.func(ana.data_3,"C8","diff")))
diff.res.C=rbind(diff.res_C1,diff.res_C2,diff.res_C3,diff.res_C4,diff.res_C5,diff.res_C6,diff.res_C7,diff.res_C8)

##########################################################################
############      Stage III: Between category comparison      ############
##########################################################################
best_C=rbind(best_C1,best_C2,best_C3,best_C4,best_C5,best_C6,best_C7,best_C8)

diff.res_3=data.frame("endpoint"="10^3",step3.func(ana.data_3))
diff.res_4=data.frame("endpoint"="10^4",step3.func(ana.data_4))
diff.res_5=data.frame("endpoint"="10^5",step3.func(ana.data_5))


