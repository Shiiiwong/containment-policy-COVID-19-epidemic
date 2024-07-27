anadata.fun <- function(studyend,event,targetcase){
  
  studyend=as.Date(studyend)
  raw.covdata<-subset(covdata,date<=studyend)
  if(event=="cases"){
    for (i in 1:length(raw.covdata$total_cases_per_million)) {
      if(is.na(raw.covdata$total_cases_per_million[i])!=1&&raw.covdata$total_cases_per_million[i]>=targetcase){raw.covdata$event[i]<-1}else{raw.covdata$event[i]<-0}
    }
  }else{}
  if(event=="deaths"){
    for (i in 1:length(raw.covdata$total_deaths_per_million)) {
      if(is.na(raw.covdata$total_deaths_per_million[i])!=1&&raw.covdata$total_deaths_per_million[i]>=targetcase){raw.covdata$event[i]<-1}else{raw.covdata$event[i]<-0}
    }
  }else{} 
  raw.covdata<-rbind(subset(raw.covdata,event==0),subset(raw.covdata,event==1)%>%distinct(country,.keep_all=T))
  NPI$date<-as.Date(0:1095,origin="2020-01-01")
  raw.data<-full_join(NPI,raw.covdata,by=c("country"="country","date"="date"))
  raw.data$event[which(is.na(raw.data$event)==1)]=0
  raw.data<-raw.data[order(raw.data$country,raw.data$date),]
  
  ana.data=data.frame("country"=sort(rep(unique(L$country),length(as.Date(studybegin:studyend)))),"id"=sort(rep(1:length(unique(L$country)),length(as.Date(studybegin:studyend)))),"tstartdate"=as.Date(studybegin:studyend))
  ana.data$tstart<-ana.data$tstartdate-studybegin-1
  ana.data$fuptime<-ana.data$tstart+1
  ana.data<-left_join(ana.data,raw.data,by=c("country"="country","tstartdate"="date"))
  
  #covariate
  ##time-invariant
  pb <- progress_bar$new(format="Analysis data preparing[:bar]:percent Time required::eta",total = length(L$country))
  ana.data$population_density=0
  ana.data$GDPpercapita=0
  ana.data$globalization=0
  ana.data$bed=0
  ana.data$aged_over65=0
  
  for (j in 1:length(L$country)) {
    for(i in 1:length(ana.data$country)){
      if(ana.data$country[i]==L$country[j]){ana.data$population_density[i]=L$population_density[j]}else{}
      if(ana.data$country[i]==L$country[j]){ana.data$GDPpercapita[i]=L$GDP_per_capita[j]}else{}
      if(ana.data$country[i]==L$country[j]){ana.data$globalization[i]=L$globalization_index[j]}else{}
      if(ana.data$country[i]==L$country[j]){ana.data$bed[i]=L$bed[j]}else{}
      if(ana.data$country[i]==L$country[j]){ana.data$aged_over65[i]=L$aged_over65[j]}else{}}
    pb$tick()
  }
  
  ##time-varying
  #total_vaccinations
  total_vaccinations$date=as.Date(total_vaccinations$date)
  ana.data<-left_join(ana.data,total_vaccinations,by=c("country"="country","tstartdate"="date"))
  ana.data$total_vaccinations[which(is.na(ana.data$total_vaccinations)==1)]=0
  #Overall NPI index
  stringency_index$date<-as.Date(0:1095,origin="2020-01-01")
  ana.data<-left_join(ana.data,stringency_index,by=c("country"="country","tstartdate"="date"))
  for(i in 1:length(ana.data$country)){
    if(i+1<=length(ana.data$country)&&ana.data$country[i+1]==ana.data$country[i]&&ana.data$event[i]==1){ana.data$event[i+1]=1}else{}
  }
  ana.data<-rbind(subset(ana.data,event==0),subset(ana.data,event==1)%>%distinct(country,.keep_all=T))
  ana.data$policy=ana.data$PolicyValue
  
  return(ana.data)
}

step1.func=function(ana.data,category,results_type){
  
  if(category=="C1"){
    ana.data$policy=ana.data$PolicyValue_C1
    ana.data$stringency_index=ana.data$stringency_index_C1
  }else{}
  
  if(category=="C2"){
    ana.data$policy=ana.data$PolicyValue_C2
    ana.data$stringency_index=ana.data$stringency_index_C2
  }else{}
  
  if(category=="C3"){
    ana.data$policy=ana.data$PolicyValue_C3
    ana.data$stringency_index=ana.data$stringency_index_C3
  }else{}
  
  if(category=="C4"){
    ana.data$policy=ana.data$PolicyValue_C4
    ana.data$stringency_index=ana.data$stringency_index_C4
  }else{}
  
  if(category=="C5"){
    ana.data$policy=ana.data$PolicyValue_C5
    ana.data$stringency_index=ana.data$stringency_index_C5
  }else{}
  
  if(category=="C6"){
    ana.data$policy=ana.data$PolicyValue_C6
    ana.data$stringency_index=ana.data$stringency_index_C6
  }else{}
  
  if(category=="C7"){
    ana.data$policy=ana.data$PolicyValue_C7
    ana.data$stringency_index=ana.data$stringency_index_C7
  }else{}
  
  if(category=="C8"){
    ana.data$policy=ana.data$PolicyValue_C8
    ana.data$stringency_index=ana.data$stringency_index_C8
  }else{}
  
  #1.CCW_func
  CCW.func=function(ana.data,policy,results_type){
    
    target_col=which(names(dur)==policy)
    Q1=quantile(dur[,target_col],0.25)
    Q2=quantile(dur[,target_col],0.50)
    Q3=quantile(dur[,target_col],0.75)
    
    #1.1 cloning
    
    copy_1 <- ana.data %>% mutate(cgroup=1) 
    copy_2 <- ana.data %>% mutate(cgroup=2)
    copy_3 <- ana.data %>% mutate(cgroup=3) 
    copy_4 <- ana.data %>% mutate(cgroup=4) 
    
    #1.2 censoring
    if(policy=="C1L1"|policy=="C2L1"|policy=="C3L1"|policy=="C4L1"|policy=="C5L1"|policy=="C6L1"|policy=="C7L1"|policy=="C8L1"){policylevel=1}else{}
    if(policy=="C1L2"|policy=="C2L2"|policy=="C3L2"|policy=="C4L2"|policy=="C5L2"|policy=="C6L2"|policy=="C7L2"|policy=="C8L2"){policylevel=2}else{}
    if(policy=="C1L3"|policy=="C2L3"|policy=="C4L3"|policy=="C6L3"|policy=="C8L3"){policylevel=3}else{}
    if(policy=="C4L4"|policy=="C8L4"){policylevel=4}else{}
    
    copy_1<-copy_1[order(copy_1$id,copy_1$tstartdate),]
    temp_1 <- copy_1 %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
    temp_2 <- copy_1 %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q1, 1, 0))
    copy_1 <- bind_rows(temp_1, temp_2)
    copy_1<-copy_1[order(copy_1$id,copy_1$tstartdate),]
    for (i in 1:length(copy_1$country)) {
      if(i+1<=length(copy_1$country)&&copy_1$country[i+1]==copy_1$country[i]&&copy_1$artificial_censor[i]==1){copy_1$artificial_censor[i+1]=1}else{}}
    copy_1<-rbind(subset(copy_1,artificial_censor==0),subset(copy_1,artificial_censor==1)%>%distinct(country,.keep_all=T))
    copy_1<-copy_1[order(copy_1$id,copy_1$tstartdate),]
    
    copy_2<-copy_2[order(copy_2$id,copy_2$tstartdate),]
    temp_1 <- copy_2 %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
    temp_2 <- copy_2 %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q1&row_number()==max(row_number()), 1, 0))
    temp_3 <- copy_2 %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q2, 1, 0))
    copy_2 <- bind_rows(temp_1, temp_2,temp_3)
    copy_2<-copy_2[order(copy_2$id,copy_2$tstartdate),]
    for (i in 1:length(copy_2$country)) {
      if(i+1<=length(copy_2$country)&&copy_2$country[i+1]==copy_2$country[i]&&copy_2$artificial_censor[i]==1){copy_2$artificial_censor[i+1]=1}else{}}
    copy_2<-rbind(subset(copy_2,artificial_censor==0),subset(copy_2,artificial_censor==1)%>%distinct(country,.keep_all=T))
    copy_2<-copy_2[order(copy_2$id,copy_2$tstartdate),]
    
    copy_3<-copy_3[order(copy_3$id,copy_3$tstartdate),]
    temp_1 <- copy_3 %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
    temp_2 <- copy_3 %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q2&row_number()==max(row_number()), 1, 0))
    temp_3 <- copy_3 %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q3, 1, 0))
    copy_3 <- bind_rows(temp_1, temp_2,temp_3)
    copy_3<-copy_3[order(copy_3$id,copy_3$tstartdate),]
    for (i in 1:length(copy_3$country)) {
      if(i+1<=length(copy_3$country)&&copy_3$country[i+1]==copy_3$country[i]&&copy_3$artificial_censor[i]==1){copy_3$artificial_censor[i+1]=1}else{}}
    copy_3<-rbind(subset(copy_3,artificial_censor==0),subset(copy_3,artificial_censor==1)%>%distinct(country,.keep_all=T))
    copy_3<-copy_3[order(copy_3$id,copy_3$tstartdate),]
    
    copy_4<-copy_4[order(copy_4$id,copy_4$tstartdate),]
    temp_1 <- copy_4 %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
    temp_2 <- copy_4 %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q3&row_number()==max(row_number()), 1, 0))
    copy_4 <- bind_rows(temp_1, temp_2)
    copy_4<-copy_4[order(copy_4$id,copy_4$tstartdate),]
    for (i in 1:length(copy_4$country)) {
      if(i+1<=length(copy_4$country)&&copy_4$country[i+1]==copy_4$country[i]&&copy_4$artificial_censor[i]==1){copy_4$artificial_censor[i+1]=1}else{}}
    copy_4<-rbind(subset(copy_4,artificial_censor==0),subset(copy_4,artificial_censor==1)%>%distinct(country,.keep_all=T))
    copy_4<-copy_4[order(copy_4$id,copy_4$tstartdate),]
    
    copy=rbind(copy_1,copy_2,copy_3,copy_4)
    
    #1.2.1 BOOT
    boot.fun.1 <- function(dat, index){
      
      copy_1<-copy[copy$cgroup==1,]
      copy_2<-copy[copy$cgroup==2,]
      copy_3<-copy[copy$cgroup==3,]
      copy_4<-copy[copy$cgroup==4,]
      copy_ref <- copy_1[index,] # allows boot to select sample
      select<-copy_ref$id
      copy_2boot<-copy_2[copy_2$id %in% select,] 
      copy_3boot<-copy_3[copy_3$id %in% select,]
      copy_4boot<-copy_4[copy_4$id %in% select,]
      copy<-rbind(copy_ref,copy_2boot,copy_3boot,copy_4boot)
      
      #weight
      copy<-data.frame("id"=copy$id,"tstart"=copy$tstart,"fuptime"=copy$fuptime,"artificial_censor"=copy$artificial_censor,"policy"=copy$policy,"event"=copy$event,"cgroup"=copy$cgroup,"total_cases_per_million"=copy$total_cases_per_million,"population_density"=copy$population_density,"bed"=copy$bed,"GDP"=copy$GDPpercapita,"stringency_index"=copy$stringency_index,"glob"=copy$globalization,"total_vaccinations"=copy$total_vaccinations,"aged_over65"=copy$aged_over65)
      copy=subset(copy,is.na(copy$stringency_index)==0&is.na(copy$GDP)==0&is.na(copy$bed)==0)
      
      weight_1 <- ipwtm(
        exposure = artificial_censor,
        family = "survival",
        numerator = ~ population_density+GDP+bed+aged_over65+glob,
        denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
        id = id,
        tstart = tstart,
        timevar = fuptime,
        type = "cens",
        trunc = 0.01,
        data = copy[copy$cgroup==1,])
      weight_2 <- ipwtm(
        exposure = artificial_censor,
        family = "survival",
        numerator = ~ population_density+GDP+bed+aged_over65+glob,
        denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
        id = id,
        tstart = tstart,
        timevar = fuptime,
        type = "cens",
        trunc = 0.01,
        data = copy[copy$cgroup==2,])
      weight_3 <- ipwtm(
        exposure = artificial_censor,
        family = "survival",
        numerator = ~ population_density+GDP+bed+aged_over65+glob,
        denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
        id = id,
        tstart = tstart,
        timevar = fuptime,
        type = "cens",
        trunc = 0.01,
        data = copy[copy$cgroup==3,])
      weight_4 <- ipwtm(
        exposure = artificial_censor,
        family = "survival",
        numerator = ~ population_density+GDP+bed+aged_over65+glob,
        denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
        id = id,
        tstart = tstart,
        timevar = fuptime,
        type = "cens",
        trunc = 0.01,
        data = copy[copy$cgroup==4,])
      
      
      #fit
      fita<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                    data=copy[copy$cgroup==1,],weights =weight_1$ipw.weights)
      fitb<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                    data=copy[copy$cgroup==2,],weights =weight_2$ipw.weights)
      fitc<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                    data=copy[copy$cgroup==3,],weights =weight_3$ipw.weights)
      fitd<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                    data=copy[copy$cgroup==4,],weights =weight_4$ipw.weights)
      
      diff.res=c(min(fita$surv)-min(fitb$surv),min(fita$surv)-min(fitc$surv),min(fita$surv)-min(fitd$surv))
      names(diff.res)=c("S1-S2","S1-S3","S1-S4")
      return(diff.res)
    }
    
    if(results_type=="diff"){bootres=boot(copy, boot.fun.1, 200)}else{} 
    
    #copy dataset
    copy<-data.frame("id"=copy$id,"tstart"=copy$tstart,"fuptime"=copy$fuptime,"artificial_censor"=copy$artificial_censor,"policy"=copy$policy,"event"=copy$event,"cgroup"=copy$cgroup,"total_cases_per_million"=copy$total_cases_per_million,"population_density"=copy$population_density,"bed"=copy$bed,"GDP"=copy$GDPpercapita,"stringency_index"=copy$stringency_index,"glob"=copy$globalization,"total_vaccinations"=copy$total_vaccinations,"aged_over65"=copy$aged_over65)
    copy=subset(copy,is.na(copy$stringency_index)==0&is.na(copy$GDP)==0&is.na(copy$bed)==0)
    
    #1.3 weight
    
    weight_1 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==1,])
    weight_2 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==2,])
    weight_3 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==3,])
    weight_4 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==4,])
    
    weights <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy)
    
    fita<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==1,],weights =weight_1$ipw.weights)
    fitb<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==2,],weights =weight_2$ipw.weights)
    fitc<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==3,],weights =weight_3$ipw.weights)
    fitd<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==4,],weights =weight_4$ipw.weights)
    
    if(results_type=="diff"){
      diff.res=data.frame("diff"=as.numeric(list(0,min(fita$surv)-min(fitb$surv),min(fita$surv)-min(fitc$surv),min(fita$surv)-min(fitd$surv))),"lower"=0,"upper"=0)
      
      diff.res[2,2:3]=boot.ci(bootres,conf = 0.95, type = "norm",index=1)$normal[1,2:3]
      diff.res[3,2:3]=boot.ci(bootres,conf = 0.95, type = "norm",index=2)$normal[1,2:3]
      diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
      
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence=data.frame("incidence"=as.numeric(list(1-min(fita$surv),1-min(fitb$surv),1-min(fitc$surv),1-min(fitd$surv))))
      return(incidence)
    }else{}
    
    if(results_type=="ptrend"){
      glm_ptrend=glm(event ~ cgroup+population_density+GDP+bed+aged_over65+glob+stringency_index, data = copy, weights = weights$ipw.weights, family = "quasibinomial")
      ptrend=coeftest(glm_ptrend, vcov=vcovHC(glm_ptrend, type="HC1")) [2,4]
      return(ptrend)
    }else{}
    
    if(results_type=="HR"){
      glm_HR <- glm(event ~ factor(cgroup)+population_density+GDP+bed+aged_over65+glob+stringency_index , data = copy, weights = weights$ipw.weights, family = "quasibinomial")
      HR=data.frame("est"=c(1,exp(coef(glm_HR))[2:4]),"lower"=c(1,exp(coefci(glm_HR, vcov=vcovHC(glm_HR, type="HC1")))[2:4,1]),"upper"=c(1,exp(coefci(glm_HR, vcov=vcovHC(glm_HR, type="HC1")))[2:4,2]))
      return(HR)
    }else{}
  }
  
  #2.output
  
  if(category=="C1"){
    if(results_type=="diff"){
      diff.res=rbind(data.frame("policy"=rep("C1L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C1L1","diff")),
                     data.frame("policy"=rep("C1L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C1L2","diff")),
                     data.frame("policy"=rep("C1L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C1L3","diff")))
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence.res=rbind(data.frame("policy"=rep("C1L1",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C1L1","incidence")),
                          data.frame("policy"=rep("C1L2",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C1L2","incidence")),
                          data.frame("policy"=rep("C1L3",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C1L3","incidence")))
      return(incidence.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=rbind(data.frame("policy"="C1L1","ptrend"=CCW.func(ana.data,"C1L1","ptrend")),
                       data.frame("policy"="C1L2","ptrend"=CCW.func(ana.data,"C1L2","ptrend")),
                       data.frame("policy"="C1L3","ptrend"=CCW.func(ana.data,"C1L3","ptrend")))
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=rbind(data.frame("policy"=rep("C1L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C1L1","HR")),
                   data.frame("policy"=rep("C1L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C1L2","HR")),
                   data.frame("policy"=rep("C1L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C1L3","HR")))
      return(HR.res)
    }else{}
  }else{}
  
  if(category=="C2"){
    if(results_type=="diff"){
      diff.res=rbind(data.frame("policy"=rep("C2L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C2L1","diff")),
                     data.frame("policy"=rep("C2L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C2L2","diff")),
                     data.frame("policy"=rep("C2L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C2L3","diff")))
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence.res=rbind(data.frame("policy"=rep("C2L1",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C2L1","incidence")),
                          data.frame("policy"=rep("C2L2",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C2L2","incidence")),
                          data.frame("policy"=rep("C2L3",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C2L3","incidence")))
      return(incidence.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=rbind(data.frame("policy"="C2L1","ptrend"=CCW.func(ana.data,"C2L1","ptrend")),
                       data.frame("policy"="C2L2","ptrend"=CCW.func(ana.data,"C2L2","ptrend")),
                       data.frame("policy"="C2L3","ptrend"=CCW.func(ana.data,"C2L3","ptrend")))
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=rbind(data.frame("policy"=rep("C2L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C2L1","HR")),
                   data.frame("policy"=rep("C2L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C2L2","HR")),
                   data.frame("policy"=rep("C2L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C2L3","HR")))
      return(HR.res)
    }else{}
  }else{}
  
  if(category=="C3"){
    if(results_type=="diff"){
      diff.res=rbind(data.frame("policy"=rep("C3L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C3L1","diff")),
                     data.frame("policy"=rep("C3L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C3L2","diff")))
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence.res=rbind(data.frame("policy"=rep("C3L1",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C3L1","incidence")),
                          data.frame("policy"=rep("C3L2",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C3L2","incidence")))
      return(incidence.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=rbind(data.frame("policy"="C3L1","ptrend"=CCW.func(ana.data,"C3L1","ptrend")),
                       data.frame("policy"="C3L2","ptrend"=CCW.func(ana.data,"C3L2","ptrend")))
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=rbind(data.frame("policy"=rep("C3L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C3L1","HR")),
                   data.frame("policy"=rep("C3L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C3L2","HR")))
      return(HR.res)
    }else{}
  }else{}
  
  if(category=="C4"){
    if(results_type=="diff"){
      diff.res=rbind(data.frame("policy"=rep("C4L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C4L1","diff")),
                     data.frame("policy"=rep("C4L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C4L2","diff")) ,
                     data.frame("policy"=rep("C4L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C4L3","diff")),
                     data.frame("policy"=rep("C4L4",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C4L4","diff")))
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence.res=rbind(data.frame("policy"=rep("C4L1",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C4L1","incidence")),
                          data.frame("policy"=rep("C4L2",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C4L2","incidence")),
                          data.frame("policy"=rep("C4L3",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C4L3","incidence")),
                          data.frame("policy"=rep("C4L4",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C4L4","incidence")))
      return(incidence.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=rbind(data.frame("policy"="C4L1","ptrend"=CCW.func(ana.data,"C4L1","ptrend")),
                       data.frame("policy"="C4L2","ptrend"=CCW.func(ana.data,"C4L2","ptrend")),
                       data.frame("policy"="C4L3","ptrend"=CCW.func(ana.data,"C4L3","ptrend")),
                       data.frame("policy"="C4L4","ptrend"=CCW.func(ana.data,"C4L4","ptrend")))
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=rbind(data.frame("policy"=rep("C4L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C4L1","HR")),
                   data.frame("policy"=rep("C4L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C4L2","HR")),
                   data.frame("policy"=rep("C4L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C4L3","HR")),
                   data.frame("policy"=rep("C4L4",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C4L4","HR")))
      return(HR.res)
    }else{}
    
  }else{}
  
  if(category=="C5"){
    if(results_type=="diff"){
      diff.res=rbind(data.frame("policy"=rep("C5L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C5L1","diff")),
                     data.frame("policy"=rep("C5L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C5L2","diff")))
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence.res=rbind(data.frame("policy"=rep("C5L1",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C5L1","incidence")),
                          data.frame("policy"=rep("C5L2",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C5L2","incidence")))
      return(incidence.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=rbind(data.frame("policy"="C5L1","ptrend"=CCW.func(ana.data,"C5L1","ptrend")),
                       data.frame("policy"="C5L2","ptrend"=CCW.func(ana.data,"C5L2","ptrend")))
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=rbind(data.frame("policy"=rep("C5L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C5L1","HR")),
                   data.frame("policy"=rep("C5L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C5L2","HR")))
      return(HR.res)
    }else{}
  }else{}
  
  if(category=="C6"){
    if(results_type=="diff"){
      diff.res=rbind(data.frame("policy"=rep("C6L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C6L1","diff")),
                     data.frame("policy"=rep("C6L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C6L2","diff")),
                     data.frame("policy"=rep("C6L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C6L3","diff")))
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence.res=rbind(data.frame("policy"=rep("C6L1",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C6L1","incidence")),
                          data.frame("policy"=rep("C6L2",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C6L2","incidence")),
                          data.frame("policy"=rep("C6L3",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C6L3","incidence")))
      return(incidence.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=rbind(data.frame("policy"="C6L1","ptrend"=CCW.func(ana.data,"C6L1","ptrend")),
                       data.frame("policy"="C6L2","ptrend"=CCW.func(ana.data,"C6L2","ptrend")),
                       data.frame("policy"="C6L3","ptrend"=CCW.func(ana.data,"C6L3","ptrend")))
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=rbind(data.frame("policy"=rep("C6L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C6L1","HR")),
                   data.frame("policy"=rep("C6L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C6L2","HR")),
                   data.frame("policy"=rep("C6L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C6L3","HR")))
      return(HR.res)
    }else{}
  }else{}
  
  if(category=="C7"){
    if(results_type=="diff"){
      diff.res=rbind(data.frame("policy"=rep("C7L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C7L1","diff")),
                     data.frame("policy"=rep("C7L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C7L2","diff")))
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence.res=rbind(data.frame("policy"=rep("C7L1",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C7L1","incidence")),
                          data.frame("policy"=rep("C7L2",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C7L2","incidence")))
      return(incidence.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=rbind(data.frame("policy"="C7L1","ptrend"=CCW.func(ana.data,"C7L1","ptrend")),
                       data.frame("policy"="C7L2","ptrend"=CCW.func(ana.data,"C7L2","ptrend")))
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=rbind(data.frame("policy"=rep("C7L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C7L1","HR")),
                   data.frame("policy"=rep("C7L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C7L2","HR")))
      return(HR.res)
    }else{}
  }else{}
  
  if(category=="C8"){
    if(results_type=="diff"){
      diff.res=rbind(data.frame("policy"=rep("C8L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C8L1","diff")),
                     data.frame("policy"=rep("C8L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C8L2","diff")) ,
                     data.frame("policy"=rep("C8L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C8L3","diff")),
                     data.frame("policy"=rep("C8L4",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),CCW.func(ana.data,"C8L4","diff")))
      return(diff.res)
    }else{}
    
    if(results_type=="incidence"){
      incidence.res=rbind(data.frame("policy"=rep("C8L1",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C8L1","incidence")),
                          data.frame("policy"=rep("C8L2",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C8L2","incidence")),
                          data.frame("policy"=rep("C8L3",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C8L3","incidence")),
                          data.frame("policy"=rep("C8L4",4),"strategy"=c("S1","S2","S3","S4"),"incidence"=CCW.func(ana.data,"C8L4","incidence")))
      return(incidence.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=rbind(data.frame("policy"="C8L1","ptrend"=CCW.func(ana.data,"C8L1","ptrend")),
                       data.frame("policy"="C8L2","ptrend"=CCW.func(ana.data,"C8L2","ptrend")),
                       data.frame("policy"="C8L3","ptrend"=CCW.func(ana.data,"C8L3","ptrend")),
                       data.frame("policy"="C8L4","ptrend"=CCW.func(ana.data,"C8L4","ptrend")))
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=rbind(data.frame("policy"=rep("C8L1",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C8L1","HR")),
                   data.frame("policy"=rep("C8L2",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C8L2","HR")),
                   data.frame("policy"=rep("C8L3",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C8L3","HR")),
                   data.frame("policy"=rep("C8L4",4),"comparison"=c("S1 (ref)","S2-S1","S3-S1","S4-S1"),"HR"=CCW.func(ana.data,"C8L4","HR")))
      return(HR.res)
    }else{}
    
  }else{}
}#step I: incidence, diff, HR and ptrend

step2.func=function(ana.data,category,results_type){
  if(category=="C1"){
    ana.data$policy=ana.data$PolicyValue_C1
    ana.data$stringency_index=ana.data$stringency_index_C1
  }else{}
  
  if(category=="C2"){
    ana.data$policy=ana.data$PolicyValue_C2
    ana.data$stringency_index=ana.data$stringency_index_C2
  }else{}
  
  if(category=="C3"){
    ana.data$policy=ana.data$PolicyValue_C3
    ana.data$stringency_index=ana.data$stringency_index_C3
  }else{}
  
  if(category=="C4"){
    ana.data$policy=ana.data$PolicyValue_C4
    ana.data$stringency_index=ana.data$stringency_index_C4
  }else{}
  
  if(category=="C5"){
    ana.data$policy=ana.data$PolicyValue_C5
    ana.data$stringency_index=ana.data$stringency_index_C5
  }else{}
  
  if(category=="C6"){
    ana.data$policy=ana.data$PolicyValue_C6
    ana.data$stringency_index=ana.data$stringency_index_C6
  }else{}
  
  if(category=="C7"){
    ana.data$policy=ana.data$PolicyValue_C7
    ana.data$stringency_index=ana.data$stringency_index_C7
  }else{}
  
  if(category=="C8"){
    ana.data$policy=ana.data$PolicyValue_C8
    ana.data$stringency_index=ana.data$stringency_index_C8
  }else{}
  
  #1.CCW_func
  #1.1 clone
  copy_L <- ana.data
  if(length(ana.data$country)==length(ana.data_3$country)){endpoint="10^3"}else{}
  if(length(ana.data$country)==length(ana.data_4$country)){endpoint="10^4"}else{}
  if(length(ana.data$country)==length(ana.data_5$country)){endpoint="10^5"}else{}
  
  
  #1.2 censor
  copy.func=function(policy,strategy){
    if(policy=="C1L1"|policy=="C2L1"|policy=="C3L1"|policy=="C4L1"|policy=="C5L1"|policy=="C6L1"|policy=="C7L1"|policy=="C8L1"){policylevel=1}else{}
    if(policy=="C1L2"|policy=="C2L2"|policy=="C3L2"|policy=="C4L2"|policy=="C5L2"|policy=="C6L2"|policy=="C7L2"|policy=="C8L2"){policylevel=2}else{}
    if(policy=="C1L3"|policy=="C2L3"|policy=="C4L3"|policy=="C6L3"|policy=="C8L3"){policylevel=3}else{}
    if(policy=="C4L4"|policy=="C8L4"){policylevel=4}else{}
    
    target_col=which(names(dur)==policy)
    Q1=quantile(dur[,target_col],0.25)
    Q2=quantile(dur[,target_col],0.50)
    Q3=quantile(dur[,target_col],0.75)
    
    if(strategy=="S1"){
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      temp_1 <- copy_L %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
      temp_2 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q1, 1, 0))
      copy_L <- bind_rows(temp_1, temp_2)
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      for (i in 1:length(copy_L$country)) {
        if(i+1<=length(copy_L$country)&&copy_L$country[i+1]==copy_L$country[i]&&copy_L$artificial_censor[i]==1){copy_L$artificial_censor[i+1]=1}else{}}
      copy_L<-rbind(subset(copy_L,artificial_censor==0),subset(copy_L,artificial_censor==1)%>%distinct(country,.keep_all=T))
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
    }else{}
    if(strategy=="S2"){
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      temp_1 <- copy_L %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
      temp_2 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q1&row_number()==max(row_number()), 1, 0))
      temp_3 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q2, 1, 0))
      copy_L <- bind_rows(temp_1, temp_2,temp_3)
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      for (i in 1:length(copy_L$country)) {
        if(i+1<=length(copy_L$country)&&copy_L$country[i+1]==copy_L$country[i]&&copy_L$artificial_censor[i]==1){copy_L$artificial_censor[i+1]=1}else{}}
      copy_L<-rbind(subset(copy_L,artificial_censor==0),subset(copy_L,artificial_censor==1)%>%distinct(country,.keep_all=T))
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
    }else{}
    if(strategy=="S3"){
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      temp_1 <- copy_L %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
      temp_2 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q2&row_number()==max(row_number()), 1, 0))
      temp_3 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q3, 1, 0))
      copy_L <- bind_rows(temp_1, temp_2,temp_3)
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      for (i in 1:length(copy_L$country)) {
        if(i+1<=length(copy_L$country)&&copy_L$country[i+1]==copy_L$country[i]&&copy_L$artificial_censor[i]==1){copy_L$artificial_censor[i+1]=1}else{}}
      copy_L<-rbind(subset(copy_L,artificial_censor==0),subset(copy_L,artificial_censor==1)%>%distinct(country,.keep_all=T))
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
    }else{} 
    if(strategy=="S4"){
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      temp_1 <- copy_L %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
      temp_2 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q3&row_number()==max(row_number()), 1, 0))
      copy_L <- bind_rows(temp_1, temp_2)
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      for (i in 1:length(copy_L$country)) {
        if(i+1<=length(copy_L$country)&&copy_L$country[i+1]==copy_L$country[i]&&copy_L$artificial_censor[i]==1){copy_L$artificial_censor[i+1]=1}else{}}
      copy_L<-rbind(subset(copy_L,artificial_censor==0),subset(copy_L,artificial_censor==1)%>%distinct(country,.keep_all=T))
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
    }else{}
    return(copy_L)
  }
  
  if(category=="C1"){
    best_C1L1=subset(best_C1L,policy=="C1L1")
    best_C1L2=subset(best_C1L,policy=="C1L2")
    best_C1L3=subset(best_C1L,policy=="C1L3")
    
    copy_1=copy.func("C1L1",best_C1L1[which(best_C1L1$endpoint==endpoint),3])%>% mutate(cgroup=1) 
    copy_2=copy.func("C1L2",best_C1L2[which(best_C1L2$endpoint==endpoint),3])%>% mutate(cgroup=2) 
    copy_3=copy.func("C1L3",best_C1L3[which(best_C1L3$endpoint==endpoint),3])%>% mutate(cgroup=3) 
    copy=rbind(copy_1,copy_2,copy_3)
  }else{}
  
  if(category=="C2"){
    best_C2L1=subset(best_C2L,policy=="C2L1")
    best_C2L2=subset(best_C2L,policy=="C2L2")
    best_C2L3=subset(best_C2L,policy=="C2L3")
    
    copy_1=copy.func("C2L1",best_C2L1[which(best_C2L1$endpoint==endpoint),3])%>% mutate(cgroup=1) 
    copy_2=copy.func("C2L2",best_C2L2[which(best_C2L2$endpoint==endpoint),3])%>% mutate(cgroup=2) 
    copy_3=copy.func("C2L3",best_C2L3[which(best_C2L3$endpoint==endpoint),3])%>% mutate(cgroup=3) 
    copy=rbind(copy_1,copy_2,copy_3)
  }else{}
  if(category=="C3"){
    best_C3L1=subset(best_C3L,policy=="C3L1")
    best_C3L2=subset(best_C3L,policy=="C3L2")
    
    copy_1=copy.func("C3L1",best_C3L1[which(best_C3L1$endpoint==endpoint),3])%>% mutate(cgroup=1) 
    copy_2=copy.func("C3L2",best_C3L2[which(best_C3L2$endpoint==endpoint),3])%>% mutate(cgroup=2) 
    
    copy=rbind(copy_1,copy_2)
  }else{}
  if(category=="C4"){
    best_C4L1=subset(best_C4L,policy=="C4L1")
    best_C4L2=subset(best_C4L,policy=="C4L2")
    best_C4L3=subset(best_C4L,policy=="C4L3")
    best_C4L4=subset(best_C4L,policy=="C4L4")
    
    copy_1=copy.func("C4L1",best_C4L1[which(best_C4L1$endpoint==endpoint),3])%>% mutate(cgroup=1) 
    copy_2=copy.func("C4L2",best_C4L2[which(best_C4L2$endpoint==endpoint),3])%>% mutate(cgroup=2) 
    copy_3=copy.func("C4L3",best_C4L3[which(best_C4L3$endpoint==endpoint),3])%>% mutate(cgroup=3) 
    copy_4=copy.func("C4L4",best_C4L4[which(best_C4L4$endpoint==endpoint),3])%>% mutate(cgroup=4)
    copy=rbind(copy_1,copy_2,copy_3,copy_4)
  }else{}
  if(category=="C5"){
    best_C5L1=subset(best_C5L,policy=="C5L1")
    best_C5L2=subset(best_C5L,policy=="C5L2")
    
    copy_1=copy.func("C5L1",best_C5L1[which(best_C5L1$endpoint==endpoint),3])%>% mutate(cgroup=1) 
    copy_2=copy.func("C5L2",best_C5L2[which(best_C5L2$endpoint==endpoint),3])%>% mutate(cgroup=2) 
    copy=rbind(copy_1,copy_2)
  }else{}
  if(category=="C6"){
    best_C6L1=subset(best_C6L,policy=="C6L1")
    best_C6L2=subset(best_C6L,policy=="C6L2")
    best_C6L3=subset(best_C6L,policy=="C6L3")
    
    copy_1=copy.func("C6L1",best_C6L1[which(best_C6L1$endpoint==endpoint),3])%>% mutate(cgroup=1) 
    copy_2=copy.func("C6L2",best_C6L2[which(best_C6L2$endpoint==endpoint),3])%>% mutate(cgroup=2) 
    copy_3=copy.func("C6L3",best_C6L3[which(best_C6L3$endpoint==endpoint),3])%>% mutate(cgroup=3) 
    copy=rbind(copy_1,copy_2,copy_3)
  }else{}
  if(category=="C7"){
    best_C7L1=subset(best_C7L,policy=="C7L1")
    best_C7L2=subset(best_C7L,policy=="C7L2")
    
    copy_1=copy.func("C7L1",best_C7L1[which(best_C7L1$endpoint==endpoint),3])%>% mutate(cgroup=1) 
    copy_2=copy.func("C7L2",best_C7L2[which(best_C7L2$endpoint==endpoint),3])%>% mutate(cgroup=2) 
    
    copy=rbind(copy_1,copy_2)
  }else{}
  if(category=="C8"){
    best_C8L1=subset(best_C8L,policy=="C8L1")
    best_C8L2=subset(best_C8L,policy=="C8L2")
    best_C8L3=subset(best_C8L,policy=="C8L3")
    best_C8L4=subset(best_C8L,policy=="C8L4")
    
    copy_1=copy.func("C8L1",best_C8L1[which(best_C8L1$endpoint==endpoint),3])%>% mutate(cgroup=1) 
    copy_2=copy.func("C8L2",best_C8L2[which(best_C8L2$endpoint==endpoint),3])%>% mutate(cgroup=2) 
    copy_3=copy.func("C8L3",best_C8L3[which(best_C8L3$endpoint==endpoint),3])%>% mutate(cgroup=3) 
    copy_4=copy.func("C8L4",best_C8L4[which(best_C8L4$endpoint==endpoint),3])%>% mutate(cgroup=4)
    copy=rbind(copy_1,copy_2,copy_3,copy_4)
  }else{}
  
  boot.fun.2 <- function(dat, index){
    #cloning
    n_group=max(copy$cgroup)
    if(n_group==2){
      copy_1<-copy[copy$cgroup==1,]
      copy_2<-copy[copy$cgroup==2,]
      copy_ref <- copy_1[index,] # allows boot to select sample
      select<-copy_ref$id
      copy_2boot<-copy_2[copy_2$id %in% select,] 
      copy<-rbind(copy_ref,copy_2boot)
    }else{}
    if(n_group==3){
      copy_1<-copy[copy$cgroup==1,]
      copy_2<-copy[copy$cgroup==2,]
      copy_3<-copy[copy$cgroup==3,]
      copy_ref <- copy_1[index,] # allows boot to select sample
      select<-copy_ref$id
      copy_2boot<-copy_2[copy_2$id %in% select,] 
      copy_3boot<-copy_3[copy_3$id %in% select,]
      copy<-rbind(copy_ref,copy_2boot,copy_3boot)}else{}
    
    if(n_group==4){
      copy_1<-copy[copy$cgroup==1,]
      copy_2<-copy[copy$cgroup==2,]
      copy_3<-copy[copy$cgroup==3,]
      copy_4<-copy[copy$cgroup==4,]
      copy_ref <- copy_1[index,] # allows boot to select sample
      select<-copy_ref$id
      copy_2boot<-copy_2[copy_2$id %in% select,] 
      copy_3boot<-copy_3[copy_3$id %in% select,]
      copy_4boot<-copy_4[copy_4$id %in% select,]
      copy<-rbind(copy_ref,copy_2boot,copy_3boot,copy_4boot)}else{}
    
    #weighting
    copy<-data.frame("id"=copy$id,"tstart"=copy$tstart,"fuptime"=copy$fuptime,"artificial_censor"=copy$artificial_censor,"policy"=copy$policy,"event"=copy$event,"cgroup"=copy$cgroup,"total_cases_per_million"=copy$total_cases_per_million,"population_density"=copy$population_density,"bed"=copy$bed,"GDP"=copy$GDPpercapita,"stringency_index"=copy$stringency_index,"glob"=copy$globalization,"total_vaccinations"=copy$total_vaccinations,"aged_over65"=copy$aged_over65)
    copy=subset(copy,is.na(copy$stringency_index)==0&is.na(copy$GDP)==0&is.na(copy$bed)==0)
    
    weight_1 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==1,])
    weight_2 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==2,])
    
    if(max(copy$cgroup)>=3){
      weight_3 <- ipwtm(
        exposure = artificial_censor,
        family = "survival",
        numerator = ~ population_density+GDP+bed+aged_over65+glob,
        denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
        id = id,
        tstart = tstart,
        timevar = fuptime,
        type = "cens",
        trunc = 0.01,
        data = copy[copy$cgroup==3,])}else{}
    if(max(copy$cgroup)==4){
      weight_4 <- ipwtm(
        exposure = artificial_censor,
        family = "survival",
        numerator = ~ population_density+GDP+bed+aged_over65+glob,
        denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
        id = id,
        tstart = tstart,
        timevar = fuptime,
        type = "cens",
        trunc = 0.01,
        data = copy[copy$cgroup==4,])}else{}
    
    
    #fit
    fita<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==1,],weights =weight_1$ipw.weights)
    fitb<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==2,],weights =weight_2$ipw.weights)
    if(max(copy$cgroup)>=3){fitc<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                                          data=copy[copy$cgroup==3,],weights =weight_3$ipw.weights)}else{}
    if(max(copy$cgroup)==4){ fitd<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                                           data=copy[copy$cgroup==4,],weights =weight_4$ipw.weights)}else{}
    
    if(max(copy$cgroup)==4){diff.res=c(min(fita$surv)-min(fitb$surv),min(fita$surv)-min(fitc$surv),min(fita$surv)-min(fitd$surv))
    names(diff.res)=c("L1-L2","L1-L3","L1-L4")}else{}
    if(max(copy$cgroup)==3){diff.res=c(min(fita$surv)-min(fitb$surv),min(fita$surv)-min(fitc$surv))
    names(diff.res)=c("L1-L2","L1-L3")}else{}
    if(max(copy$cgroup)==2){diff.res=min(fita$surv)-min(fitb$surv)
    names(diff.res)=c("L1-L2")}else{}
    
    return(diff.res)
  }
  
  if(results_type=="diff"){bootres=boot(copy, boot.fun.2, 200)}else{} #little longer
  
  #1.3 weighting
  copy<-data.frame("id"=copy$id,"tstart"=copy$tstart,"fuptime"=copy$fuptime,"artificial_censor"=copy$artificial_censor,"policy"=copy$policy,"event"=copy$event,"cgroup"=copy$cgroup,"total_cases_per_million"=copy$total_cases_per_million,"population_density"=copy$population_density,"bed"=copy$bed,"GDP"=copy$GDPpercapita,"stringency_index"=copy$stringency_index,"glob"=copy$globalization,"total_vaccinations"=copy$total_vaccinations,"aged_over65"=copy$aged_over65)
  copy=subset(copy,is.na(copy$stringency_index)==0&is.na(copy$GDP)==0&is.na(copy$bed)==0)
  
  weight_1 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==1,])
  weight_2 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==2,])
  
  if(max(copy$cgroup)>=3){
    weight_3 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==3,])}else{}
  if(max(copy$cgroup)==4){
    weight_4 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==4,])}else{}
  
  weights <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy)
  
  #2 output 
  #fit
  
  fita<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==1,],weights =weight_1$ipw.weights)
  fitb<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==2,],weights =weight_2$ipw.weights)
  
  if(max(copy$cgroup)>=3){fitc<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                                        data=copy[copy$cgroup==3,],weights =weight_3$ipw.weights)}else{}
  if(max(copy$cgroup)==4){ fitd<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                                         data=copy[copy$cgroup==4,],weights =weight_4$ipw.weights)}else{}
  
  if(results_type=="diff"){
    if(max(copy$cgroup)==4){diff.res=data.frame("diff"=as.numeric(list(0,min(fita$surv)-min(fitb$surv),min(fita$surv)-min(fitc$surv),min(fita$surv)-min(fitd$surv))),"lower"=0,"upper"=0)
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    }else{}
    
    if(max(copy$cgroup)==3){diff.res=data.frame("diff"=as.numeric(list(0,min(fita$surv)-min(fitb$surv),min(fita$surv)-min(fitc$surv))),"lower"=0,"upper"=0)
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    }else{}
    
    if(max(copy$cgroup)==2){diff.res=data.frame("diff"=as.numeric(list(0,min(fita$surv)-min(fitb$surv))),"lower"=0,"upper"=0)
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    }else{}
  }else{}
  
  if(results_type=="ptrend"){
    glm_ptrend=glm(event ~ cgroup+population_density+GDP+bed+aged_over65+glob+stringency_index, data = copy, weights = weights$ipw.weights, family = "quasibinomial")
    ptrend=coeftest(glm_ptrend, vcov=vcovHC(glm_ptrend, type="HC1")) [2,4]
    
  }else{}
  if(results_type=="HR"){
    
    max_cgroup=length(unique(copy$cgroup))
    
    glm_HR <- glm(event ~ factor(cgroup)+population_density+GDP+bed+aged_over65+glob+stringency_index , data = copy, weights = weights$ipw.weights, family = "quasibinomial")
    HR=data.frame("est"=c(1,exp(coef(glm_HR))[2:max_cgroup]),"lower"=c(1,exp(coefci(glm_HR, vcov=vcovHC(glm_HR, type="HC1")))[2:max_cgroup,1]),"upper"=c(1,exp(coefci(glm_HR, vcov=vcovHC(glm_HR, type="HC1")))[2:max_cgroup,2]))
    
    
  }else{}
  
  if(category=="C1"){
    if(results_type=="diff"){
      diff.res=data.frame("policy"=c("C1L1","C1L2","C1L3"),"comparison"=c("L1 (ref)","L1-L2","L1-L3"),diff.res)
      return(diff.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=data.frame("policy"="C1",ptrend)
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=data.frame("policy"=c("C1L1","C1L2","C1L3"),"comparison"=c("L1 (ref)","L1-L2","L1-L3"),HR)
      return(HR.res)
    }else{}
    
  }else{}
  
  if(category=="C2"){
    if(results_type=="diff"){
      diff.res=data.frame("policy"=c("C2L1","C2L2","C2L3"),"comparison"=c("L1 (ref)","L1-L2","L1-L3"),diff.res)
      return(diff.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=data.frame("policy"="C2",ptrend)
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=data.frame("policy"=c("C2L1","C2L2","C2L3"),"comparison"=c("L1 (ref)","L1-L2","L1-L3"),HR)
      return(HR.res)
    }else{}
    
  }else{}
  
  if(category=="C3"){
    if(results_type=="diff"){
      diff.res=data.frame("policy"=c("C3L1","C3L2"),"comparison"=c("L1 (ref)","L1-L2"),diff.res)
      return(diff.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=data.frame("policy"="C3",ptrend)
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=data.frame("policy"=c("C3L1","C3L2"),"comparison"=c("L1 (ref)","L1-L2"),HR)
      return(HR.res)
    }else{}
    
  }else{}
  
  
  if(category=="C4"){
    if(results_type=="diff"){
      diff.res=data.frame("policy"=c("C4L1","C4L2","C4L3","C4L4"),"comparison"=c("L1 (ref)","L1-L2","L1-L3","L1-L4"),diff.res)
      return(diff.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=data.frame("policy"="C4",ptrend)
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=data.frame("policy"=c("C4L1","C4L2","C4L3","C4L4"),"comparison"=c("L1 (ref)","L1-L2","L1-L3","L1-L4"),HR)
      return(HR.res)
    }else{}
    
  }else{}
  
  if(category=="C5"){
    if(results_type=="diff"){
      diff.res=data.frame("policy"=c("C5L1","C5L2"),"comparison"=c("L1 (ref)","L1-L2"),diff.res)
      return(diff.res)
    }else{}
    if(results_type=="ptrend"){
      ptrend.res=data.frame("policy"="C5",ptrend)
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=data.frame("policy"=c("C5L1","C5L2"),"comparison"=c("L1 (ref)","L1-L2"),HR)
      return(HR.res)
    }else{}
  }else{}
  
  
  if(category=="C6"){
    if(results_type=="diff"){
      diff.res=data.frame("policy"=c("C6L1","C6L2","C6L3"),"comparison"=c("L1 (ref)","L1-L2","L1-L3"),diff.res)
      return(diff.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=data.frame("policy"="C6",ptrend)
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=data.frame("policy"=c("C6L1","C6L2","C6L3"),"comparison"=c("L1 (ref)","L1-L2","L1-L3"),HR)
      return(HR.res)
    }else{}
  }else{}
  
  if(category=="C7"){
    if(results_type=="diff"){
      diff.res=data.frame("policy"=c("C7L1","C7L2"),"comparison"=c("L1 (ref)","L1-L2"),diff.res)
      return(diff.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=data.frame("policy"="C7",ptrend)
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=data.frame("policy"=c("C7L1","C7L2"),"comparison"=c("L1 (ref)","L1-L2"),HR)
      return(HR.res)
    }else{}
    
  }else{}
  
  if(category=="C8"){
    if(results_type=="diff"){
      diff.res=data.frame("policy"=c("C8L1","C8L2","C8L3","C8L4"),"comparison"=c("L1 (ref)","L1-L2","L1-L3","L1-L4"),diff.res)
      return(diff.res)
    }else{}
    
    if(results_type=="ptrend"){
      ptrend.res=data.frame("policy"="C8",ptrend)
      return(ptrend.res)
    }else{}
    
    if(results_type=="HR"){
      HR.res=data.frame("policy"=c("C8L1","C8L2","C8L3","C8L4"),"comparison"=c("L1 (ref)","L1-L2","L1-L3","L1-L4"),HR)
      return(HR.res)
    }else{}
    
  }else{}
}#step II: diff, HR and ptrend

step3.func=function(ana.data){
  
  #1.CCW_func
  #1.1 clone
  
  if(length(ana.data$country)==length(ana.data_3$country)){Endpoint="10^3"}else{}
  if(length(ana.data$country)==length(ana.data_4$country)){Endpoint="10^4"}else{}
  if(length(ana.data$country)==length(ana.data_5$country)){Endpoint="10^5"}else{}
  
  
  #1.2 censor
  copy.func=function(policy,strategy){
    
    if(policy=="C1L1"|policy=="C1L2"|policy=="C1L3"){
      ana.data$policy=ana.data$PolicyValue_C1
      ana.data$stringency_index=ana.data$stringency_index_C1
    }else{}
    
    if(policy=="C2L1"|policy=="C2L2"|policy=="C2L3"){
      ana.data$policy=ana.data$PolicyValue_C2
      ana.data$stringency_index=ana.data$stringency_index_C2
    }else{}
    
    if(policy=="C3L1"|policy=="C3L2"){
      ana.data$policy=ana.data$PolicyValue_C3
      ana.data$stringency_index=ana.data$stringency_index_C3
    }else{}
    
    if(policy=="C4L1"|policy=="C4L2"|policy=="C4L3"|policy=="C4L4"){
      ana.data$policy=ana.data$PolicyValue_C4
      ana.data$stringency_index=ana.data$stringency_index_C4
    }else{}
    
    if(policy=="C5L1"|policy=="C5L2"){
      ana.data$policy=ana.data$PolicyValue_C5
      ana.data$stringency_index=ana.data$stringency_index_C5
    }else{}
    
    if(policy=="C6L1"|policy=="C6L2"|policy=="C6L3"){
      ana.data$policy=ana.data$PolicyValue_C6
      ana.data$stringency_index=ana.data$stringency_index_C6
    }else{}
    
    if(policy=="C7L1"|policy=="C7L2"|policy=="C7L3"){
      ana.data$policy=ana.data$PolicyValue_C7
      ana.data$stringency_index=ana.data$stringency_index_C7
    }else{}
    
    if(policy=="C8L1"|policy=="C8L2"|policy=="C8L3"|policy=="C8L4"){
      ana.data$policy=ana.data$PolicyValue_C8
      ana.data$stringency_index=ana.data$stringency_index_C8
    }else{}
    
    if(policy=="C1L1"|policy=="C2L1"|policy=="C3L1"|policy=="C4L1"|policy=="C5L1"|policy=="C6L1"|policy=="C7L1"|policy=="C8L1"){policylevel=1}else{}
    if(policy=="C1L2"|policy=="C2L2"|policy=="C3L2"|policy=="C4L2"|policy=="C5L2"|policy=="C6L2"|policy=="C7L2"|policy=="C8L2"){policylevel=2}else{}
    if(policy=="C1L3"|policy=="C2L3"|policy=="C4L3"|policy=="C6L3"|policy=="C8L3"){policylevel=3}else{}
    if(policy=="C4L4"|policy=="C8L4"){policylevel=4}else{}
    
    target_col=which(names(dur)==policy)
    Q1=quantile(dur[,target_col],0.25)
    Q2=quantile(dur[,target_col],0.50)
    Q3=quantile(dur[,target_col],0.75)
    
    copy_L <- ana.data
    if(strategy=="S1"){
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      temp_1 <- copy_L %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
      temp_2 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q1, 1, 0))
      copy_L <- bind_rows(temp_1, temp_2)
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      for (i in 1:length(copy_L$country)) {
        if(i+1<=length(copy_L$country)&&copy_L$country[i+1]==copy_L$country[i]&&copy_L$artificial_censor[i]==1){copy_L$artificial_censor[i+1]=1}else{}}
      copy_L<-rbind(subset(copy_L,artificial_censor==0),subset(copy_L,artificial_censor==1)%>%distinct(country,.keep_all=T))
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
    }else{}
    if(strategy=="S2"){
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      temp_1 <- copy_L %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
      temp_2 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q1&row_number()==max(row_number()), 1, 0))
      temp_3 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q2, 1, 0))
      copy_L <- bind_rows(temp_1, temp_2,temp_3)
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      for (i in 1:length(copy_L$country)) {
        if(i+1<=length(copy_L$country)&&copy_L$country[i+1]==copy_L$country[i]&&copy_L$artificial_censor[i]==1){copy_L$artificial_censor[i+1]=1}else{}}
      copy_L<-rbind(subset(copy_L,artificial_censor==0),subset(copy_L,artificial_censor==1)%>%distinct(country,.keep_all=T))
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
    }else{}
    if(strategy=="S3"){
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      temp_1 <- copy_L %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
      temp_2 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q2&row_number()==max(row_number()), 1, 0))
      temp_3 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(row_number()>=Q3, 1, 0))
      copy_L <- bind_rows(temp_1, temp_2,temp_3)
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      for (i in 1:length(copy_L$country)) {
        if(i+1<=length(copy_L$country)&&copy_L$country[i+1]==copy_L$country[i]&&copy_L$artificial_censor[i]==1){copy_L$artificial_censor[i+1]=1}else{}}
      copy_L<-rbind(subset(copy_L,artificial_censor==0),subset(copy_L,artificial_censor==1)%>%distinct(country,.keep_all=T))
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
    }else{} 
    if(strategy=="S4"){
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      temp_1 <- copy_L %>% group_by(country) %>% filter(policy!=policylevel)%>%mutate(artificial_censor =0)
      temp_2 <- copy_L %>% group_by(country) %>% filter(policy==policylevel)%>%mutate(artificial_censor = ifelse(max(row_number())<Q3&row_number()==max(row_number()), 1, 0))
      copy_L <- bind_rows(temp_1, temp_2)
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
      for (i in 1:length(copy_L$country)) {
        if(i+1<=length(copy_L$country)&&copy_L$country[i+1]==copy_L$country[i]&&copy_L$artificial_censor[i]==1){copy_L$artificial_censor[i+1]=1}else{}}
      copy_L<-rbind(subset(copy_L,artificial_censor==0),subset(copy_L,artificial_censor==1)%>%distinct(country,.keep_all=T))
      copy_L<-copy_L[order(copy_L$id,copy_L$tstartdate),]
    }else{}
    return(copy_L)
  }
  
  copy_1=copy.func(best_C1$policy[which(best_C1$endpoint==Endpoint)],best_C1$strategy[which(best_C1$endpoint==Endpoint)])%>% mutate(cgroup=1) 
  copy_2=copy.func(best_C2$policy[which(best_C2$endpoint==Endpoint)],best_C2$strategy[which(best_C2$endpoint==Endpoint)])%>% mutate(cgroup=2) 
  copy_3=copy.func(best_C3$policy[which(best_C3$endpoint==Endpoint)],best_C3$strategy[which(best_C3$endpoint==Endpoint)])%>% mutate(cgroup=3) 
  copy_4=copy.func(best_C4$policy[which(best_C4$endpoint==Endpoint)],best_C4$strategy[which(best_C4$endpoint==Endpoint)])%>% mutate(cgroup=4)
  copy_5=copy.func(best_C5$policy[which(best_C5$endpoint==Endpoint)],best_C5$strategy[which(best_C5$endpoint==Endpoint)])%>% mutate(cgroup=5) 
  copy_6=copy.func(best_C6$policy[which(best_C6$endpoint==Endpoint)],best_C6$strategy[which(best_C6$endpoint==Endpoint)])%>% mutate(cgroup=6) 
  copy_7=copy.func(best_C7$policy[which(best_C7$endpoint==Endpoint)],best_C7$strategy[which(best_C7$endpoint==Endpoint)])%>% mutate(cgroup=7) 
  copy_8=copy.func(best_C8$policy[which(best_C8$endpoint==Endpoint)],best_C8$strategy[which(best_C8$endpoint==Endpoint)])%>% mutate(cgroup=8)
  
  copy<-rbind(copy_1,copy_2,copy_3,copy_4,copy_5,copy_6,copy_7,copy_8)
  
  best_C=rbind(best_C1[which(best_C1$endpoint==Endpoint),],
               best_C2[which(best_C2$endpoint==Endpoint),],
               best_C3[which(best_C3$endpoint==Endpoint),],
               best_C4[which(best_C4$endpoint==Endpoint),],
               best_C5[which(best_C5$endpoint==Endpoint),],
               best_C6[which(best_C6$endpoint==Endpoint),],
               best_C7[which(best_C7$endpoint==Endpoint),],
               best_C8[which(best_C8$endpoint==Endpoint),])
  
  ref=which(best_C$incidence==max(best_C$incidence))
  
  boot.fun.3 <- function(dat, index){
    #cloning
    copy_1<-copy[copy$cgroup==1,]
    copy_2<-copy[copy$cgroup==2,]
    copy_3<-copy[copy$cgroup==3,]
    copy_4<-copy[copy$cgroup==4,]
    copy_5<-copy[copy$cgroup==5,]
    copy_6<-copy[copy$cgroup==6,]
    copy_7<-copy[copy$cgroup==7,]
    copy_8<-copy[copy$cgroup==8,]
    
    if(best_C$policy[ref]=="C1L1"|best_C$policy[ref]=="C1L2"|best_C$policy[ref]=="C1L3"){copy_ref <- copy_1[index,]}else{}
    if(best_C$policy[ref]=="C2L1"|best_C$policy[ref]=="C2L2"|best_C$policy[ref]=="C2L3"){copy_ref <- copy_2[index,]}else{}
    if(best_C$policy[ref]=="C3L1"|best_C$policy[ref]=="C3L2"){copy_ref <- copy_3[index,]}else{}
    if(best_C$policy[ref]=="C4L1"|best_C$policy[ref]=="C4L2"|best_C$policy[ref]=="C4L3"|best_C$policy[ref]=="C4L4"){copy_ref <- copy_4[index,]}else{}
    if(best_C$policy[ref]=="C5L1"|best_C$policy[ref]=="C5L2"){copy_ref <- copy_5[index,]}else{}
    if(best_C$policy[ref]=="C6L1"|best_C$policy[ref]=="C6L2"|best_C$policy[ref]=="C6L3"){copy_ref <- copy_6[index,]}else{}
    if(best_C$policy[ref]=="C7L1"|best_C$policy[ref]=="C7L2"){copy_ref <- copy_7[index,]}else{}    
    if(best_C$policy[ref]=="C8L1"|best_C$policy[ref]=="C8L2"|best_C$policy[ref]=="C8L3"|best_C$policy[ref]=="C8L4"){copy_ref <- copy_8[index,]}else{}
    
    select<-copy_ref$id
    
    if(best_C$policy[ref]!="C1L1"&best_C$policy[ref]!="C1L2"&best_C$policy[ref]!="C1L3"){    copy_1boot<-copy_1[copy_1$id %in% select,]}else{copy_1boot<-copy_ref}
    if(best_C$policy[ref]!="C2L1"&best_C$policy[ref]!="C2L2"&best_C$policy[ref]!="C2L3"){    copy_2boot<-copy_2[copy_2$id %in% select,]}else{copy_2boot<-copy_ref}
    if(best_C$policy[ref]!="C3L1"&best_C$policy[ref]!="C3L2"){    copy_3boot<-copy_3[copy_3$id %in% select,]}else{copy_3boot<-copy_ref}
    if(best_C$policy[ref]!="C4L1"&best_C$policy[ref]!="C4L2"&best_C$policy[ref]!="C4L3"&best_C$policy[ref]!="C4L4"){    copy_4boot<-copy_4[copy_4$id %in% select,]}else{copy_4boot<-copy_ref}
    if(best_C$policy[ref]!="C5L1"&best_C$policy[ref]!="C5L2"){    copy_5boot<-copy_5[copy_5$id %in% select,]}else{copy_5boot<-copy_ref}
    if(best_C$policy[ref]!="C6L1"&best_C$policy[ref]!="C6L2"&best_C$policy[ref]!="C6L3"){    copy_6boot<-copy_6[copy_6$id %in% select,]}else{copy_6boot<-copy_ref}
    if(best_C$policy[ref]!="C7L1"&best_C$policy[ref]!="C7L2"){    copy_7boot<-copy_7[copy_7$id %in% select,]}else{copy_7boot<-copy_ref}    
    if(best_C$policy[ref]!="C8L1"&best_C$policy[ref]!="C8L2"&best_C$policy[ref]!="C8L3"&best_C$policy[ref]!="C8L4"){    copy_8boot<-copy_8[copy_8$id %in% select,]}else{copy_8boot<-copy_ref}
    
    copy<-rbind(copy_1boot,copy_2boot,copy_3boot,copy_4boot,copy_5boot,copy_6boot,copy_7boot,copy_8boot)
    
    #weighting
    copy<-data.frame("id"=copy$id,"tstart"=copy$tstart,"fuptime"=copy$fuptime,"artificial_censor"=copy$artificial_censor,"policy"=copy$policy,"event"=copy$event,"cgroup"=copy$cgroup,"total_cases_per_million"=copy$total_cases_per_million,"population_density"=copy$population_density,"bed"=copy$bed,"GDP"=copy$GDPpercapita,"stringency_index"=copy$stringency_index,"glob"=copy$globalization,"total_vaccinations"=copy$total_vaccinations,"aged_over65"=copy$aged_over65)
    copy=subset(copy,is.na(copy$stringency_index)==0&is.na(copy$GDP)==0&is.na(copy$bed)==0)
    
    weight_1 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==1,])
    weight_2 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==2,])
    weight_3 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==3,])
    weight_4 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==4,])
    weight_5 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==5,])
    weight_6 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==6,])
    weight_7 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==7,])
    weight_8 <- ipwtm(
      exposure = artificial_censor,
      family = "survival",
      numerator = ~ population_density+GDP+bed+aged_over65+glob,
      denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "cens",
      trunc = 0.01,
      data = copy[copy$cgroup==8,])
    
    #fit
    fita<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==1,],weights =weight_1$ipw.weights)
    fitb<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==2,],weights =weight_2$ipw.weights)
    fitc<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==3,],weights =weight_3$ipw.weights)
    fitd<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==4,],weights =weight_4$ipw.weights)
    fite<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==5,],weights =weight_5$ipw.weights)
    fitf<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==6,],weights =weight_6$ipw.weights)
    fitg<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==7,],weights =weight_7$ipw.weights)
    fith<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                  data=copy[copy$cgroup==8,],weights =weight_8$ipw.weights)
    
    if(best_C[ref,2]=="C1L1"|best_C[ref,2]=="C1L2"|best_C[ref,2]=="C1L3"){
      diff.res=c(min(fita$surv)-min(fitb$surv),min(fita$surv)-min(fitc$surv),min(fita$surv)-min(fitd$surv),min(fita$surv)-min(fite$surv),min(fita$surv)-min(fitf$surv),min(fita$surv)-min(fitg$surv),min(fita$surv)-min(fith$surv))
      names(diff.res)=c("C2-C1","C3-C1","C4-C1","C5-C1","C6-C1","C7-C1","C8-C1")
    }else{}
    if(best_C[ref,2]=="C2L1"|best_C[ref,2]=="C2L2"|best_C[ref,2]=="C2L3"){
      diff.res=c(min(fitb$surv)-min(fita$surv),min(fitb$surv)-min(fitc$surv),min(fitb$surv)-min(fitd$surv),min(fitb$surv)-min(fite$surv),min(fitb$surv)-min(fitf$surv),min(fitb$surv)-min(fitg$surv),min(fitb$surv)-min(fith$surv))
      names(diff.res)=c("C1-C2","C3-C2","C4-C2","C5-C2","C6-C2","C7-C2","C8-C2")
    }else{}
    if(best_C[ref,2]=="C3L1"|best_C[ref,2]=="C3L2"){
      diff.res=c(min(fitc$surv)-min(fita$surv),min(fitc$surv)-min(fitb$surv),min(fitc$surv)-min(fitd$surv),min(fitc$surv)-min(fite$surv),min(fitc$surv)-min(fitf$surv),min(fitc$surv)-min(fitg$surv),min(fitc$surv)-min(fith$surv))
      names(diff.res)=c("C1-C3","C2-C3","C4-C3","C5-C3","C6-C3","C7-C3","C8-C3")
    }else{}
    if(best_C[ref,2]=="C4L1"|best_C[ref,2]=="C4L2"|best_C[ref,2]=="C4L3"|best_C[ref,2]=="C4L4"){
      diff.res=c(min(fitd$surv)-min(fita$surv),min(fitd$surv)-min(fitb$surv),min(fitd$surv)-min(fitc$surv),min(fitd$surv)-min(fite$surv),min(fitd$surv)-min(fitf$surv),min(fitd$surv)-min(fitg$surv),min(fitd$surv)-min(fith$surv))
      names(diff.res)=c("C1-C4","C2-C4","C3-C4","C5-C4","C6-C4","C7-C4","C8-C4")
    }else{}
    if(best_C[ref,2]=="C5L1"|best_C[ref,2]=="C5L2"){
      diff.res=c(min(fite$surv)-min(fita$surv),min(fite$surv)-min(fitb$surv),min(fite$surv)-min(fitc$surv),min(fite$surv)-min(fitd$surv),min(fite$surv)-min(fitf$surv),min(fite$surv)-min(fitg$surv),min(fite$surv)-min(fith$surv))
      names(diff.res)=c("C1-C5","C2-C5","C3-C5","C4-C5","C6-C5","C7-C5","C8-C5")
    }else{}
    if(best_C[ref,2]=="C6L1"|best_C[ref,2]=="C6L2"|best_C[ref,2]=="C6L3"){
      diff.res=c(min(fitf$surv)-min(fita$surv),min(fitf$surv)-min(fitb$surv),min(fitf$surv)-min(fitc$surv),min(fitf$surv)-min(fitd$surv),min(fitf$surv)-min(fite$surv),min(fitf$surv)-min(fitg$surv),min(fitf$surv)-min(fith$surv))
      names(diff.res)=c("C1-C6","C2-C6","C3-C6","C4-C6","C5-C6","C7-C6","C8-C6")
    }else{}
    if(best_C[ref,2]=="C7L1"|best_C[ref,2]=="C7L2"){
      diff.res=c(min(fitg$surv)-min(fita$surv),min(fitg$surv)-min(fitb$surv),min(fitg$surv)-min(fitc$surv),min(fitg$surv)-min(fitd$surv),min(fitg$surv)-min(fite$surv),min(fitg$surv)-min(fitf$surv),min(fitg$surv)-min(fith$surv))
      names(diff.res)=c("C1-C7","C2-C7","C3-C7","C4-C7","C5-C7","C6-C7","C8-C7")
    }else{}    
    if(best_C[ref,2]=="C8L1"|best_C[ref,2]=="C8L2"|best_C[ref,2]=="C8L3"|best_C[ref,2]=="C8L4"){
      diff.res=c(min(fith$surv)-min(fita$surv),min(fith$surv)-min(fitb$surv),min(fith$surv)-min(fitc$surv),min(fith$surv)-min(fitd$surv),min(fith$surv)-min(fite$surv),min(fith$surv)-min(fitf$surv),min(fith$surv)-min(fitg$surv))
      names(diff.res)=c("C1-C8","C2-C8","C3-C8","C4-C8","C5-C8","C6-C8","C7-C8")
    }else{}
    
    return(diff.res)
  }
  
  bootres=boot(copy, boot.fun.3, 200)
  
  #1.3 weight
  copy<-data.frame("id"=copy$id,"tstart"=copy$tstart,"fuptime"=copy$fuptime,"artificial_censor"=copy$artificial_censor,"policy"=copy$policy,"event"=copy$event,"cgroup"=copy$cgroup,"total_cases_per_million"=copy$total_cases_per_million,"population_density"=copy$population_density,"bed"=copy$bed,"GDP"=copy$GDPpercapita,"stringency_index"=copy$stringency_index,"glob"=copy$globalization,"total_vaccinations"=copy$total_vaccinations,"aged_over65"=copy$aged_over65)
  copy=subset(copy,is.na(copy$stringency_index)==0&is.na(copy$GDP)==0&is.na(copy$bed)==0)
  weight_1 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==1,])
  weight_2 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==2,])
  weight_3 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==3,])
  weight_4 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==4,])
  weight_5 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==5,])
  weight_6 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==6,])
  weight_7 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==7,])
  weight_8 <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy[copy$cgroup==8,])
  
  weights <- ipwtm(
    exposure = artificial_censor,
    family = "survival",
    numerator = ~ population_density+GDP+bed+aged_over65+glob,
    denominator = ~population_density+GDP+bed+aged_over65+glob+stringency_index,
    id = id,
    tstart = tstart,
    timevar = fuptime,
    type = "cens",
    trunc = 0.01,
    data = copy)
  
  
  #2 output 
  #fit
  
  fita<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==1,],weights =weight_1$ipw.weights)
  fitb<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==2,],weights =weight_2$ipw.weights)
  fitc<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==3,],weights =weight_3$ipw.weights)
  fitd<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==4,],weights =weight_4$ipw.weights)
  fite<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==5,],weights =weight_5$ipw.weights)
  fitf<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==6,],weights =weight_6$ipw.weights)
  fitg<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==7,],weights =weight_7$ipw.weights)
  fith<-survfit(Surv(tstart, fuptime, event) ~1+ cluster(id),
                data=copy[copy$cgroup==8,],weights =weight_8$ipw.weights)
  
  
  if(best_C[ref,2]=="C1L1"|best_C[ref,2]=="C1L2"|best_C[ref,2]=="C1L3"){
    
    diff.res=data.frame("diff"=as.numeric(list(0,min(fita$surv)-min(fitb$surv),min(fita$surv)-min(fitc$surv),min(fita$surv)-min(fitd$surv),min(fita$surv)-min(fite$surv),min(fita$surv)-min(fitf$surv),min(fita$surv)-min(fitg$surv),min(fita$surv)-min(fith$surv))),
                        "lower"=0,"upper"=0)
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95, type = "norm",index=2)$normal[1,2:3]
    diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    diff.res[5,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=4)$normal[1,2:3]
    diff.res[6,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=5)$normal[1,2:3]
    diff.res[7,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=6)$normal[1,2:3]
    diff.res[8,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=7)$normal[1,2:3]
    comparison_C=c("C1 (ref)","C2-C1","C3-C1","C4-C1","C5-C1","C6-C1","C7-C1","C8-C1")
  }else{}
  
  if(best_C[ref,2]=="C2L1"|best_C[ref,2]=="C2L2"|best_C[ref,2]=="C2L3"){
    diff.res=data.frame("diff"=as.numeric(list(min(fitb$surv)-min(fita$surv),0,min(fitb$surv)-min(fitc$surv),min(fitb$surv)-min(fitd$surv),min(fitb$surv)-min(fite$surv),min(fitb$surv)-min(fitf$surv),min(fitb$surv)-min(fitg$surv),min(fitb$surv)-min(fith$surv))),
                        "lower"=0,"upper"=0)
    diff.res[1,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    diff.res[5,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=4)$normal[1,2:3]
    diff.res[6,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=5)$normal[1,2:3]
    diff.res[7,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=6)$normal[1,2:3]
    diff.res[8,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=7)$normal[1,2:3]
    
    comparison_C=c("C1-C2","C2 (ref)","C3-C2","C4-C2","C5-C2","C6-C2","C7-C2","C8-C2")
    
  }else{}
  
  if(best_C[ref,2]=="C3L1"|best_C[ref,2]=="C3L2"){
    diff.res=data.frame("diff"=as.numeric(list(min(fitc$surv)-min(fita$surv),min(fitc$surv)-min(fitb$surv),0,min(fitc$surv)-min(fitd$surv),min(fitc$surv)-min(fite$surv),min(fitc$surv)-min(fitf$surv),min(fitc$surv)-min(fitg$surv),min(fitc$surv)-min(fith$surv))),
                        "lower"=0,"upper"=0)
    diff.res[1,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    diff.res[5,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=4)$normal[1,2:3]
    diff.res[6,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=5)$normal[1,2:3]
    diff.res[7,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=6)$normal[1,2:3]
    diff.res[8,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=7)$normal[1,2:3]
    
    comparison_C=c("C1-C3","C2-C3","C3 (ref)","C4-C3","C5-C3","C6-C3","C7-C3","C8-C3")
    
  }else{}
  
  if(best_C[ref,2]=="C4L1"|best_C[ref,2]=="C4L2"|best_C[ref,2]=="C4L3"|best_C[ref,2]=="C4L4"){
    diff.res=data.frame("diff"=as.numeric(list(min(fitd$surv)-min(fita$surv),min(fitd$surv)-min(fitb$surv),min(fitd$surv)-min(fitc$surv),0,min(fitd$surv)-min(fite$surv),min(fitd$surv)-min(fitf$surv),min(fitd$surv)-min(fitg$surv),min(fitd$surv)-min(fith$surv))),
                        "lower"=0,"upper"=0)
    diff.res[1,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    diff.res[5,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=4)$normal[1,2:3]
    diff.res[6,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=5)$normal[1,2:3]
    diff.res[7,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=6)$normal[1,2:3]
    diff.res[8,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=7)$normal[1,2:3]
    
    comparison_C=c("C1-C4","C2-C4","C3-C4","C4 (ref)","C5-C4","C6-C4","C7-C4","C8-C4")
  }else{}
  
  if(best_C[ref,2]=="C5L1"|best_C[ref,2]=="C5L2"){
    diff.res=data.frame("diff"=as.numeric(list(min(fite$surv)-min(fita$surv),min(fite$surv)-min(fitb$surv),min(fite$surv)-min(fitc$surv),min(fite$surv)-min(fitd$surv),0,min(fite$surv)-min(fitf$surv),min(fite$surv)-min(fitg$surv),min(fite$surv)-min(fith$surv))),
                        "lower"=0,"upper"=0)
    diff.res[1,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=4)$normal[1,2:3]
    diff.res[6,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=5)$normal[1,2:3]
    diff.res[7,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=6)$normal[1,2:3]
    diff.res[8,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=7)$normal[1,2:3]
    
    comparison_C=c("C1-C5","C2-C5","C3-C5","C4-C5","C5 (ref)","C6-C5","C7-C5","C8-C5")
    
  }else{}
  
  if(best_C[ref,2]=="C6L1"|best_C[ref,2]=="C6L2"|best_C[ref,2]=="C6L3"){
    
    diff.res=data.frame("diff"=as.numeric(list(min(fitf$surv)-min(fita$surv),min(fitf$surv)-min(fitb$surv),min(fitf$surv)-min(fitc$surv),min(fitf$surv)-min(fitd$surv),min(fitf$surv)-min(fite$surv),0,min(fitf$surv)-min(fitg$surv),min(fitf$surv)-min(fith$surv))),
                        "lower"=0,"upper"=0)
    diff.res[1,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=4)$normal[1,2:3]
    diff.res[5,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=5)$normal[1,2:3]
    diff.res[7,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=6)$normal[1,2:3]
    diff.res[8,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=7)$normal[1,2:3]
    
    comparison_C=c("C1-C6","C2-C6","C3-C6","C4-C6","C5-C6","C6 (ref)","C7-C6","C8-C6")
  }else{}
  
  if(best_C[ref,2]=="C7L1"|best_C[ref,2]=="C7L2"){
    diff.res=data.frame("diff"=as.numeric(list(min(fitg$surv)-min(fita$surv),min(fitg$surv)-min(fitb$surv),min(fitg$surv)-min(fitc$surv),min(fitg$surv)-min(fitd$surv),min(fitg$surv)-min(fite$surv),min(fitg$surv)-min(fitf$surv),0,min(fitg$surv)-min(fith$surv))),
                        "lower"=0,"upper"=0)
    diff.res[1,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=4)$normal[1,2:3]
    diff.res[5,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=5)$normal[1,2:3]
    diff.res[6,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=6)$normal[1,2:3]
    diff.res[8,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=7)$normal[1,2:3]
    
    comparison_C=c("C1-C7","C2-C7","C3-C7","C4-C7","C5-C7","C6-C7","C7 (ref)","C8-C7")
  }else{}
  
  if(best_C[ref,2]=="C8L1"|best_C[ref,2]=="C8L2"|best_C[ref,2]=="C8L3"|best_C[ref,2]=="C8L4"){
    diff.res=data.frame("diff"=as.numeric(list(min(fith$surv)-min(fita$surv),min(fith$surv)-min(fitb$surv),min(fith$surv)-min(fitc$surv),min(fith$surv)-min(fitd$surv),min(fith$surv)-min(fite$surv),min(fith$surv)-min(fitf$surv),min(fith$surv)-min(fitg$surv),0)),
                        "lower"=0,"upper"=0)
    diff.res[1,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=1)$normal[1,2:3]
    diff.res[2,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=2)$normal[1,2:3]
    diff.res[3,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=3)$normal[1,2:3]
    diff.res[4,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=4)$normal[1,2:3]
    diff.res[5,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=5)$normal[1,2:3]
    diff.res[6,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=6)$normal[1,2:3]
    diff.res[7,2:3]=boot.ci(bootres, conf = 0.95,type = "norm",index=7)$normal[1,2:3]
    
    comparison_C=c("C1-C8","C2-C8","C3-C8","C4-C8","C5-C8","C6-C8","C7-C8","C8 (ref)")
    
  }else{}
  
  Res=data.frame("policy"=best_C$policy,
                 "strategy"=best_C$strategy,
                 "comparison"=comparison_C,
                 diff.res)
  
  return(Res)
}#step III: diff