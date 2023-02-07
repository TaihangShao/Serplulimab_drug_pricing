#library
library(survHE)
library(survival)
library(survminer)
library(tidyverse)
library(ggplot2)
library("discSurv")    # Create life tables
# library(ggsci)

###########################-----###
#For OS
###########################-----###
cam <- read.delim("pfs_che.txt")
cam<-data.frame(cam$event,cam$time)
cam<-rename(cam,"censrec"="cam.event","recyrs"="cam.time")
cam$recyrs<-cam$recyrs/12
cam$recyrs2<-cam$recyrs*365.25/21
cam$rectime<-as.integer(cam$recyrs*365.25)

table(cam$censrec)
cam %>%
  summarise(Min_surv = min(recyrs, na.rm = TRUE),
            Max_surv = max(recyrs, na.rm = TRUE),
            Mean_surv = mean(recyrs, na.rm = TRUE),
            SD_surv = sd(recyrs, na.rm = TRUE))

##-----niv data 3Y Monthly life table estimates of hazard-----###
cam$rectime2 <- as.integer(cam$recyrs2) + 1
#+1 above as integer rounds down, we want to round up
ltBC <- lifeTable(cam, timeColumn = "rectime2", eventColumn = "censrec")
ltHaz0 <- data.frame(hazKM = ltBC$Output$hazard, Time = (seq(1:length(ltBC$Output[,1])))*21/365,
                     AtRisk = ltBC$Output$atRisk, Events = ltBC$Output$events)
# The above hazard is the product-limit (KM) estimate. Also calculate the life-table (acturial) estimate
ltHaz0$hazLT = ltHaz0$Events / (ltHaz0$AtRisk - ltHaz0$Events/2)
# Generate log-time
ltHaz0$lnTime <- log(ltHaz0$Time)
# For random effects add an ID for each time period
ltHaz0$MyId <- 1:dim(ltHaz0)[1] # Generate id variable 
# For AR(1) model get outcomes lagged by one.
ltHaz0$EventsL <- lag(ltHaz0$Events)
# Set first lagged value = 0 (usually would discard, but retain so IC are comparable. Can be justified as a prior value)
ltHaz0$EventsL[1] <- 0
#Set surv data
ltHaz0$surv <- ltBC$Output$S
#timedelta
ltHaz0$timedelta<-ltHaz0$Time[2]-ltHaz0$Time[1]
ltHaz0<-ltHaz0[1:39,]
###########################-----###
###########################-----###
###########################-----###
###
pfs1 <- read.delim("pfs_che.txt")
# pfs1$recys=pfs1$time/12
fit_pfs1<-survfit(Surv(pfs1$time,pfs1$event)~1,data=pfs1)
data_fit_pfs1<-data.frame(time=fit_pfs1$time/12, surv=fit_pfs1$surv,lower=fit_pfs1$lower,upper=fit_pfs1$upper)

reference0<-data.frame(Time=data_fit_pfs1$time,Surv=data_fit_pfs1$surv)
prime2<-c(0,1)
reference0<-rbind(prime2,reference0)

###
pfs0 <- read.delim("pfs_c.txt")
# pfs0$recys=pfs0$time/12
fit_pfs0<-survfit(Surv(pfs0$time,pfs0$event)~1,data=pfs0)
data_fit_pfs0<-data.frame(time=fit_pfs0$time/12, surv=fit_pfs0$surv,lower=fit_pfs0$lower,upper=fit_pfs0$upper)

reference1<-data.frame(Time=data_fit_pfs0$time,Surv=data_fit_pfs0$surv)
prime2<-c(0,1)
reference1<-rbind(prime2,reference1)

###
dfsurv_pfs_ser<-read.csv("pfs_che_surv.csv")
dfsurv_pfs_ser<-dfsurv_pfs_ser[1:39,]
dfsurv_pfs_ser<-dfsurv_pfs_ser[,-1]
dfFigSurv1 = dfsurv_pfs_ser %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))

f_surv_mod1= ggplot() +
  geom_ribbon(data = data_fit_pfs1,aes(x=time,ymin=lower,ymax=upper),alpha = 0.3,show.legend = FALSE,fill="gray50")+
  geom_line(data = data_fit_pfs1,aes(x=time,y=lower),linetype = 2,colour="gray50")+
  geom_line(data=data_fit_pfs1,aes(x=time,y=upper),linetype = 2,colour="gray50")+
  geom_ribbon(data = data_fit_pfs0,aes(x=time,ymin=lower,ymax=upper),alpha = 0.5,show.legend = FALSE,fill="pink")+
  geom_line(data = data_fit_pfs0,aes(x=time,y=lower),linetype = 2,colour="red")+
  geom_line(data=data_fit_pfs0,aes(x=time,y=upper),linetype = 2,colour="red")+
  geom_line(data=dfFigSurv1, aes(x=Time, y=survProp, group=Model, colour=Model), size=1.2,show.legend = FALSE) +
  geom_line(data=reference0, aes(x=Time, y=Surv,colour="KM"), size=0.75,colour="Black")+
  geom_line(data=reference1, aes(x=Time, y=Surv,colour="KM"), size=0.75,colour="red")+
  scale_color_discrete(name="Model")+
  # scale_color_lancet()+
  expand_limits(y=c(0,1),x=c(0,3)) + 
  facet_wrap(~Model,nrow=4)+
  scale_x_continuous(breaks = c(seq(from=0, to=3,by = 1))) +
  ylab("Ser_PFS") +
  xlab("Time(Years)") +
  guides(color = guide_legend(ncol = 1))  +
  theme_bw()+
  theme(legend.position = "bottom") 
f_surv_mod1

###########################-----###
#test original
###########################-----###

ser <- read.delim("pfs_c.txt")
ser<-data.frame(ser$event,ser$time)
ser<-rename(ser,"censrec"="ser.event","recyrs"="ser.time")
ser$recyrs<-ser$recyrs/12
ser$recyrs2<-ser$recyrs*365.25/21
ser$rectime<-as.integer(ser$recyrs*365.25)

table(ser$censrec)
ser %>%
  summarise(Min_surv = min(recyrs, na.rm = TRUE),
            Max_surv = max(recyrs, na.rm = TRUE),
            Mean_surv = mean(recyrs, na.rm = TRUE),
            SD_surv = sd(recyrs, na.rm = TRUE))

##-----niv data 3Y Monthly life table estimates of hazard-----###
ser$rectime2 <- as.integer(ser$recyrs2) + 1
#+1 above as integer rounds down, we want to round up
ltBC <- lifeTable(ser, timeColumn = "rectime2", eventColumn = "censrec")
ltHaz1 <- data.frame(hazKM = ltBC$Output$hazard, Time = (seq(1:length(ltBC$Output[,1])))*21/365,
                     AtRisk = ltBC$Output$atRisk, Events = ltBC$Output$events)
# The above hazard is the product-limit (KM) estimate. Also calculate the life-table (acturial) estimate
ltHaz1$hazLT = ltHaz1$Events / (ltHaz1$AtRisk - ltHaz1$Events/2)
# Generate log-time
ltHaz1$lnTime <- log(ltHaz1$Time)
# For random effects add an ID for each time period
ltHaz1$MyId <- 1:dim(ltHaz1)[1] # Generate id variable 
# For AR(1) model get outcomes lagged by one.
ltHaz1$EventsL <- lag(ltHaz1$Events)
# Set first lagged value = 0 (usually would discard, but retain so IC are comparable. Can be justified as a prior value)
ltHaz1$EventsL[1] <- 0
#Set surv data
ltHaz1$surv <- ltBC$Output$S
#timedelta
ltHaz1$timedelta<-ltHaz1$Time[2]-ltHaz1$Time[1]

cal_df<-dfsurv_pfs_ser[1:35,]
cal_lthaz<-ltHaz1
cal_lthaz<-data.frame(Time=cal_lthaz$Time,Surv=cal_lthaz$surv)

md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1-1","FP1-2",
      "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")

MSE0<-as.data.frame(matrix(nrow = 1,ncol = 19))
colnames(MSE0) <- md
for (i in 1:19) {
  MSE0[i]<-mean((cal_df[,i+1]-cal_lthaz$Surv)^2)
}
MSE0<-MSE0*10000

dfplotMSE01<-data.frame(name=colnames(MSE0),MSE0=t(MSE0[1,]))
colnames(dfplotMSE01)<-c("name","MSE")
MSE0_legend<-as.character(round(MSE0[1,],digits=2))
fig_MSE0_1<-ggplot(data=dfplotMSE01,mapping=aes(x=name,y=MSE,fill=name,group=factor(1)))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label = MSE0_legend, vjust = -0.8, hjust = 0.5, color = name), show.legend = TRUE)
fig_MSE0_1



###########################-----###
#test 12:39
###########################-----###
cal_df<-dfsurv_pfs_ser[12:39,]
cal_lthaz<-ltHaz0[12:39,]
cal_lthaz<-data.frame(Time=cal_lthaz$Time,Surv=cal_lthaz$surv)

md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1-1","FP1-2",
      "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")

MSE<-as.data.frame(matrix(nrow = 1,ncol = 19))
colnames(MSE) <- md
for (i in 1:19) {
  MSE[i]<-mean((cal_df[,i+1]-cal_lthaz$Surv)^2)
}
MSE<-MSE*10000

dfplotmse1<-data.frame(name=colnames(MSE),MSE=t(MSE[1,]))
colnames(dfplotmse1)<-c("name","MSE")
MSE_legend<-as.character(round(MSE[1,],digits=2))
fig_mse_1<-ggplot(data=dfplotmse1,mapping=aes(x=name,y=MSE,fill=name,group=factor(1)))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label = MSE_legend, vjust = -0.8, hjust = 0.5, color = name), show.legend = TRUE)
  # scale_fill_jama()
fig_mse_1

ggsave(f_surv_mod1,width = 12,height = 8,dpi = 800,filename = "et_pfs_che.png",path = "C:/Users/Administrator/Desktop/斯鲁利单抗/external validation/pfs_che")
ggsave(fig_MSE0_1,width = 12,height = 8,dpi = 800,filename = "pfs_che_mse.png",path = "C:/Users/Administrator/Desktop/斯鲁利单抗/external validation/pfs_che")
ggsave(fig_mse_1,width = 12,height = 8,dpi = 800,filename = "et_mse_pfs_che.png",path = "C:/Users/Administrator/Desktop/斯鲁利单抗/external validation/pfs_che")
res_mse=rbind(MSE,MSE0)
write.csv(res_mse,file = "C:/Users/Administrator/Desktop/斯鲁利单抗/external validation/pfs_che/MSE.csv",row.names=FALSE)



