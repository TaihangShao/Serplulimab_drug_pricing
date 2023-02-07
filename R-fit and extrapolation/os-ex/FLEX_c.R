library("ggplot2")   # Mainly for plotting
library("tidyverse")   # Mainly for plotting
library("flexsurv")    # RPs (also loads survival) and BC case-study
library("gridExtra")   # Plotting
library("survHE")   # Plotting
library("discSurv")    # Create life tables
library("mgcv")        # GAM/RCS
library("survminer")   #KM curve


theme_set(theme_light())  # GGplot theme

## Load data base

che <- read.delim("os_c.txt")
che<-data.frame(che$event,che$time)
che<-rename(che,"censrec"="che.event","recyrs"="che.time")
che$recyrs<-che$recyrs/12
che$recyrs2<-che$recyrs*365.25/21
che$rectime<-as.integer(che$recyrs*365.25) 


###----------------------------------------------------------- Process data ------------------------------------------------------------------------------###
table(che$censrec)
che %>%
  summarise(Min_surv = min(recyrs, na.rm = TRUE),
            Max_surv = max(recyrs, na.rm = TRUE),
            Mean_surv = mean(recyrs, na.rm = TRUE),
            SD_surv = sd(recyrs, na.rm = TRUE))


##-----ipi data 3Y Monthly life table estimates of hazard-----##
che$rectime2 <- as.integer(che$recyrs2) + 1
#+1 above as integer rounds down, we want to round up
ltBC1 <- lifeTable(che, timeColumn = "rectime2", eventColumn = "censrec")
ltHaz1 <- data.frame(hazKM = ltBC1$Output$hazard, Time = (seq(1:length(ltBC1$Output[,1]))-0.5)*21/365,
                    AtRisk = ltBC1$Output$atRisk, Events = ltBC1$Output$events)
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
ltHaz1$surv <- ltBC1$Output$S
#timedelta
ltHaz1$timedelta<-ltHaz1$Time[2]-ltHaz1$Time[1]

###--------------------------------------------------------------------------------------------------------------------------------------------------------####
###--------------------------------------------------------------------------------------------------------------------------------------------------------####
###########################################                                 CHE                                     ###########################################
###--------------------------------------------------------------------------------------------------------------------------------------------------------####
###--------------------------------------------------------------------------------------------------------------------------------------------------------####
# ipi ##
ltHaz<-ltHaz1
sd_bc<-data.frame("recyrs"=che$recyrs,"censrec"=che$censrec)

####----New Data----####
follow_up <- 34
numMod <- 19 # Models considered
MyTH <- 20 # Time Horizon (years)
MyStep <- 17.4 # Number of obs. per year
MyN <- MyTH*MyStep # Total time points (observed & extrapolated)
dfHazEst <- array(dim=c(numMod, MyN))
Newtime <- data.frame(Time = seq(from=1/MyStep, to=MyTH, by=1/MyStep), AtRisk = 1)
Newtime$MyId <- 1:dim(Newtime)[1]
Newtime$MyId <- ifelse(Newtime$MyId > follow_up, follow_up, Newtime$MyId)  # Random effects: Using last observed ID for extrapolation
Newtime$EventsL <- 0
Newtime$EventsL[1:follow_up] <- lag(ltHaz$Events)
Newtime$EventsL[1] <- 0
Newtime$EventsL <- ifelse(Newtime$MyId > follow_up, 0, Newtime$EventsL) # AR: Using last observed event count for extrapolation
Newtime$timedelta<-Newtime$Time[2]-Newtime$Time[1]
Newtime$lnTime<-log(Newtime$Time)
# Also have 1x GOF matrix. Rows = Methods, Columns = Method, LL, AIC
dfGOF <- data.frame(matrix(, nrow=19, ncol=4))
colnames(dfGOF) <- c("Model","LnL","Params","AIC")
# Below is constant for when have to derive log-likelihood
llCons <- sum(ltHaz$Events*log(ltHaz$AtRisk) - log(factorial(ltHaz$Events)))
# Names of models to consider
modnames <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1-1","FP1-2",
                 "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")
md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1-1","FP1-2",
      "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")
#=============================================================================================================#
#=============================================================================================================#

#                                       CHE fit and Extrapolation

#=============================================================================================================#
#=============================================================================================================#

####################################
#######    Standard Dist     #######
####################################
MODi <- 1 # Model index
MyDists <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
for (i in 1:7){
  glmTemp <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = sd_bc, dist = MyDists[[i]])
  dfHazEst[MODi,] <- summary(glmTemp, t=Newtime$Time, type="hazard")[[1]]$est/17.4
  ltHaz[[MyDists[[i]]]] <- summary(glmTemp, t=ltHaz$Time, type="hazard")[[1]]$est/17.4
  dfGOF[MODi,1] <- MyDists[[i]]
  dfGOF[MODi,2] <- sum(ltHaz$Events*log(ltHaz[[MyDists[[i]]]]) - ltHaz[[MyDists[[i]]]]*ltHaz$AtRisk) + llCons
  dfGOF[MODi,3] <- glmTemp$npars
  MODi<-MODi+1
}

#########################
#######    FP     #######
#########################

#-----FP1 -----
myLnL <- array(dim=8)
myAIC <- array(dim=8)
MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
for (i in 1:7){
  glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)
  myLnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[i] <- extractAIC(glmTemp)[2]
}
### run for 0
glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)
myLnL[8] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
myAIC[8] <- extractAIC(glmTemp)[2]

FP1res <- data.frame(c("-2","-1","-0.5","0","0.5","1","2","3"))
FP1res <- cbind(FP1res,myLnL,myAIC)
colnames(FP1res) <- c("Powers","LnL","AIC")
FP1res <-arrange(FP1res,AIC)
FP1res[1,]
FP1res[2,]

#-----FP2 -----
myLnL <- array(dim=36)
myAIC <- array(dim=36)
MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
index <- 1
for (i in 1:7){
  for (j in 1:7){
    if (j > i) {
      glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][j])+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)# 
      myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      myAIC[index] <- extractAIC(glmTemp)[2]
      index <- index + 1
    }
  }
}
for (i in 1:7) {
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][i]*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)# 
  myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[index] <- extractAIC(glmTemp)[2]
  index <- index + 1
}

for (i in 1:7) {
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)# 
  myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[index] <- extractAIC(glmTemp)[2]
  index <- index + 1
}

glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)# 
myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
myAIC[index] <- extractAIC(glmTemp)[2]

FP2res <- data.frame(c("-2,-1","-2,-0.5","-2,0.5","-2,1","-2,2","-2,3","-1,-0.5","-1,0.5","-1,1","-1,2","-1,3","-0.5,0.5","-0.5,1","-0.5,2","-0.5,3",
                       "0.5,1","0.5,2","0.5,3","1,2","1,3","2,3","-2,-2","-1,-1","-0.5,-0.5","0.5,0.5","1,1","2,2","3,3",
                       "-2,0","-1,0","-0.5,0","0.5,0","1,0","2,0","3,0","0,0"))
FP2res <- cbind(FP2res,myLnL,myAIC)
colnames(FP2res) <- c("Powers","LnL","AIC")
FP2res <-arrange(FP2res,AIC)
FP2res[1,]
FP2res[2,]



###combine FP results####
#(FP2-best)

modFP1_1 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(3))+ offset(log(timedelta)) , family=binomial(link=cloglog), data=ltHaz)
modFP1_2 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(-0.5))+ offset(log(timedelta)) , family=binomial(link=cloglog), data=ltHaz)
modFP2_1 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(0.5)) + I(Time^(0.5)*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)
modFP2_2 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^0.5) + I(Time^1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)

MODi <- 8
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(modFP1_1)[2] - 2*extractAIC(modFP1_1)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(modFP1_1)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(modFP1_1, newdata=Newtime, type="response") # Extrapolated
ltHaz[modnames[[MODi]]] <- predict(modFP1_1, newdata=ltHaz, type="response")  # Within-sample



MODi<-9
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(modFP1_2)[2] - 2*extractAIC(modFP1_2)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(modFP1_2)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(modFP1_2, newdata=Newtime, type="response") # Extrapolated
ltHaz[modnames[[MODi]]] <- predict(modFP1_2, newdata=ltHaz, type="response")  # Within-sample

MODi<-10
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(modFP2_1)[2] - 2*extractAIC(modFP2_1)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(modFP2_1)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(modFP2_1, newdata=Newtime, type="response") # Extrapolated
ltHaz[modnames[[MODi]]] <- predict(modFP2_1, newdata=ltHaz, type="response")  # Within-sample

MODi<-11
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(modFP2_2)[2] - 2*extractAIC(modFP2_2)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(modFP2_2)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(modFP2_2, newdata=Newtime, type="response") # Extrapolated
ltHaz[modnames[[MODi]]] <- predict(modFP2_2, newdata=ltHaz, type="response")  # Within-sample


##########################
#######    RCS     #######
##########################
MODi <- 12
# First need knot locations for up to 5 internal knots.
# Basing these on equally-spaced percentiles of the observed (uncensored) death times.
bc2 <- subset(che, censrec==1)
myLnL <- array(dim=5)
myAIC <- array(dim=5)
for (i in 1:5){
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz)
  myLnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[i] <- extractAIC(glmTemp)[2]
}
RCSres <- data.frame(c("One","Two","Three","Four","Five"))
RCSres <- cbind(RCSres,myLnL,myAIC)
colnames(RCSres) <- c("Int.Knots","LnL","AIC")
RCSres<-arrange(RCSres)
RCSres[1,]
RCSres[2,]
i<-1
glmTemp1 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz)
i<-2
glmTemp2 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz)

dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(glmTemp1)[2] - 2*extractAIC(glmTemp1)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(glmTemp1)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(glmTemp1, newdata=Newtime, type="response")
ltHaz[modnames[[MODi]]] <- predict(glmTemp1, newdata=ltHaz, type="response")

MODi<-13
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(glmTemp2)[2] - 2*extractAIC(glmTemp2)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(glmTemp2)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(glmTemp2, newdata=Newtime, type="response")
ltHaz[modnames[[MODi]]] <- predict(glmTemp2, newdata=ltHaz, type="response")

#########################
#######    RP     #######
#########################
MODi <- 14
MyAIC <- array(dim=c(6,3))
MyScale <- list("hazard","odds","normal")
for (i in 1:3){
  for (j in 0:5){
    fit<-try(MyTemp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc, k = j, scale = MyScale[[i]]))
    if("try-error" %in% class(fit)) {
      MyAIC[[(i-1)*6+1+j]] <- "error"
    }
    else{
      flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc, k = j, scale = MyScale[[i]])
      MyAIC[[(i-1)*6+1+j]] <- (-2*MyTemp$loglik+2*MyTemp$npars)
    }
  }
}

MyAICResults <- as.data.frame(cbind(seq(1:6)-1,MyAIC))
colnames(MyAICResults) <- c("Int.Knots","Hazard","Odds","Normal")
best_rp_hazard<-data.frame(MyAICResults$Int.Knots,MyAICResults$Hazard)
colnames(best_rp_hazard)<-c("knots","AIC")
best_rp_odds<-data.frame(MyAICResults$Int.Knots,MyAICResults$Odds)
colnames(best_rp_odds)<-c("knots","AIC")
best_rp_normal<-data.frame(MyAICResults$Int.Knots,MyAICResults$Normal)
colnames(best_rp_normal)<-c("knots","AIC")
best_rp_hazard<-arrange(best_rp_hazard,AIC)
best_rp_odds<-arrange(best_rp_odds,AIC)
best_rp_normal<-arrange(best_rp_normal,AIC)

rp_input<-as.data.frame(array(dim=c(6,2)))
rp_input$V1<- c("hazard","hazard","odds","odds","normal","normal")
rp_input$V2<-c(best_rp_hazard[1,1],
               best_rp_hazard[2,1],
               best_rp_odds[1,1],
               best_rp_odds[2,1],
               best_rp_normal[1,1],
               best_rp_normal[2,1])
colnames(rp_input)<-c("scale","knots")
rp_input$knots<-as.numeric(rp_input$knots)

for (i in 1:6) {
  rp_scale<-rp_input[i,1]
  rp_k<-rp_input[i,2]
  rpTemp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc, k = rp_k, scale = rp_scale)
  rpAIC <- (-2*rpTemp$loglik+2*rpTemp$npars)
  dfHazEst[MODi,] <- summary(rpTemp, t=Newtime$Time, type="hazard")[[1]]$est/17.4
  ltHaz[modnames[[MODi]]] <- summary(rpTemp, t=ltHaz$Time, type="hazard")[[1]]$est/17.4
  
  dfGOF[MODi,1] <- modnames[[MODi]]
  dfGOF[MODi,2] <- sum(ltHaz$Events*log(ltHaz[modnames[[MODi]]]) - ltHaz[[modnames[[MODi]]]]*ltHaz$AtRisk) + llCons
  dfGOF[MODi,3] <- rpTemp$npars  
  
  MODi<-MODi+1
}  


#######################################
#             AIC & export            # 
#######################################
dfGOF$AIC <- -2*dfGOF$LnL + 2*dfGOF$Params
write.csv(dfGOF,"che-os-AIC.csv")

##################------------surv rate calculate  ###  lthaz   ###---------------##############
lthaz_plot<-data.frame(ltHaz$Time,ltHaz[,11:29])
dfsurv_lthaz<-as.data.frame(array(dim=c(37,20)))
colnames(dfsurv_lthaz)<-c("Time",modnames)
dfsurv_lthaz$Time<-lthaz_plot$ltHaz.Time
for (i in 1:19) {
  temp_dfhaz<-data.frame(lthaz_plot$ltHaz.Time,lthaz_plot[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurv<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurv_lthaz[,i+1]<-temp_dfsurv$survProp
}

dfFigSurv_lthaz = dfsurv_lthaz %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
##plot##
f_surv2= ggplot() +
  geom_line(data=dfFigSurv_lthaz, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
  geom_line(data = ltHaz,aes(x=Time,y=surv),size=0.8)+
  scale_color_discrete(name="Model")+
  expand_limits(y=c(0,1),x=c(0,2)) + 
  facet_wrap(~Model,nrow=4)+
  scale_x_continuous(breaks = c(seq(from=0, to=2,by = 0.5))) +
  ylab("Overall survival") +
  xlab("Time(years)") +
  guides(color = guide_legend(ncol = 1))  +
  theme(axis.title.x= element_text(size=20,color="black"))+
  theme(axis.title.y= element_text(size=20,color="black"))+
  theme(axis.text.x= element_text(size=16, color="black"))+
  theme(axis.text.y= element_text(size=16, color="black"))+
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 14),strip.text = element_text(size = 14)) 
f_surv2

ggsave("fig_surv_che_os_fit.png",f_surv2,width = 12,height = 8,dpi = 600)


##################------------surv rate calculate  ###  dfhazest   ###---------------##############
#calculate survival over time
dfhar<-t(dfHazEst)
colnames(dfhar) <- dfGOF[,1]
dfhar <- cbind(data.frame(Newtime$Time), dfhar)
dfhar <- rename(dfhar,"Time"="Newtime.Time")

dfsurv<-as.data.frame(array(dim=c(348,20)))
colnames(dfsurv)<-c("Time",modnames)
dfsurv$Time<-dfhar$Time
###  20 model hazard  ###
for (i in 1:19) {
  temp_dfhaz<-data.frame(dfhar$Time,dfhar[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurv<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurv[,i+1]<-temp_dfsurv$survProp
}
write.csv(dfsurv,"surv.refpemche.csv")

best<-flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc, k = 0, scale = "odds")
best$coefficients

