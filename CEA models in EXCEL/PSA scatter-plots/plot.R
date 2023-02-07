library("ggplot2")
US <- read.csv("US.csv")
China<- read.csv("China.csv")

US$Group<-factor(US$Group,levels=c("0","200","400","600","800","1000","1200","1400","1600","1800","2000"))

summary(US)

wtp1<-100000
wtp2<-150000
y1<-function(x){wtp1*x}
y2<-function(x){wtp2*x}
ros<-ggplot(data = US)+geom_point(aes(x=Delta_QALYs,y=Delta_Costs,group=Group,colour=Group))+
  scale_color_discrete(name="Price of Serplulimab")+
    coord_cartesian(ylim=c(0,160000),xlim = c(0,1))+
    stat_function(fun = y1,geom = "line",xlim = c(0,1),size=1,linetype="dashed",colour="black")+
    stat_function(fun = y2,geom = "line",xlim = c(0,1),size=1,linetype="dashed",colour="black")+
  geom_text(aes(x= 0.85, y = 70000, label = "WTP = 100,000"),size=7)+
  geom_text(aes(x= 0.85, y = 110000, label = "WTP = 150,000"),size=7)+
    labs(x="Incremental QALYs",y="Incremental Costs", title = "Scatter plot of PSA (US)")+
  theme(axis.title.x= element_text(size=20,color="black"))+
  theme(axis.title.y= element_text(size=20,color="black"))+
  theme(panel.grid= element_blank())+
  theme(axis.line = element_line(colour = "black",size = 1))+
  theme(axis.text.x= element_text(size=16, color="black"))+
  theme(axis.text.y= element_text(size=16, color="black"))+
  theme(title = element_text(size=20, color="black"))+
  theme(panel.background = element_blank())+
  theme(panel.border=element_blank())+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),legend.title = element_text(size = 14),strip.text = element_text(size = 14)) +
  guides(color = guide_legend(nrow = 2))  
  
ros
ggsave("PSA_us.png",ros,width = 18,height = 12,dpi = 300)


summary(China)
China$Group<-factor(China$Group,levels=c("0","200","400","600","800","1000","1200","1400","1600","1800","2000"))

wtp1<-12728
wtp2<-38184
y1<-function(x){wtp1*x}
y2<-function(x){wtp2*x}
ros<-ggplot(data = China)+geom_point(aes(x=Delta_QALYs,y=Delta_Costs,group=Group,colour=Group))+
  scale_color_discrete(name="Price of Serplulimab")+
  coord_cartesian(ylim=c(0,160000),xlim = c(0,1))+
  stat_function(fun = y1,geom = "line",xlim = c(0,1),size=1,linetype="dashed",colour="black")+
  stat_function(fun = y2,geom = "line",xlim = c(0,1),size=1,linetype="dashed",colour="black")+
  geom_text(aes(x= 0.85, y = 5000, label = "WTP = 12,728"),size=7)+
  geom_text(aes(x= 0.85, y = 25000, label = "WTP = 38,184"),size=7)+
  labs(x="Incremental QALYs",y="Incremental Costs", title = "Scatter plot of PSA (China)")+
  theme(axis.title.x= element_text(size=20,color="black"))+
  theme(axis.title.y= element_text(size=20,color="black"))+
  theme(panel.grid= element_blank())+
  theme(axis.line = element_line(colour = "black",size = 1))+
  theme(axis.text.x= element_text(size=16, color="black"))+
  theme(axis.text.y= element_text(size=16, color="black"))+
  theme(title = element_text(size=20, color="black"))+
  theme(panel.background = element_blank())+
  theme(panel.border=element_blank())+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),legend.title = element_text(size = 14),strip.text = element_text(size = 14)) +
  guides(color = guide_legend(nrow = 2))  

ros
ggsave("PSA_ch.png",ros,width = 18,height = 12,dpi = 300)

