######################################
#         PM learning larvae         #
######################################

########## PM impact on life history traits ##########

#### Survival data
h_t<- read.table("survie_isa.txt",header=TRUE,sep="")
head(h_t)


library(survival)
survdiff(Surv(day_before_death,status)~treatment,data=h_t)
#p-value=0.2 --> no =/= of survival proba bw groups

#### Mortality data
M<-glm(status~treatment,data=h_t,family=binomial) #glm binomial bc var "Status" binomial (0/1)
summary(M)
Anova(M, test="F") #model significative
plot(M)

library(dunn.test)
dunn.test(h_t$status,h_t$treatment)
#=/= bw sub and env ; sub and ctrl

#plot mortality
nb1<-100-(39/39)*100
nb2<-100-(37/37)*100
nb3<-100-(26/56)*100
anova(nb1, nb2, nb3)
matobsnb<-matrix(c(nb1,nb2,nb3),nrow=1,dimnames=list(c("")))

barplot(matobsnb
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Mortality (%)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,110)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) # X axis
axis(side=2,at=c(0,25,50,75,100),cex.axis=1,las=2) # Y axis
abline(h= 0, col = "black") # add horizontal line

text(0.5,50,"a")  #add letter representing significant differences bw groups
text(2.5,50,"a")
text(4.5,60,"b")

########## PM impact on phenotypes ##########

#### Growth data
g<- read.table("size.txt",header=TRUE,sep="")
head(g)
colnames(g)  
names(g)[names(g) == "growth_factor_j1j7.cm."] <- "g_f"

shapiro.test(g$g_f) # do not follow normality

kruskal.test(g$g_f~g$treatment, data=g)
#p-value <0.05 --> =/= bw treatments

dunn.test(g$g_f,g$treatment) #=/= bw sub and ctrl only

#barplot --> better plot type
meanlc1<-mean((g$g_f[g$treatment=="control"])) #get mean
errorlc1<-(sd((g$g_f[g$treatment=="control"])))/(sqrt(length((g$g_f[g$treatment=="control"])))) #get errors

meanlc2<-mean((g$g_f[g$treatment=="environmental"]))
errorlc2<-(sd((g$g_f[g$treatment=="environmental"])))/(sqrt(length((g$g_f[g$treatment=="environmental"]))))

meanlc3<-mean((g$g_f[g$treatment=="sublethal"]))
errorlc3<-(sd((g$g_f[g$treatment=="sublethal"])))/(sqrt(length((g$g_f[g$treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c(""))) #matrix() creates a matrix from the given set of values


barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="growth factor"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,1.26)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) # X axis
axis(side=2,at=c(0,0.25,0.50,0.75,1,1.25),cex.axis=1,las=2) # Y axis
abline(h= 0, col = "black") # add horizontal line

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1.23,"a")  #add letter representing significant differences bw groups
text(2.5,1.22,"a")
text(4.5,1.18,"b")

########## PM impact on behavior ##########

#############################   LIGHT "CUMULATIVE" 
lc<- read.table("light_cumule.txt",header=TRUE,sep="")
head(lc)

#### mean mobility
lc$mobility_mean<-lc$mobility_mean/100

#model 1 : with interaction of treatment and condition
mod1<-glm(mobility_mean~treatment*condition,family="quasibinomial",data=lc) #family quasibiomial bc var bw 0 and 1
summary(mod1)
library(car) 
Anova(mod1,test="F") #test F -->"quasi" family
#no effect of the interaction --> remove it

#model 2 : w/o interaction of treatment and condition
mod2<-glm(mobility_mean~treatment+condition,family="quasibinomial",data=lc) #on enlève l'interaction
summary(mod2)
anova(mod2,mod1,test="F") # no significative =/= bw both model --> keep the simpler model
Anova(mod2,test="F") #treatment effect + condition effect (light/dark)

dunn.test(lc$mobility_mean,lc$treatment)
#/!\ =/= significative bw sub and ctrl and env but p-value =0.000... 

#plots treatment effect on mean mobility
meanlc1<-mean((lc$mobility_mean[lc$treatment=="control"])) #get mean
errorlc1<-(sd((lc$mobility_mean[lc$treatment=="control"])))/(sqrt(length((lc$mobility_mean[lc$treatment=="control"])))) #get errors

meanlc2<-mean((lc$mobility_mean[lc$treatment=="env"]))
errorlc2<-(sd((lc$mobility_mean[lc$treatment=="env"])))/(sqrt(length((lc$mobility_mean[lc$treatment=="env"]))))

meanlc3<-mean((lc$mobility_mean[lc$treatment=="sub"]))
errorlc3<-(sd((lc$mobility_mean[lc$treatment=="sub"])))/(sqrt(length((lc$mobility_mean[lc$treatment=="sub"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c(""))) #matrix() creates a matrix from the given set of values

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Mean mobility"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,0.20)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) # X axis
axis(side=2,at=c(0,0.05,0.10,0.15,0.20),cex.axis=1,las=2) # Y axis
abline(h= 0, col = "black") # add horizontal line

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,0.18,"a")  #add letter representing significant differences bw groups
text(2.5,0.18,"a")
text(4.5,0.13,"b")

#plots condition effect on mean mobility
meanlc11<-mean((lc$mobility_mean[lc$condition=="dark"]))
errorlc11<-(sd((lc$mobility_mean[lc$condition=="dark"])))/(sqrt(length((lc$mobility_mean[lc$condition=="dark"]))))

meanlc12<-mean((lc$mobility_mean[lc$condition=="light"]))
errorlc12<-(sd((lc$mobility_mean[lc$condition=="light"])))/(sqrt(length((lc$mobility_mean[lc$condition=="light"]))))

matobslc1 <-matrix(c(meanlc11,meanlc12),nrow=1,dimnames=list(c("")))

barplot(matobslc1
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Test condition"
        ,ylab="Mean mobility"
        ,cex.lab=1.2
        ,xlim=c(0,4.5)
        ,ylim=c(0,0.20)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0.5,1)
        ,col=c("white","grey"))

axis(side=1.5,at=c(1.5,3.5),labels=c("Light","Dark"),tick=FALSE,cex.axis=1)  
axis(side=2,at=c(0,0.05,0.10,0.15,0.20),cex.axis=1,las=2)  
abline(h= 0, col = "black")


arrows(1.5, meanlc11 - errorlc11, 1.5,meanlc11 + errorlc11, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(3.5, meanlc12 - errorlc12, 3.5,meanlc12 + errorlc12, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(1.5,0.145,"a")
text(3.5,0.17,"b")

#### duration in zone

#model 1 : with interaction of treatment and condition
mod3<-glm(duration_inzone~treatment*condition,family="gaussian",data=lc)
summary(mod3)
Anova(mod3,test="LR")
#no effect of the interaction --> remove it

#model 2 : w/o interaction of treatment and condition
mod4<-glm(duration_inzone~treatment+condition,family="gaussian",data=lc) 
summary(mod4)
anova(mod4,mod3,test="F") # no significative =/= bw both model --> keep the simpler model
Anova(mod4,test="F") #treatment effect but no condition effect (light/dark)

dunn.test(lc$duration_inzone,lc$treatment)
#/!\ =/= significative bw sub and ctrl only

#plots treatment effect on duration_inzone
meanlc1<-mean((lc$duration_inzone[lc$treatment=="control"]))
errorlc1<-(sd((lc$duration_inzone[lc$treatment=="control"])))/(sqrt(length((lc$duration_inzone[lc$treatment=="control"]))))

meanlc2<-mean((lc$duration_inzone[lc$treatment=="env"]))
errorlc2<-(sd((lc$duration_inzone[lc$treatment=="env"])))/(sqrt(length((lc$duration_inzone[lc$treatment=="env"]))))

meanlc3<-mean((lc$duration_inzone[lc$treatment=="sub"]))
errorlc3<-(sd((lc$duration_inzone[lc$treatment=="sub"])))/(sqrt(length((lc$duration_inzone[lc$treatment=="sub"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Duration in zone (sec)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,110)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,10,20,30,40, 50, 60, 70, 80, 90, 100),cex.axis=1,las=2) 
abline(h= 0, col = "black") 

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,60,"a") 
text(2.5,90,"a,b")
text(4.5,108,"b")

#plots condition effect on duration in zone
meanlc11<-mean((lc$duration_inzone[lc$condition=="dark"]))
errorlc11<-(sd((lc$duration_inzone[lc$condition=="dark"])))/(sqrt(length((lc$duration_inzone[lc$condition=="dark"]))))

meanlc12<-mean((lc$duration_inzone[lc$condition=="light"]))
errorlc12<-(sd((lc$duration_inzone[lc$condition=="light"])))/(sqrt(length((lc$duration_inzone[lc$condition=="light"]))))

matobslc1 <-matrix(c(meanlc11,meanlc12),nrow=1,dimnames=list(c("")))

barplot(matobslc1
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Test condition"
        ,ylab="duration in zone (sec)"
        ,cex.lab=1.2
        ,xlim=c(0,4.5)
        ,ylim=c(0,100)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0.5,1)
        ,col=c("white","grey"))

axis(side=1.5,at=c(1.5,3.5),labels=c("Light","Dark"),tick=FALSE,cex.axis=1)  
axis(side=2,at=c(0,10,20, 30, 40, 50, 60, 70, 80),cex.axis=1,las=2)  
abline(h= 0, col = "black")


arrows(1.5, meanlc11 - errorlc11, 1.5,meanlc11 + errorlc11, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(3.5, meanlc12 - errorlc12, 3.5,meanlc12 + errorlc12, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(1.5,95,"a")
text(3.5,75,"a")

#### frequency in zone

#model 1 : with interaction of treatment and condition
mod5<-glm(frequency_inzone~treatment*condition,family="gaussian",data=lc)
summary(mod5)
Anova(mod5,test="LR")
#no effect of the interaction --> remove it

#model 2 : w/o interaction of treatment and condition
mod6<-glm(frequency_inzone~treatment+condition,family="gaussian",data=lc) 
summary(mod6)
anova(mod6,mod5,test="F") # no significative =/= bw both model --> keep the simpler model
Anova(mod6,test="F") #treatment effect + condition effect (light/dark)

dunn.test(lc$frequency_inzone,lc$treatment)
#/!\ =/= significative bw sub and ctrl and env

#plots treatment effect on frequency_inzone
meanlc1<-mean((lc$frequency_inzone[lc$treatment=="control"]))
errorlc1<-(sd((lc$frequency_inzone[lc$treatment=="control"])))/(sqrt(length((lc$frequency_inzone[lc$treatment=="control"]))))

meanlc2<-mean((lc$frequency_inzone[lc$treatment=="env"]))
errorlc2<-(sd((lc$frequency_inzone[lc$treatment=="env"])))/(sqrt(length((lc$frequency_inzone[lc$treatment=="env"]))))

meanlc3<-mean((lc$frequency_inzone[lc$treatment=="sub"]))
errorlc3<-(sd((lc$frequency_inzone[lc$treatment=="sub"])))/(sqrt(length((lc$frequency_inzone[lc$treatment=="sub"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Frequency in zone"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,40)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,10,20,30,40),cex.axis=1,las=2) 
abline(h= 0, col = "black") 

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,35,"a") 
text(2.5,35,"a")
text(4.5,25,"b")

#plots condition effect on frequency in zone
meanlc11<-mean((lc$frequency_inzone[lc$condition=="dark"]))
errorlc11<-(sd((lc$frequency_inzone[lc$condition=="dark"])))/(sqrt(length((lc$frequency_inzone[lc$condition=="dark"]))))

meanlc12<-mean((lc$frequency_inzone[lc$condition=="light"]))
errorlc12<-(sd((lc$frequency_inzone[lc$condition=="light"])))/(sqrt(length((lc$frequency_inzone[lc$condition=="light"]))))

matobslc1 <-matrix(c(meanlc11,meanlc12),nrow=1,dimnames=list(c("")))

barplot(matobslc1
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Test condition"
        ,ylab="frequency in zone"
        ,cex.lab=1.2
        ,xlim=c(0,4.5)
        ,ylim=c(0,40)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0.5,1)
        ,col=c("white","grey"))

axis(side=1.5,at=c(1.5,3.5),labels=c("Light","Dark"),tick=FALSE,cex.axis=1)  
axis(side=2,at=c(0,10,20, 30, 40, 50, 60, 70, 80),cex.axis=1,las=2)  
abline(h= 0, col = "black")


arrows(1.5, meanlc11 - errorlc11, 1.5,meanlc11 + errorlc11, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(3.5, meanlc12 - errorlc12, 3.5,meanlc12 + errorlc12, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(1.5,35,"a")
text(3.5,30,"b")

#### velocity mean

#model 1 : with interaction of treatment and condition
mod7<-glm(velocity_mean~treatment*condition,family="gaussian",data=lc)
summary(mod7)
Anova(mod7,test="LR")
#no effect of the interaction for alpha=0.05 --> remove it

#model 2 : w/o interaction of treatment and condition
mod8<-glm(velocity_mean~treatment+condition,family="gaussian",data=lc) 
summary(mod8)
anova(mod8,mod7,test="F") # no significative =/= bw both model --> keep the simpler model
Anova(mod8,test="F") #treatment effect BUT not condition effect (light/dark)

dunn.test(lc$velocity_mean,lc$treatment)
#/!\ =/= significative bw sub and ctrl and env

#plots treatment effect on velocity_mean
meanlc1<-mean((lc$velocity_mean[lc$treatment=="control"]))
errorlc1<-(sd((lc$velocity_mean[lc$treatment=="control"])))/(sqrt(length((lc$velocity_mean[lc$treatment=="control"]))))

meanlc2<-mean((lc$velocity_mean[lc$treatment=="env"]))
errorlc2<-(sd((lc$velocity_mean[lc$treatment=="env"])))/(sqrt(length((lc$velocity_mean[lc$treatment=="env"]))))

meanlc3<-mean((lc$velocity_mean[lc$treatment=="sub"]))
errorlc3<-(sd((lc$velocity_mean[lc$treatment=="sub"])))/(sqrt(length((lc$velocity_mean[lc$treatment=="sub"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Mean velocity (cm/s)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,0.7)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.10,0.20,0.30,0.40, 0.50),cex.axis=1,las=2) 
abline(h= 0, col = "black") 

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,0.55,"a") 
text(2.5,0.55,"a")
text(4.5,0.4,"b")

#plots condition effect on velocity_mean
meanlc11<-mean((lc$velocity_mean[lc$condition=="dark"]))
errorlc11<-(sd((lc$velocity_mean[lc$condition=="dark"])))/(sqrt(length((lc$velocity_mean[lc$condition=="dark"]))))

meanlc12<-mean((lc$velocity_mean[lc$condition=="light"]))
errorlc12<-(sd((lc$velocity_mean[lc$condition=="light"])))/(sqrt(length((lc$velocity_mean[lc$condition=="light"]))))

matobslc1 <-matrix(c(meanlc11,meanlc12),nrow=1,dimnames=list(c("")))

barplot(matobslc1
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Test condition"
        ,ylab="Mean velocity (cm/s)"
        ,cex.lab=1.2
        ,xlim=c(0,4.5)
        ,ylim=c(0,0.6)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0.5,1)
        ,col=c("white","grey"))

axis(side=1.5,at=c(1.5,3.5),labels=c("Light","Dark"),tick=FALSE,cex.axis=1)  
axis(side=2,at=c(0,0.10,0.20, 0.30, 0.40, 0.50),cex.axis=1,las=2)  
abline(h= 0, col = "black")


arrows(1.5, meanlc11 - errorlc11, 1.5,meanlc11 + errorlc11, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(3.5, meanlc12 - errorlc12, 3.5,meanlc12 + errorlc12, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(1.5,0.45,"a")
text(3.5,0.50,"a")

#### velocity max
#model 1 : with interaction of treatment and condition
mod9<-glm(velocity_max~treatment*condition,family="gaussian",data=lc)
summary(mod9)
Anova(mod9,test="LR")
#no effect of the interaction for alpha=0.05 --> remove it

#model 2 : w/o interaction of treatment and condition
mod10<-glm(velocity_max~treatment+condition,family="gaussian",data=lc) 
summary(mod10)
anova(mod10,mod9,test="F") # no significative =/= bw both model --> keep the simpler model
Anova(mod10,test="F") #treatment effect BUT not condition effect (light/dark)

dunn.test(lc$velocity_max,lc$treatment)
#/!\ =/= significative bw sub and ctrl and env

#plots treatment effect on velocity_mean
meanlc1<-mean((lc$velocity_max[lc$treatment=="control"]))
errorlc1<-(sd((lc$velocity_max[lc$treatment=="control"])))/(sqrt(length((lc$velocity_max[lc$treatment=="control"]))))

meanlc2<-mean((lc$velocity_max[lc$treatment=="env"]))
errorlc2<-(sd((lc$velocity_max[lc$treatment=="env"])))/(sqrt(length((lc$velocity_max[lc$treatment=="env"]))))

meanlc3<-mean((lc$velocity_max[lc$treatment=="sub"]))
errorlc3<-(sd((lc$velocity_max[lc$treatment=="sub"])))/(sqrt(length((lc$velocity_max[lc$treatment=="sub"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Maximum velocity (cm/s)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,9)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,1,2,3,4,5,6,7,8),cex.axis=1,las=2) 
abline(h= 0, col = "black") 

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,8,"a") 
text(2.5,8,"a")
text(4.5,6,"b")

#plots condition effect on velocity_meax
meanlc11<-mean((lc$velocity_max[lc$condition=="dark"]))
errorlc11<-(sd((lc$velocity_max[lc$condition=="dark"])))/(sqrt(length((lc$velocity_max[lc$condition=="dark"]))))

meanlc12<-mean((lc$velocity_max[lc$condition=="light"]))
errorlc12<-(sd((lc$velocity_max[lc$condition=="light"])))/(sqrt(length((lc$velocity_max[lc$condition=="light"]))))

matobslc1 <-matrix(c(meanlc11,meanlc12),nrow=1,dimnames=list(c("")))

barplot(matobslc1
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Test condition"
        ,ylab="Maximum velocity (cm/s)"
        ,cex.lab=1.2
        ,xlim=c(0,4.5)
        ,ylim=c(0,8)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0.5,1)
        ,col=c("white","grey"))

axis(side=1.5,at=c(1.5,3.5),labels=c("Light","Dark"),tick=FALSE,cex.axis=1)  
axis(side=2,at=c(0,1,2,3,4,5,6,7),cex.axis=1,las=2)  
abline(h= 0, col = "black")


arrows(1.5, meanlc11 - errorlc11, 1.5,meanlc11 + errorlc11, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(3.5, meanlc12 - errorlc12, 3.5,meanlc12 + errorlc12, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(1.5,7,"a")
text(3.5,7,"a")

#### Acceleration max
#model 1 : with interaction of treatment and condition
mod11<-glm(acceleration_max~treatment*condition,family="gaussian",data=lc)
summary(mod9)
Anova(mod11,test="LR")
#no effect of the interaction for alpha=0.05 --> remove it

#model 2 : w/o interaction of treatment and condition
mod12<-glm(acceleration_max~treatment+condition,family="gaussian",data=lc) 
summary(mod10)
anova(mod12,mod11,test="F") # no significative =/= bw both model --> keep the simpler model
Anova(mod12,test="F") #treatment effect BUT not condition effect (light/dark)

dunn.test(lc$acceleration_max,lc$treatment)
#/!\ =/= significative bw sub and ctrl and env

#plots treatment effect on velocity_mean
meanlc1<-mean((lc$acceleration_max[lc$treatment=="control"]))
errorlc1<-(sd((lc$acceleration_max[lc$treatment=="control"])))/(sqrt(length((lc$acceleration_max[lc$treatment=="control"]))))

meanlc2<-mean((lc$acceleration_max[lc$treatment=="env"]))
errorlc2<-(sd((lc$acceleration_max[lc$treatment=="env"])))/(sqrt(length((lc$acceleration_max[lc$treatment=="env"]))))

meanlc3<-mean((lc$acceleration_max[lc$treatment=="sub"]))
errorlc3<-(sd((lc$acceleration_max[lc$treatment=="sub"])))/(sqrt(length((lc$acceleration_max[lc$treatment=="sub"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Accelaration maximum (cm/s²)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,125)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,10,20,30,40,50,60,70,80,90,100,110,120),cex.axis=1,las=2) 
abline(h= 0, col = "black") 

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,122,"a") 
text(2.5,115,"a")
text(4.5,80,"b")

#plots condition effect on acceleration_max
meanlc11<-mean((lc$acceleration_max[lc$condition=="dark"]))
errorlc11<-(sd((lc$acceleration_max[lc$condition=="dark"])))/(sqrt(length((lc$acceleration_max[lc$condition=="dark"]))))

meanlc12<-mean((lc$acceleration_max[lc$condition=="light"]))
errorlc12<-(sd((lc$acceleration_max[lc$condition=="light"])))/(sqrt(length((lc$acceleration_max[lc$condition=="light"]))))

matobslc1 <-matrix(c(meanlc11,meanlc12),nrow=1,dimnames=list(c("")))

barplot(matobslc1
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Test condition"
        ,ylab="Maximum acceleration (cm/s²)"
        ,cex.lab=1.2
        ,xlim=c(0,4.5)
        ,ylim=c(0,120)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0.5,1)
        ,col=c("white","grey"))

axis(side=1.5,at=c(1.5,3.5),labels=c("Light","Dark"),tick=FALSE,cex.axis=1)  
axis(side=2,at=c(0,10,20,30,40,50,60,70,80,90,100,110),cex.axis=1,las=2)  
abline(h= 0, col = "black")


arrows(1.5, meanlc11 - errorlc11, 1.5,meanlc11 + errorlc11, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(3.5, meanlc12 - errorlc12, 3.5,meanlc12 + errorlc12, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(1.5,110,"a")
text(3.5,105,"a")

####Total distance swimmed
#model 1 : with interaction of treatment and condition
mod13<-glm(relative_tot_dist~treatment*condition,family="gaussian",data=lc)
summary(mod13)
Anova(mod13,test="LR")
#no effect of the interaction for alpha=0.05 --> remove it

#model 2 : w/o interaction of treatment and condition
mod14<-glm(relative_tot_dist~treatment+condition,family="gaussian",data=lc) 
summary(mod14)
anova(mod14,mod13,test="F") # no significative =/= bw both model --> keep the simpler model
Anova(mod14,test="F") #treatment effect BUT not condition effect (light/dark)

dunn.test(lc$relative_tot_dist,lc$treatment)
#/!\ =/= significative bw sub and ctrl and env

#plots treatment effect on relative_tot_dist
meanlc1<-mean((lc$relative_tot_dist[lc$treatment=="control"]))
errorlc1<-(sd((lc$relative_tot_dist[lc$treatment=="control"])))/(sqrt(length((lc$relative_tot_dist[lc$treatment=="control"]))))

meanlc2<-mean((lc$relative_tot_dist[lc$treatment=="env"]))
errorlc2<-(sd((lc$relative_tot_dist[lc$treatment=="env"])))/(sqrt(length((lc$relative_tot_dist[lc$treatment=="env"]))))

meanlc3<-mean((lc$relative_tot_dist[lc$treatment=="sub"]))
errorlc3<-(sd((lc$relative_tot_dist[lc$treatment=="sub"])))/(sqrt(length((lc$relative_tot_dist[lc$treatment=="sub"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="relative total distance swim (cm)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,60)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,10,20,30,40,50,60),cex.axis=1,las=2) 
abline(h= 0, col = "black") 

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,55,"a") 
text(2.5,55,"a")
text(4.5,40,"b")

#plots condition effect on relative_tot_dist
meanlc11<-mean((lc$relative_tot_dist[lc$condition=="dark"]))
errorlc11<-(sd((lc$relative_tot_dist[lc$condition=="dark"])))/(sqrt(length((lc$relative_tot_dist[lc$condition=="dark"]))))

meanlc12<-mean((lc$relative_tot_dist[lc$condition=="light"]))
errorlc12<-(sd((lc$relative_tot_dist[lc$condition=="light"])))/(sqrt(length((lc$relative_tot_dist[lc$condition=="light"]))))

matobslc1 <-matrix(c(meanlc11,meanlc12),nrow=1,dimnames=list(c("")))

barplot(matobslc1
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Test condition"
        ,ylab="relative total distance swim (cm)"
        ,cex.lab=1.2
        ,xlim=c(0,4.5)
        ,ylim=c(0,50)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0.5,1)
        ,col=c("white","grey"))

axis(side=1.5,at=c(1.5,3.5),labels=c("Light","Dark"),tick=FALSE,cex.axis=1)  
axis(side=2,at=c(0,10,20,30,40,50),cex.axis=1,las=2)  
abline(h= 0, col = "black")


arrows(1.5, meanlc11 - errorlc11, 1.5,meanlc11 + errorlc11, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(3.5, meanlc12 - errorlc12, 3.5,meanlc12 + errorlc12, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(1.5,45,"a")
text(3.5,48,"a")

########## PM impact on gene expression ##########
####GSS
GSS<-read.table("GSS.txt", header=T, sep="", fill=T)
head(GSS)

shapiro.test(GSS$GE)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(GE~Treatment, data=GSS) #no effect of treatment on GSS gene expression


#plot results
meanlc1<-mean((GSS$GE[GSS$Treatment=="control"]))
errorlc1<-(sd((GSS$GE[GSS$Treatment=="control"])))/(sqrt(length((GSS$GE[GSS$Treatment=="control"]))))

meanlc2<-mean((GSS$GE[GSS$Treatment=="environmental"]))
errorlc2<-(sd((GSS$GE[GSS$Treatment=="environmental"])))/(sqrt(length((GSS$GE[GSS$Treatment=="environmental"]))))

meanlc3<-mean((GSS$GE[GSS$Treatment=="sublethal"]))
errorlc3<-(sd((GSS$GE[GSS$Treatment=="sublethal"])))/(sqrt(length((GSS$GE[GSS$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="GSS gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,1.8)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.4,0.8,1.2, 1.6),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1,"a")  
text(2.5,1.6,"a")
text(4.5,1.7,"a")

####MAOA
MAOA<-read.table("MAOA.txt", header=T, sep="", fill=T)
MAOA<- na.omit(MAOA)
head(MAOA)

shapiro.test(MAOA$GE)  #follow normal distrib

#anova
M1<-aov(GE~Treatment, data=MAOA) 
summary(M1) #no effect of treatment

#plot results
meanlc1<-mean((MAOA$GE[MAOA$Treatment=="control"]))
errorlc1<-(sd((MAOA$GE[MAOA$Treatment=="control"])))/(sqrt(length((MAOA$GE[MAOA$Treatment=="control"]))))

meanlc2<-mean((MAOA$GE[MAOA$Treatment=="environmental"]))
errorlc2<-(sd((MAOA$GE[MAOA$Treatment=="environmental"])))/(sqrt(length((MAOA$GE[MAOA$Treatment=="environmental"]))))

meanlc3<-mean((MAOA$GE[MAOA$Treatment=="sublethal"]))
errorlc3<-(sd((MAOA$GE[MAOA$Treatment=="sublethal"])))/(sqrt(length((MAOA$GE[MAOA$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="MAOA gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,1.3)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.4,0.8,1.2),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,0.9,"a")  
text(2.5,1.15,"a")
text(4.5,1.22,"a")

####Mecp2
mecp2<-read.table("Mecp2.txt", header=T, sep="", fill=T)
mecp2<- na.omit(mecp2)
head(mecp2)

shapiro.test(mecp2$GE)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(GE~Treatment, data=mecp2) #no effect of treatment


#plot results
meanlc1<-mean((mecp2$GE[mecp2$Treatment=="control"]))
errorlc1<-(sd((mecp2$GE[mecp2$Treatment=="control"])))/(sqrt(length((mecp2$GE[mecp2$Treatment=="control"]))))

meanlc2<-mean((mecp2$GE[mecp2$Treatment=="environmental"]))
errorlc2<-(sd((mecp2$GE[mecp2$Treatment=="environmental"])))/(sqrt(length((mecp2$GE[mecp2$Treatment=="environmental"]))))

meanlc3<-mean((mecp2$GE[mecp2$Treatment=="sublethal"]))
errorlc3<-(sd((mecp2$GE[mecp2$Treatment=="sublethal"])))/(sqrt(length((mecp2$GE[mecp2$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Mecp2 gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,2.5)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.5,1,1.5,2, 2.5),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,2.1,"a")  
text(2.5,1.6,"a")
text(4.5,2.4,"a")

####Nipbl
nipbl<-read.table("Nipbl.txt", header=T, sep="", fill=T)
nipbl<- na.omit(nipbl)
head(nipbl)

shapiro.test(nipbl$GE)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(GE~Treatment, data=nipbl) #no effect of treatment


#plot results
meanlc1<-mean((nipbl$GE[nipbl$Treatment=="control"]))
errorlc1<-(sd((nipbl$GE[nipbl$Treatment=="control"])))/(sqrt(length((nipbl$GE[nipbl$Treatment=="control"]))))

meanlc2<-mean((nipbl$GE[nipbl$Treatment=="environmental"]))
errorlc2<-(sd((nipbl$GE[nipbl$Treatment=="environmental"])))/(sqrt(length((nipbl$GE[nipbl$Treatment=="environmental"]))))

meanlc3<-mean((nipbl$GE[nipbl$Treatment=="sublethal"]))
errorlc3<-(sd((nipbl$GE[nipbl$Treatment=="sublethal"])))/(sqrt(length((nipbl$GE[nipbl$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Nipbl gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,1.8)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.3,0.6,0.9,1.2,1.5,1.8),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1.75,"a")  
text(2.5,1.2,"a")
text(4.5,1.6,"a")

#Dnmt3a
db3a<-read.table("D3a.txt", header=T, sep="", fill=T)
db3a<- na.omit(db3a)
head(db3a)

shapiro.test(db3a$GE)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(GE~Treatment, data=db3a) #no effect of treatment


#plot results
meanlc1<-mean((db3a$GE[db3a$Treatment=="control"]))
errorlc1<-(sd((db3a$GE[db3a$Treatment=="control"])))/(sqrt(length((db3a$GE[db3a$Treatment=="control"]))))

meanlc2<-mean((db3a$GE[db3a$Treatment=="environmental"]))
errorlc2<-(sd((db3a$GE[db3a$Treatment=="environmental"])))/(sqrt(length((db3a$GE[db3a$Treatment=="environmental"]))))

meanlc3<-mean((db3a$GE[db3a$Treatment=="sublethal"]))
errorlc3<-(sd((db3a$GE[db3a$Treatment=="sublethal"])))/(sqrt(length((db3a$GE[db3a$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Dnmt3a gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,2.5)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.5,1,1.5,2,2.5),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1.4,"a")  
text(2.5,2.3,"a")
text(4.5,1.8,"a")

####Tollip
tollip<-read.table("Tollip.txt", header=T, sep="", fill=T)
tollip<- na.omit(tollip)
head(tollip)

shapiro.test(tollip$GE)  #follow normal distrib

M2<-aov(GE~Treatment, data=tollip) 
summary(M2) #no effect of the treatment

#plot results
meanlc1<-mean((tollip$GE[tollip$Treatment=="control"]))
errorlc1<-(sd((tollip$GE[tollip$Treatment=="control"])))/(sqrt(length((tollip$GE[tollip$Treatment=="control"]))))

meanlc2<-mean((tollip$GE[tollip$Treatment=="environmental"]))
errorlc2<-(sd((tollip$GE[tollip$Treatment=="environmental"])))/(sqrt(length((tollip$GE[tollip$Treatment=="environmental"]))))

meanlc3<-mean((tollip$GE[tollip$Treatment=="sublethal"]))
errorlc3<-(sd((tollip$GE[tollip$Treatment=="sublethal"])))/(sqrt(length((tollip$GE[tollip$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Tollip gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,1.4)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.4,0.8,1.2),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1.2,"a")  
text(2.5,1.3,"a")
text(4.5,1.2,"a")

####DB3ba
db3ba<-read.table("Dnmt3ba.txt", header=T, sep="", fill=T)
db3ba<- na.omit(db3ba)
head(db3ba)

shapiro.test(db3ba$GE)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(GE~Treatment, data=db3ba)

#plot results
meanlc1<-mean((db3ba$GE[db3ba$Treatment=="control"]))
errorlc1<-(sd((db3ba$GE[db3ba$Treatment=="control"])))/(sqrt(length((db3ba$GE[db3ba$Treatment=="control"]))))

meanlc2<-mean((db3ba$GE[db3ba$Treatment=="environmental"]))
errorlc2<-(sd((db3ba$GE[db3ba$Treatment=="environmental"])))/(sqrt(length((db3ba$GE[db3ba$Treatment=="environmental"]))))

meanlc3<-mean((db3ba$GE[db3ba$Treatment=="sublethal"]))
errorlc3<-(sd((db3ba$GE[db3ba$Treatment=="sublethal"])))/(sqrt(length((db3ba$GE[db3ba$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Dnmt3ba gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,1.3)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.4,0.8,1.2),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1.17,"a")  
text(2.5,1,"a")
text(4.5,1.17,"a")


########## PM impact on gene methylation ##########
####Nipbl
pyronipbl<-read.table("nipbl-pyro.txt", header=T, sep="", fill=T)
head(pyronipbl)

shapiro.test(pyronipbl$methylation)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(methylation~Treatment, data=pyronipbl)

#plot results
meanlc1<-mean((pyronipbl$methylation[pyronipbl$Treatment=="control"]))
errorlc1<-(sd((pyronipbl$methylation[pyronipbl$Treatment=="control"])))/(sqrt(length((pyronipbl$methylation[pyronipbl$Treatment=="control"]))))

meanlc2<-mean((pyronipbl$methylation[pyronipbl$Treatment=="environmental"]))
errorlc2<-(sd((pyronipbl$methylation[pyronipbl$Treatment=="environmental"])))/(sqrt(length((pyronipbl$methylation[pyronipbl$Treatment=="environmental"]))))

meanlc3<-mean((pyronipbl$methylation[pyronipbl$Treatment=="sublethal"]))
errorlc3<-(sd((pyronipbl$methylation[pyronipbl$Treatment=="sublethal"])))/(sqrt(length((pyronipbl$methylation[pyronipbl$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Nipbl methylation (%)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,33)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,10,20,30),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,27,"a")  
text(2.5,25,"a")
text(4.5,32,"a")

####Dnmt3a
pyrod3a<-read.table("D3a-pyro.txt", header=T, sep="", fill=T)
head(pyrod3a)

shapiro.test(pyrod3a$methylation)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(methylation~Treatment, data=pyrod3a)

#plot results
meanlc1<-mean((pyrod3a$methylation[pyrod3a$Treatment=="control"]))
errorlc1<-(sd((pyrod3a$methylation[pyrod3a$Treatment=="control"])))/(sqrt(length((pyrod3a$methylation[pyrod3a$Treatment=="control"]))))

meanlc2<-mean((pyrod3a$methylation[pyrod3a$Treatment=="environmental"]))
errorlc2<-(sd((pyrod3a$methylation[pyrod3a$Treatment=="environmental"])))/(sqrt(length((pyrod3a$methylation[pyrod3a$Treatment=="environmental"]))))

meanlc3<-mean((pyrod3a$methylation[pyrod3a$Treatment=="sublethal"]))
errorlc3<-(sd((pyrod3a$methylation[pyrod3a$Treatment=="sublethal"])))/(sqrt(length((pyrod3a$methylation[pyrod3a$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Dnmt3a methylation (%)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,57)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,15,25,35, 45, 55),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,54,"a")  
text(2.5,50,"a")
text(4.5,55,"a")

####GSS
pyrogss<-read.table("gss-pyro.txt", header=T, sep="", fill=T)
head(pyrogss)

shapiro.test(pyrogss$methylation)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(methylation~Treatment, data=pyrogss)

#plot results
meanlc1<-mean((pyrogss$methylation[pyrogss$Treatment=="control"]))
errorlc1<-(sd((pyrogss$methylation[pyrogss$Treatment=="control"])))/(sqrt(length((pyrogss$methylation[pyrogss$Treatment=="control"]))))

meanlc2<-mean((pyrogss$methylation[pyrogss$Treatment=="environmental"]))
errorlc2<-(sd((pyrogss$methylation[pyrogss$Treatment=="environmental"])))/(sqrt(length((pyrogss$methylation[pyrogss$Treatment=="environmental"]))))

meanlc3<-mean((pyrogss$methylation[pyrogss$Treatment=="sublethal"]))
errorlc3<-(sd((pyrogss$methylation[pyrogss$Treatment=="sublethal"])))/(sqrt(length((pyrogss$methylation[pyrogss$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="GSS methylation (%)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,96)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,10,20,30,40,50,60,70,80,90),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,93,"a")  
text(2.5,94,"a")
text(4.5,93,"a")

####Mecp2
pyromecp2<-read.table("mecp2-pyro.txt", header=T, sep="", fill=T)
head(pyromecp2)

shapiro.test(pyromecp2$methylation)  #do not follow normal distrib

#let's use kruskal wallis
kruskal.test(methylation~Treatment, data=pyromecp2)

#plot results
meanlc1<-mean((pyromecp2$methylation[pyromecp2$Treatment=="control"]))
errorlc1<-(sd((pyromecp2$methylation[pyromecp2$Treatment=="control"])))/(sqrt(length((pyromecp2$methylation[pyromecp2$Treatment=="control"]))))

meanlc2<-mean((pyromecp2$methylation[pyromecp2$Treatment=="environmental"]))
errorlc2<-(sd((pyromecp2$methylation[pyromecp2$Treatment=="environmental"])))/(sqrt(length((pyromecp2$methylation[pyromecp2$Treatment=="environmental"]))))

meanlc3<-mean((pyromecp2$methylation[pyromecp2$Treatment=="sublethal"]))
errorlc3<-(sd((pyromecp2$methylation[pyromecp2$Treatment=="sublethal"])))/(sqrt(length((pyromecp2$methylation[pyromecp2$Treatment=="sublethal"]))))

matobslc0 <-matrix(c(meanlc1,meanlc2,meanlc3),nrow=1,dimnames=list(c("")))

barplot(matobslc0
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (µg/L)"
        ,ylab="Mecp2 methylation (%)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,6)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,1.5, 2.5, 3.5, 4.5, 5.5),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc1 - errorlc1, 0.5,meanlc1 + errorlc1, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc2 - errorlc2, 2.5,meanlc2 + errorlc2, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc3 - errorlc3, 4.5,meanlc3 + errorlc3, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,4.7,"a")  
text(2.5,5.4,"a")
text(4.5,5.8,"a")
