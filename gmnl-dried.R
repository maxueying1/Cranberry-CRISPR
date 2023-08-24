############################################################################
#############################   TREATMENT 1   ##############################
############################################################################
library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

##read data
setwd("D:/RA/Cranberry-gene/data clean/for R")
mydata_tr1<-read.csv('Dried_Total_1.csv', fileEncoding = 'UTF-8-BOM')
head(mydata_tr1)
##mean scale
mean(mydata_tr1$Scale)
sum(mydata_tr1$Scale <= 7, na.rm=TRUE)
sum(mydata_tr1$Scale < 7, na.rm=TRUE)

##delete unrelated columns
mydata_tr1$cset <- NULL
mydata_tr1$Block <- NULL
mydata_tr1$Option1 <- NULL
mydata_tr1$Option2 <- NULL
mydata_tr1$Option3 <- NULL

##generate certainty dummy
mydata_tr1$certain=c()
for (i in 1:4500){
  if (mydata_tr1$Scale[i]>=7){mydata_tr1$certain[i]=1}
  else if (mydata_tr1$Scale[i]<7){mydata_tr1$certain[i]=0}
} 
##check uncertainty dummy
mydata_tr1$uncertain=c()
for (i in 1:4500){
  if (mydata_tr1$Scale[i]<7){mydata_tr1$uncertain[i]=1}
  else if (mydata_tr1$Scale[i]>=7){mydata_tr1$uncertain[i]=0}
} 
##check 
mean(mydata_tr1$certain)
sum(mydata_tr1$Scale>=7)/4500

##Generate a column named id, here id means group, which needs to be distinct 
##across all individuals, not just within each individual.
##This id(group) column will be identified as "chid.id" in the gmnl command.
mydata_tr1$id=c()
for (i in 1:4500){
  mydata_tr1$id[i]=ceiling(i/3)
}
##change name of alternatives
for (i in 1:4500){
  if (mydata_tr1$altern[i]==1){mydata_tr1$altern[i]="Option1"}
  else if (mydata_tr1$altern[i]==2){mydata_tr1$altern[i]="Option2"}
  else if (mydata_tr1$altern[i]==3){mydata_tr1$altern[i]="Option3"}}
##view the cleaned data
head(mydata_tr1)

##labeled the alternative with reduced sugar content as 1
mydata_tr1$Reduced=c()
for (i in 1:4500){
  if (mydata_tr1$TSugar[i]==14){mydata_tr1$Reduced[i]=1}
  else if (mydata_tr1$TSugar[i]!=14){mydata_tr1$Reduced[i]=0}
}
head(mydata_tr1)

##labeled the alternative with regular sugar content as 1
mydata_tr1$Regular=c()
for (i in 1:4500){
  if (mydata_tr1$TSugar[i]==29){mydata_tr1$Regular[i]=1}
  else if (mydata_tr1$TSugar[i]!=29){mydata_tr1$Regular[i]=0}
}
head(mydata_tr1)

##generate bland flavor variable
mydata_tr1$BlandFlavor=1-mydata_tr1$FullFlav
for (i in 1:1500){
  mydata_tr1$BlandFlavor[3*i]=0
}
head(mydata_tr1)

mean(mydata_tr1$FullFlav)
mean(mydata_tr1$BlandFlavor)

###################################################################################
##########################   PREFERENCE-SPACE MODELS    ###  MIXED LOGIT  #########
###################################################################################

##first, change dataset into gmnl format by using mlogit.data
##this works
TM_tr1 <- mlogit.data(mydata_tr1,choice="Choice",shape="long",
                      alt.var="altern",
                      chid.var="id",id.var="PersonID")
head(TM_tr1)

#####################################################################
###########################  Simple MNL  ############################
#####################################################################

mnl_tr1 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                data=TM_tr1,model="mnl")
summary(mnl_tr1)
getSummary.gmnl(mnl_tr1,alpha=0.05)
##wtp
-1*wtp.gmnl(mnl_tr1,wrt="price")

#####################################################################
#######################  RPL or MIXED LOGIT  ########################
#####################################################################

##Bopt is random, correlated parameters
rplc_or_tr1 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                    data=TM_tr1,
                    ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                    R=500, 
                    #halton=NA,
                    panel=TRUE,
                    correlation=TRUE,
                    model="mixl")
summary(rplc_or_tr1)
-1*wtp.gmnl(rplc_or_tr1,wrt="price")
cor.gmnl(rplc_or_tr1)
getSummary.gmnl(rplc_or_tr1,alpha=0.05)

###################################################################################
###############################    WTP-SPACE MODELS    ############################
###################################################################################

##models in wtp space start with a model in preference space
##then make inference on the appropriate ratio that represents wtp

##first we need to compute the negative of the price attribute
wtp_data_tr1<-mlogit.data(mydata_tr1,choice="Choice",shape="long",alt.var="altern",
                          chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data_tr1)

#######################################################################
#########################  G-MNL-II  ##################################
#######################################################################

####Bopt is random
wtp_rplc_or_tr1<- -1*wtp.gmnl(rplc_or_tr1,wrt="price")
wtp_rplc_or_tr1<-data.frame(wtp_rplc_or_tr1)
start_or_tr1_1 <-c(1,wtp_rplc_or_tr1$Estimate[1:4],rep(0.1,11),0.2,0)
start_or_tr1_2 <- c(1,wtp_rplc_or_tr1$Estimate[1:4],rep(0.1,11),0,0)
start_or_tr1_3 <- c(1,rep(0.1,15),0.1,0)
start_or_tr1_4 <- c(1,rep(1,15),1,0)
start_or_tr1_5 <- c(1,rep(-1,15),-1,0)
start_tr1=rbind(start_or_tr1_1,start_or_tr1_2,start_or_tr1_3,start_or_tr1_4,start_or_tr1_5)

wtpc_or_tr1 <-list()
ll_tr1<- list()
for (i in 1:5){
  wtpc_or_tr1[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt
                           |0|0|0|certain-1,
                           data=wtp_data_tr1,
                           ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                           R=500,
                           print.init = TRUE, 
                           panel=TRUE,
                           start=start_tr1[i,],
                           fixed=c(TRUE,rep(FALSE,16),TRUE), 
                           model="gmnl",
                           correlation=TRUE,
                           notscale=c(0,0,0,0,1),
                           method="bhhh",
                           iterlim=500)
  ll_tr1[[i]]=getSummary.gmnl(wtpc_or_tr1[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_tr1

summary(wtpc_or_tr1[[5]])
cor.gmnl(wtpc_or_tr1[[2]])
getSummary.gmnl(wtpc_or_tr1[[2]],alpha=0.05)


#---------------------------------------------------------------------------
#-----------------------------  TREATMENT 2  -------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")
mydata_tr2<-read.csv('Dried_Total_2.csv', fileEncoding = 'UTF-8-BOM')
head(mydata_tr2)

mean(mydata_tr2$Scale)
sum(mydata_tr2$Scale <= 7, na.rm=TRUE)
sum(mydata_tr2$Scale < 7, na.rm=TRUE)

mydata_tr2$cset <- NULL
mydata_tr2$Block <- NULL
mydata_tr2$Option1 <- NULL
mydata_tr2$Option2 <- NULL
mydata_tr2$Option3 <- NULL

##generate certainty dummy
mydata_tr2$certain=c()
for (i in 1:4500){
  if (mydata_tr2$Scale[i]>=7){mydata_tr2$certain[i]=1}
  else if (mydata_tr2$Scale[i]<7){mydata_tr2$certain[i]=0}
} 
##check uncertainty dummy
mydata_tr2$uncertain=c()
for (i in 1:4500){
  if (mydata_tr2$Scale[i]<7){mydata_tr2$uncertain[i]=1}
  else if (mydata_tr2$Scale[i]>=7){mydata_tr2$uncertain[i]=0}
} 
##check 
mean(mydata_tr2$certain)
sum(mydata_tr2$Scale>=7)/4500

mydata_tr2$id=c()
for (i in 1:4500){
  mydata_tr2$id[i]=ceiling(i/3)
}

for (i in 1:4500){
  if (mydata_tr2$altern[i]==1){mydata_tr2$altern[i]="Option1"}
  else if (mydata_tr2$altern[i]==2){mydata_tr2$altern[i]="Option2"}
  else if (mydata_tr2$altern[i]==3){mydata_tr2$altern[i]="Option3"}}
head(mydata_tr2)

##labeled the alternative with reduced sugar content as 1
mydata_tr2$Reduced=c()
for (i in 1:4500){
  if (mydata_tr2$TSugar[i]==14){mydata_tr2$Reduced[i]=1}
  else if (mydata_tr2$TSugar[i]!=14){mydata_tr2$Reduced[i]=0}
}
head(mydata_tr2)

##labeled the alternative with regular sugar content as 1
mydata_tr2$Regular=c()
for (i in 1:4500){
  if (mydata_tr2$TSugar[i]==29){mydata_tr2$Regular[i]=1}
  else if (mydata_tr2$TSugar[i]!=29){mydata_tr2$Regular[i]=0}
}
head(mydata_tr2)

##generate bland flavor variable
mydata_tr2$BlandFlavor=1-mydata_tr2$FullFlav
for (i in 1:1500){
  mydata_tr2$BlandFlavor[3*i]=0
}
head(mydata_tr2)

mean(mydata_tr2$FullFlav)
mean(mydata_tr2$BlandFlavor)

##############################################################################
#####################   PREFERENCE SPACE MODELS    ###########################
##############################################################################

TM_tr2 <- mlogit.data(mydata_tr2,choice="Choice",shape="long",alt.var="altern",
                      chid.var="id",id.var="PersonID")
head(TM_tr2)

##MNL
mnl_tr2 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,data=TM_tr2,model="mnl")
summary(mnl_tr2)
getSummary.gmnl(mnl_tr2,alpha=0.05)
-1*wtp.gmnl(mnl_tr2,wrt="price")

#####################################################################
#######################  RPL or MIXED LOGIT  ########################
#####################################################################

##Bopt is random, corr
rplc_or_tr2 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                    data=TM_tr2,
                    ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                    R=500, 
                    panel=TRUE,
                    correlation=TRUE,
                    model="mixl")
summary(rplc_or_tr2)
-1*wtp.gmnl(rplc_or_tr2,wrt="price")
cor.gmnl(rplc_or_tr2)
getSummary.gmnl(rplc_or_tr2,alpha=0.05)

##############################################################################
#############################    WTP-SPACE MODELS    #########################
##############################################################################

wtp_data_tr2<-mlogit.data(mydata_tr2,choice="Choice",shape="long",alt.var="altern",
                          chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data_tr2)

#######################################################################
#########################  G-MNL-II  ##################################
#######################################################################

####Bopt is random, corr
wtp_rplc_or_tr2<- -1*wtp.gmnl(rplc_or_tr2,wrt="price")
wtp_rplc_or_tr2<-data.frame(wtp_rplc_or_tr2)

start_or_tr2_1 <-c(1,wtp_rplc_or_tr2$Estimate[1:4],rep(0.1,11),0.2,0)
start_or_tr2_2 <- c(1,wtp_rplc_or_tr2$Estimate[1:4],rep(0.1,11),0,0)
start_or_tr2_3 <- c(1,rep(0.1,15),0.1,0)
start_or_tr2_4 <- c(1,rep(1,15),1,0)
start_or_tr2_5 <- c(1,rep(-1,15),-1,0)

start=rbind(start_or_tr2_1,start_or_tr2_2,start_or_tr2_3,
            start_or_tr2_4,start_or_tr2_5)
wtpc_or_tr2 <-list()
ll_tr2 <- list()
for (i in 1:5){
  wtpc_or_tr2[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt
                           |0|0|0|certain-1,
                           data=wtp_data_tr2,
                           ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                           R=500,
                           #halton = NA,
                           print.init = TRUE, 
                           panel=TRUE,
                           start=start[i,],
                           fixed=c(TRUE,rep(FALSE,16),TRUE), 
                           model="gmnl",
                           correlation=TRUE,
                           notscale=c(0,0,0,0,1),
                           method="bhhh",
                           iterlim=500)
  ll_tr2[[i]]=getSummary.gmnl(wtpc_or_tr2[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_tr2

summary(wtpc_or_tr2[[3]])
cor.gmnl(wtpc_or_tr2[[3]])
getSummary.gmnl(wtpc_or_tr2[[3]],alpha=0.05)

#### conditional WTPs
indwtp_tr2<-effect.gmnl(wtpc_or_tr2[[1]],effect="ce",
                        par=c("Regular","Gene","BlandFlavor"))
#indwtp_tr2


#---------------------------------------------------------------------------
#-----------------------------  TREATMENT 3  -------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")
mydata_tr3<-read.csv('Dried_Total_3.csv', fileEncoding = 'UTF-8-BOM')
head(mydata_tr3)

mean(mydata_tr3$Scale)
sum(mydata_tr3$Scale <= 7, na.rm=TRUE)
sum(mydata_tr3$Scale < 7, na.rm=TRUE)
sum(mydata_tr3$Scale >= 7, na.rm=TRUE)/4500

mydata_tr3$cset <- NULL
mydata_tr3$Block <- NULL
mydata_tr3$Option1 <- NULL
mydata_tr3$Option2 <- NULL
mydata_tr3$Option3 <- NULL

##generate certainty dummy
mydata_tr3$certain=c()
for (i in 1:4500){
  if (mydata_tr3$Scale[i]>=7){mydata_tr3$certain[i]=1}
  else if (mydata_tr3$Scale[i]<7){mydata_tr3$certain[i]=0}
} 
##check uncertainty dummy
mydata_tr3$uncertain=c()
for (i in 1:4500){
  if (mydata_tr3$Scale[i]<7){mydata_tr3$uncertain[i]=1}
  else if (mydata_tr3$Scale[i]>=7){mydata_tr3$uncertain[i]=0}
} 
##check 
mean(mydata_tr3$certain)
sum(mydata_tr3$Scale>=7)/4500

mydata_tr3$id=c()
for (i in 1:4500){
  mydata_tr3$id[i]=ceiling(i/3)
}

for (i in 1:4500){
  if (mydata_tr3$altern[i]==1){mydata_tr3$altern[i]="Option1"}
  else if (mydata_tr3$altern[i]==2){mydata_tr3$altern[i]="Option2"}
  else if (mydata_tr3$altern[i]==3){mydata_tr3$altern[i]="Option3"}}
head(mydata_tr3)


##labeled the alternative with reduced sugar content as 1
mydata_tr3$Reduced=c()
for (i in 1:4500){
  if (mydata_tr3$TSugar[i]==14){mydata_tr3$Reduced[i]=1}
  else if (mydata_tr3$TSugar[i]!=14){mydata_tr3$Reduced[i]=0}
}
head(mydata_tr3)

##labeled the alternative with regular sugar content as 1
mydata_tr3$Regular=c()
for (i in 1:4500){
  if (mydata_tr3$TSugar[i]==29){mydata_tr3$Regular[i]=1}
  else if (mydata_tr3$TSugar[i]!=29){mydata_tr3$Regular[i]=0}
}
head(mydata_tr3)

##generate bland flavor variable
mydata_tr3$BlandFlavor=1-mydata_tr3$FullFlav
for (i in 1:1500){
  mydata_tr3$BlandFlavor[3*i]=0
}
head(mydata_tr3)

mean(mydata_tr3$FullFlav)
mean(mydata_tr3$BlandFlavor)

##############################################################################
#####################   PREFERENCE SPACE MODELS    ###########################
##############################################################################

TM_tr3 <- mlogit.data(mydata_tr3,choice="Choice",shape="long",alt.var="altern",
                      chid.var="id",id.var="PersonID")
head(TM_tr3)

##MNL
mnl_tr3 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,data=TM_tr3,model="mnl")
summary(mnl_tr3)
getSummary.gmnl(mnl_tr3,alpha=0.05)
-1*wtp.gmnl(mnl_tr3,wrt="price")

#####################################################################
#######################  RPL or MIXED LOGIT  ########################
#####################################################################

##Bopt is random, corr
rplc_or_tr3 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                    data=TM_tr3,
                    ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                    R=500, 
                    panel=TRUE,
                    correlation=TRUE,
                    model="mixl")
summary(rplc_or_tr3)
-1*wtp.gmnl(rplc_or_tr3,wrt="price")
cor.gmnl(rplc_or_tr3)
getSummary.gmnl(rplc_or_tr3,alpha=0.05)

##############################################################################
#############################    WTP-SPACE MODELS    #########################
##############################################################################

wtp_data_tr3<-mlogit.data(mydata_tr3,choice="Choice",shape="long",alt.var="altern",
                          chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data_tr3)

#######################################################################
#########################  G-MNL-II  ##################################
#######################################################################

####Bopt is random, corr
wtp_rplc_or_tr3<- -1*wtp.gmnl(rplc_or_tr3,wrt="price")
wtp_rplc_or_tr3<-data.frame(wtp_rplc_or_tr3)
start_or_tr3_1 <-c(1,wtp_rplc_or_tr3$Estimate[1:4],rep(0.1,11),0.2,0)
start_or_tr3_2 <- c(1,wtp_rplc_or_tr3$Estimate[1:4],rep(0.1,11),0,0)
start_or_tr3_3 <- c(1,rep(0.1,15),0.1,0)
start_or_tr3_4 <- c(1,rep(1,15),1,0)
start_or_tr3_5 <- c(1,rep(-1,15),-1,0)

start_tr3=rbind(start_or_tr3_1,start_or_tr3_2,start_or_tr3_3,start_or_tr3_4,start_or_tr3_5)



wtpc_or_tr3 <-list()
ll_tr3 <- list()

for (i in 1:5){
  wtpc_or_tr3[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt
                           |0|0|0|certain-1,
                           data=wtp_data_tr3,
                           ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                           R=500,
                           print.init = TRUE, 
                           panel=TRUE,
                           start=start_tr3[i,],
                           fixed=c(TRUE,rep(FALSE,16),TRUE), 
                           model="gmnl",
                           correlation=TRUE,
                           notscale=c(0,0,0,0,1),
                           method="bhhh",
                           iterlim=500)
  ll_tr3[[i]]=getSummary.gmnl(wtpc_or_tr3[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_tr3

summary(wtpc_or_tr3[[5]])
cor.gmnl(wtpc_or_tr3[[5]])
getSummary.gmnl(wtpc_or_tr3[[5]],alpha=0.05)

#### conditional WTPs
indwtp_tr3<-effect.gmnl(wtpc_or_tr3[[1]],effect="ce",
                        par=c("Regular","Gene","BlandFlavor"))
#indwtp_tr3


#---------------------------------------------------------------------------
#-----------------------------  TREATMENT 4  -------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")
mydata_tr4<-read.csv('Dried_Total_4.csv', fileEncoding = 'UTF-8-BOM')
head(mydata_tr4)

mean(mydata_tr4$Scale)
sum(mydata_tr4$Scale <= 7, na.rm=TRUE)
sum(mydata_tr4$Scale < 7, na.rm=TRUE)

mydata_tr4$cset <- NULL
mydata_tr4$Block <- NULL
mydata_tr4$Option1 <- NULL
mydata_tr4$Option2 <- NULL
mydata_tr4$Option3 <- NULL

##generate certainty dummy
mydata_tr4$certain=c()
for (i in 1:4500){
  if (mydata_tr4$Scale[i]>=7){mydata_tr4$certain[i]=1}
  else if (mydata_tr4$Scale[i]<7){mydata_tr4$certain[i]=0}
} 
##check uncertainty dummy
mydata_tr4$uncertain=c()
for (i in 1:4500){
  if (mydata_tr4$Scale[i]<7){mydata_tr4$uncertain[i]=1}
  else if (mydata_tr4$Scale[i]>=7){mydata_tr4$uncertain[i]=0}
} 
##check 
mean(mydata_tr4$certain)
sum(mydata_tr4$Scale>=7)/4500

mydata_tr4$id=c()
for (i in 1:4500){
  mydata_tr4$id[i]=ceiling(i/3)
}

for (i in 1:4500){
  if (mydata_tr4$altern[i]==1){mydata_tr4$altern[i]="Option1"}
  else if (mydata_tr4$altern[i]==2){mydata_tr4$altern[i]="Option2"}
  else if (mydata_tr4$altern[i]==3){mydata_tr4$altern[i]="Option3"}}
head(mydata_tr4)

##labeled the alternative with reduced sugar content as 1
mydata_tr4$Reduced=c()
for (i in 1:4500){
  if (mydata_tr4$TSugar[i]==14){mydata_tr4$Reduced[i]=1}
  else if (mydata_tr4$TSugar[i]!=14){mydata_tr4$Reduced[i]=0}
}
head(mydata_tr4)


##labeled the alternative with regular sugar content as 1
mydata_tr4$Regular=c()
for (i in 1:4500){
  if (mydata_tr4$TSugar[i]==29){mydata_tr4$Regular[i]=1}
  else if (mydata_tr4$TSugar[i]!=29){mydata_tr4$Regular[i]=0}
}
head(mydata_tr4)

##generate bland flavor variable
mydata_tr4$BlandFlavor=1-mydata_tr4$FullFlav
for (i in 1:1500){
  mydata_tr4$BlandFlavor[3*i]=0
}
head(mydata_tr4)

mean(mydata_tr4$FullFlav)
mean(mydata_tr4$BlandFlavor)

##############################################################################
#####################   PREFERENCE SPACE MODELS    ###########################
##############################################################################

TM_tr4 <- mlogit.data(mydata_tr4,choice="Choice",shape="long",alt.var="altern",
                      chid.var="id",id.var="PersonID")
head(TM_tr4)

##MNL
mnl_tr4 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,data=TM_tr4,model="mnl")
summary(mnl_tr4)
getSummary.gmnl(mnl_tr4,alpha=0.05)
-1*wtp.gmnl(mnl_tr4,wrt="price")

#####################################################################
#######################  RPL or MIXED LOGIT  ########################
#####################################################################

##Bopt is random, corr
rplc_or_tr4 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                    data=TM_tr4,
                    ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                    R=500, 
                    #halton=NA,
                    panel=TRUE,
                    correlation=TRUE,
                    method="bhhh",
                    model="mixl")
summary(rplc_or_tr4)
-1*wtp.gmnl(rplc_or_tr4,wrt="price")
cor.gmnl(rplc_or_tr4)
getSummary.gmnl(rplc_or_tr4,alpha=0.05)

##############################################################################
#############################    WTP-SPACE MODELS    #########################
##############################################################################

wtp_data_tr4<-mlogit.data(mydata_tr4,choice="Choice",shape="long",alt.var="altern",
                          chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data_tr4)

#######################################################################
#########################  G-MNL-II  ##################################
#######################################################################

####Bopt is random, corr
wtp_rplc_or_tr4<- -1*wtp.gmnl(rplc_or_tr4,wrt="price")
wtp_rplc_or_tr4<-data.frame(wtp_rplc_or_tr4)

start_or_tr4_1 <-c(1,wtp_rplc_or_tr4$Estimate[1:4],rep(0.1,11),0.2,0)
start_or_tr4_2 <- c(1,wtp_rplc_or_tr4$Estimate[1:4],rep(0.1,11),0,0)
start_or_tr4_3 <- c(1,rep(0.1,15),0.1,0)
start_or_tr4_4 <- c(1,rep(1,15),1,0)
start_or_tr4_5 <- c(1,rep(-1,15),-1,0)

start_tr4=rbind(start_or_tr4_1,start_or_tr4_2,start_or_tr4_3,
                start_or_tr4_4,start_or_tr4_5)

wtpc_or_tr4 <-list()
ll_tr4 <- list()

for (i in 1:5){
  wtpc_or_tr4[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt
                           |0|0|0|certain-1,
                           data=wtp_data_tr4,
                           ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                           R=1000,
                           #halton = NA,
                           print.init = TRUE, 
                           panel=TRUE,
                           start=start_tr4[i,],
                           fixed=c(TRUE,rep(FALSE,16),TRUE), 
                           model="gmnl",
                           correlation=TRUE,
                           notscale=c(0,0,0,0,1),
                           method="bhhh",
                           iterlim=500)
  ll_tr4[[i]]=getSummary.gmnl(wtpc_or_tr4[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_tr4

summary(wtpc_or_tr4[[5]])
cor.gmnl(wtpc_or_tr4[[5]])
getSummary.gmnl(wtpc_or_tr4[[5]],alpha=0.05)

#### conditional WTPs
indwtp_tr4<-effect.gmnl(wtpc_or_tr4[[1]],effect="ce",
                        par=c("Regular","Gene","BlandFlavor"))
#indwtp_tr4


set.seed(1118)
### Unconditional WTP's ###################################
##total sugar
un_indwtp_Regular_tr1 <- wtpc_or_tr1[[2]]$coefficients["Regular"] + (abs(wtpc_or_tr1[[2]]$coefficients["sd.Regular.Regular"]))*rnorm(250)
mean(un_indwtp_Regular_tr1)
un_indwtp_Regular_tr2 <- wtpc_or_tr2[[3]]$coefficients["Regular"] + (abs(wtpc_or_tr2[[3]]$coefficients["sd.Regular.Regular"]))*rnorm(250)
mean(un_indwtp_Regular_tr2)
un_indwtp_Regular_tr3 <- wtpc_or_tr3[[5]]$coefficients["Regular"] + (abs(wtpc_or_tr3[[5]]$coefficients["sd.Regular.Regular"]))*rnorm(250)
mean(un_indwtp_Regular_tr3)
un_indwtp_Regular_tr4 <- wtpc_or_tr4[[5]]$coefficients["Regular"] + (abs(wtpc_or_tr4[[5]]$coefficients["sd.Regular.Regular"]))*rnorm(250)
mean(un_indwtp_Regular_tr4)
##gene
un_indwtp_gene_tr1 <- wtpc_or_tr1[[2]]$coefficients["Gene"] + (abs(wtpc_or_tr1[[2]]$coefficients["sd.Gene.Gene"]))*rnorm(250)
mean(un_indwtp_gene_tr1)
un_indwtp_gene_tr2 <- wtpc_or_tr2[[3]]$coefficients["Gene"] + (abs(wtpc_or_tr2[[3]]$coefficients["sd.Gene.Gene"]))*rnorm(250)
mean(un_indwtp_gene_tr2)
un_indwtp_gene_tr3 <- wtpc_or_tr3[[5]]$coefficients["Gene"] + (abs(wtpc_or_tr3[[5]]$coefficients["sd.Gene.Gene"]))*rnorm(250)
mean(un_indwtp_gene_tr3)
un_indwtp_gene_tr4 <- wtpc_or_tr4[[5]]$coefficients["Gene"] + (abs(wtpc_or_tr4[[5]]$coefficients["sd.Gene.Gene"]))*rnorm(250)
mean(un_indwtp_gene_tr4)
##fullflav
un_indwtp_blandflavor_tr1 <- wtpc_or_tr1[[2]]$coefficients["BlandFlavor"] + (abs(wtpc_or_tr1[[2]]$coefficients["sd.BlandFlavor.BlandFlavor"]))*rnorm(250)
mean(un_indwtp_blandflavor_tr1)
un_indwtp_blandflavor_tr2 <- wtpc_or_tr2[[3]]$coefficients["BlandFlavor"] + (abs(wtpc_or_tr2[[3]]$coefficients["sd.BlandFlavor.BlandFlavor"]))*rnorm(250)
mean(un_indwtp_blandflavor_tr2)
un_indwtp_blandflavor_tr3 <- wtpc_or_tr3[[5]]$coefficients["BlandFlavor"] + (abs(wtpc_or_tr3[[5]]$coefficients["sd.BlandFlavor.BlandFlavor"]))*rnorm(250)
mean(un_indwtp_blandflavor_tr3)
un_indwtp_blandflavor_tr4 <- wtpc_or_tr4[[5]]$coefficients["BlandFlavor"] + (abs(wtpc_or_tr4[[5]]$coefficients["sd.BlandFlavor.BlandFlavor"]))*rnorm(250)
mean(un_indwtp_blandflavor_tr4)


###################t test
##regular sugar content
t.test(un_indwtp_Regular_tr1,un_indwtp_Regular_tr2,alternative="less")
t.test(un_indwtp_Regular_tr1,un_indwtp_Regular_tr3,alternative="greater")
t.test(un_indwtp_Regular_tr1,un_indwtp_Regular_tr4,alternative="two.sided")

##gene
t.test(un_indwtp_gene_tr1,un_indwtp_gene_tr2,alternative="less")
t.test(un_indwtp_gene_tr1,un_indwtp_gene_tr3,alternative="greater")
t.test(un_indwtp_gene_tr1,un_indwtp_gene_tr4,alternative="two.sided")

##bland flavor
t.test(un_indwtp_blandflavor_tr1,un_indwtp_blandflavor_tr2,alternative="less")
t.test(un_indwtp_blandflavor_tr1,un_indwtp_blandflavor_tr3,alternative="greater")
t.test(un_indwtp_blandflavor_tr1,un_indwtp_blandflavor_tr4,alternative="two.sided")








################################################################################
#------------------------------------------------------------------------------#
#-------------------------------  pooled data  --------------------------------#
#------------------------------------------------------------------------------#
################################################################################

#---------------------------------------------------------------------------
#---------------------------------  POOLED ---------------------------------
#---------------------------------------------------------------------------
library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

##read data
setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")
mydata1<-read.csv('Dried_Total_1.csv', fileEncoding = 'UTF-8-BOM')
mydata2<-read.csv('Dried_Total_2.csv', fileEncoding = 'UTF-8-BOM')
mydata3<-read.csv('Dried_Total_3.csv', fileEncoding = 'UTF-8-BOM')
mydata4<-read.csv('Dried_Total_4.csv', fileEncoding = 'UTF-8-BOM')

##combined 4 data 
data_pool <- rbind(mydata1,mydata2,mydata3,mydata4)
dim(data_pool)

data_pool$cset <- NULL
data_pool$Block <- NULL
data_pool$Option1 <- NULL
data_pool$Option2 <- NULL
data_pool$Option3 <- NULL

##generate certainty dummy
data_pool$certain=c()
for (i in 1:18000){
  if (data_pool$Scale[i]>=7){data_pool$certain[i]=1}
  else if (data_pool$Scale[i]<7){data_pool$certain[i]=0}
} 
##check uncertainty dummy
data_pool$uncertain=c()
for (i in 1:18000){
  if (data_pool$Scale[i]<7){data_pool$uncertain[i]=1}
  else if (data_pool$Scale[i]>=7){data_pool$uncertain[i]=0}
} 
##check 
mean(data_pool$certain)
sum(data_pool$Scale>=7)/18000

data_pool$id=c()
for (i in 1:18000){
  data_pool$id[i]=ceiling(i/3)
}
##change name of alternatives
for (i in 1:18000){
  if (data_pool$altern[i]==1){data_pool$altern[i]="Option1"}
  else if (data_pool$altern[i]==2){data_pool$altern[i]="Option2"}
  else if (data_pool$altern[i]==3){data_pool$altern[i]="Option3"}}
##view the cleaned data
head(data_pool)
tail(data_pool) 


##labeled the alternative with reduced sugar content as 1
data_pool$Reduced=c()
for (i in 1:18000){
  if (data_pool$TSugar[i]==14){data_pool$Reduced[i]=1}
  else if (data_pool$TSugar[i]!=14){data_pool$Reduced[i]=0}
}
head(data_pool)

##############################################################################
#####################   PREFERENCE SPACE MODELS    ###########################
##############################################################################

TM_pool <- mlogit.data(data_pool,choice="Choice",shape="long",alt.var="altern",
                       chid.var="id",id.var="PersonID")
head(TM_pool)

##MNL
mnl_pool <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,data=TM_pool,model="mnl")
summary(mnl_pool)
getSummary.gmnl(mnl_pool,alpha=0.05)
-1*wtp.gmnl(mnl_pool,wrt="price")

#####################################################################
#######################  RPL or MIXED LOGIT  ########################
#####################################################################

##Bopt is random, corr
rplc_or_pool <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                     data=TM_pool,
                     ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                     R=500, 
                     panel=TRUE,
                     correlation=TRUE,
                     model="mixl")
summary(rplc_or_pool)
-1*wtp.gmnl(rplc_or_pool,wrt="price")
cor.gmnl(rplc_or_pool)
getSummary.gmnl(rplc_or_pool,alpha=0.05)

##############################################################################
#############################    WTP-SPACE MODELS    #########################
##############################################################################

wtp_data_pool<-mlogit.data(data_pool,choice="Choice",shape="long",alt.var="altern",
                           chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data_pool)

#######################################################################
#########################  G-MNL-II  ##################################
#######################################################################

####Bopt is random, corr
wtp_rplc_or_pool<- -1*wtp.gmnl(rplc_or_pool,wrt="price")
wtp_rplc_or_pool<-data.frame(wtp_rplc_or_pool)

start_or_pool_1 <-c(1,wtp_rplc_or_pool$Estimate[1:4],rep(0.1,11),0.2,0)
start_or_pool_2 <- c(1,wtp_rplc_or_pool$Estimate[1:4],rep(0.1,11),0,0)
start_or_pool_3 <- c(1,rep(0.1,15),0.1,0)
start_or_pool_4 <- c(1,rep(1,15),1,0)
start_or_pool_5 <- c(1,rep(-1,15),-1,0)

start=rbind(start_or_pool_1,start_or_pool_2,start_or_pool_3,
            start_or_pool_4,start_or_pool_5)
wtpc_or_pool <-list()
ll_pool <- list()
for (i in 1:5){
  wtpc_or_pool[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|0|0|0|certain-1,
                            data=wtp_data_pool,
                            ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                            R=500,
                            #halton = NA,
                            print.init = TRUE, 
                            panel=TRUE,
                            start=start[i,],
                            fixed=c(TRUE,rep(FALSE,16),TRUE), 
                            model="gmnl",
                            correlation=TRUE,
                            notscale=c(0,0,0,0,1),
                            method="bhhh",
                            iterlim=500)
  ll_pool[[i]]=getSummary.gmnl(wtpc_or_pool[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_pool

summary(wtpc_or_pool[[5]])
cor.gmnl(wtpc_or_pool[[5]])
getSummary.gmnl(wtpc_or_pool[[5]],alpha=0.05)

######################################################################
#############################  GMNL-I  ###############################
######################################################################

####Bopt is random corr
start_or1_pool <- c(1,wtp_rplc_or_pool$Estimate[1:4],rep(0.1,11),0.2,1)  
wtpc_or1_pool <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|0|0|0|certain-1,
                      data=wtp_data_pool,
                      ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                      R=500,
                      print.init = TRUE, 
                      panel=TRUE,
                      start=start_or1_pool,
                      fixed=c(TRUE,rep(FALSE,16),TRUE), 
                      model="gmnl",
                      correlation=TRUE,
                      notscale=c(0,0,0,0,1),
                      method="bhhh",
                      iterlim=500)
summary(wtpc_or1_pool)
cor.gmnl(wtpc_or1_pool)
getSummary.gmnl(wtpc_or1_pool)

###############################G-MNL
###################################################
gmnl_pool <-list()
gmnl_ll_pool<- list()
for (i in 1:5){
  gmnl_pool[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt
                         |0|0|0|certain-1,
                         data=wtp_data_pool,
                         ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                         R=500,
                         print.init = TRUE, 
                         panel=TRUE,
                         start=start[i,],
                         fixed=c(TRUE,rep(FALSE,16),FALSE), 
                         model="gmnl",
                         correlation=TRUE,
                         notscale=c(0,0,0,0,1),
                         method="bhhh",
                         iterlim=500)
  gmnl_ll_pool[[i]]=getSummary.gmnl(gmnl_pool[[i]],alpha=0.05)$sumstat["logLik"]
}
gmnl_ll_pool
summary(gmnl_pool[[1]])
cor.gmnl(gmnl_pool[[1]])
getSummary.gmnl(gmnl_pool[[1]],alpha=0.05)

#---------------------------------------------------------------------------
#-----------------------------  POOLED 1 & 2  ------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")

##read treatment 1 data
mydata1<-read.csv('Dried_Total_1.csv', fileEncoding = 'UTF-8-BOM')
head(mydata1)
dim(mydata1)
##read treatment 2 data
mydata2<-read.csv('Dried_Total_2.csv', fileEncoding = 'UTF-8-BOM')
head(mydata2)
dim(mydata2)

##combine data1 and data2
data12 <- rbind(mydata1,mydata2)
dim(data12)

##data clean, delete the irrelevant columns
data12$cset <- NULL
data12$Block <- NULL
data12$Option1 <- NULL
data12$Option2 <- NULL
data12$Option3 <- NULL

##generate certainty dummy
data12$certain=c()
for (i in 1:9000){
  if (data12$Scale[i]>=7){data12$certain[i]=1}
  else if (data12$Scale[i]<7){data12$certain[i]=0}
} 
##check uncertainty dummy
data12$uncertain=c()
for (i in 1:9000){
  if (data12$Scale[i]<7){data12$uncertain[i]=1}
  else if (data12$Scale[i]>=7){data12$uncertain[i]=0}
} 
##check 
mean(data12$certain)
sum(data12$Scale>=7)/9000

data12$id=c()
for (i in 1:9000){
  data12$id[i]=ceiling(i/3)
}
##change name of alternatives
for (i in 1:9000){
  if (data12$altern[i]==1){data12$altern[i]="Option1"}
  else if (data12$altern[i]==2){data12$altern[i]="Option2"}
  else if (data12$altern[i]==3){data12$altern[i]="Option3"}}
##view the cleaned data
head(data12)
tail(data12)  ##number of obs: 3000, there are 3000 group id

TM_12 <- mlogit.data(data12,choice="Choice",shape="long",alt.var="altern",
                     chid.var="id",id.var="PersonID")
head(TM_12)

rplc_or_12 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                   data=TM_12,
                   ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                   R=500, 
                   panel=TRUE,
                   correlation=TRUE,
                   model="mixl")
summary(rplc_or_12)
-1*wtp.gmnl(rplc_or_12,wrt="price")
cor.gmnl(rplc_or_12)

wtp_data12<-mlogit.data(data12,choice="Choice",shape="long",alt.var="altern",
                        chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data12)

####Bopt is random
wtp_rplc_or_12<- -1*wtp.gmnl(rplc_or_12,wrt="price")
wtp_rplc_or_12<-data.frame(wtp_rplc_or_12)

start_or_12_1 <- c(1,wtp_rplc_or_12$Estimate[1:4],rep(0.1,11),0.2,0) 
start_or_12_2 <- c(1,rep(0.1,15),0.1,0) 
start_or_12_3 <- c(1,wtp_rplc_or_12$Estimate[1:4],rep(0.1,11),0,0) 
start_or_12_4 <- c(1,rep(1,15),1,0)
start_or_12_5 <- c(1,rep(-1,15),-1,0)

start_12=rbind(start_or_12_1,start_or_12_2,start_or_12_3,
               start_or_12_4,start_or_12_5)

wtpc_or_12 <-list()
ll_12 <- list()

for (i in 1:5){
  wtpc_or_12[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|0|0|0|certain-1,
                          data=wtp_data12,
                          ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                          R=500,
                          #halton = NA,
                          print.init = TRUE, 
                          panel=TRUE,
                          start=start_12[i,],
                          fixed=c(TRUE,rep(FALSE,16),TRUE), 
                          model="gmnl",
                          correlation=TRUE,
                          notscale=c(0,0,0,0,1),
                          method="bhhh",
                          iterlim=500)
  ll_12[[i]]=getSummary.gmnl(wtpc_or_12[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_12

summary(wtpc_or_12[[1]])
cor.gmnl(wtpc_or_12[[1]])
getSummary.gmnl(wtpc_or_12[[1]],alpha=0.05)


#---------------------------------------------------------------------------
#-----------------------------  POOLED 1 & 3  ------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")

##read treatment 1 data
mydata1<-read.csv('Dried_Total_1.csv', fileEncoding = 'UTF-8-BOM')
head(mydata1)
dim(mydata1)
##read treatment 3 data
mydata3<-read.csv('Dried_Total_3.csv', fileEncoding = 'UTF-8-BOM')
head(mydata3)
dim(mydata3)

##combine data1 and data3
data13 <- rbind(mydata1,mydata3)
dim(data13)

data13$cset <- NULL
data13$Block <- NULL
data13$Option1 <- NULL
data13$Option2 <- NULL
data13$Option3 <- NULL

##generate certainty dummy
data13$certain=c()
for (i in 1:9000){
  if (data13$Scale[i]>=7){data13$certain[i]=1}
  else if (data13$Scale[i]<7){data13$certain[i]=0}
} 
##check uncertainty dummy
data13$uncertain=c()
for (i in 1:9000){
  if (data13$Scale[i]<7){data13$uncertain[i]=1}
  else if (data13$Scale[i]>=7){data13$uncertain[i]=0}
} 
##check 
mean(data13$certain)
sum(data13$Scale>=7)/9000


data13$id=c()
for (i in 1:9000){
  data13$id[i]=ceiling(i/3)
}
for (i in 1:9000){
  if (data13$altern[i]==1){data13$altern[i]="Option1"}
  else if (data13$altern[i]==2){data13$altern[i]="Option2"}
  else if (data13$altern[i]==3){data13$altern[i]="Option3"}}
head(data13)
tail(data13)  ##number of obs: 3000, there are 3000 group id

TM_13 <- mlogit.data(data13,choice="Choice",shape="long",alt.var="altern",
                     chid.var="id",id.var="PersonID")
head(TM_13)

rplc_or_13 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                   data=TM_13,
                   ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                   R=500, 
                   #halton=NA,
                   panel=TRUE,
                   correlation=TRUE,
                   notscale=c(0,0,0,0,1),
                   model="mixl")
summary(rplc_or_13)
-1*wtp.gmnl(rplc_or_13,wrt="price")
cor.gmnl(rplc_or_13)

wtp_data13<-mlogit.data(data13,choice="Choice",shape="long",alt.var="altern",
                        chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data13)

####Bopt is random
wtp_rplc_or_13<- -1*wtp.gmnl(rplc_or_13,wrt="price")
wtp_rplc_or_13<-data.frame(wtp_rplc_or_13)

start_or_13_1 <- c(1,wtp_rplc_or_13$Estimate[1:4],rep(0.1,11),0.2,0) 
start_or_13_2 <- c(1,rep(0.1,15),0.1,0) 
start_or_13_3 <- c(1,wtp_rplc_or_13$Estimate[1:4],rep(0.1,11),0,0) 
start_or_13_4 <- c(1,rep(1,15),1,0)
start_or_13_5 <- c(1,rep(-1,15),-1,0)

start_13=rbind(start_or_13_1,start_or_13_2,start_or_13_3,
               start_or_13_4,start_or_13_5)

wtpc_or_13 <-list()
ll_13 <- list()

for (i in 1:5){
  wtpc_or_13[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|0|0|0|certain-1,
                          data=wtp_data13,
                          ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                          R=500,
                          #halton = NA,
                          print.init = TRUE, 
                          panel=TRUE,
                          start=start_13[i,],
                          fixed=c(TRUE,rep(FALSE,16),TRUE), 
                          model="gmnl",
                          correlation=TRUE,
                          notscale=c(0,0,0,0,1),
                          method="bhhh",
                          iterlim=500)
  ll_13[[i]]=getSummary.gmnl(wtpc_or_13[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_13

summary(wtpc_or_13[[1]])
cor.gmnl(wtpc_or_13[[1]])
getSummary.gmnl(wtpc_or_13[[1]],alpha=0.05)


#---------------------------------------------------------------------------
#-----------------------------  POOLED 1 & 4  ------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")
##read treatment 1 data
mydata1<-read.csv('Dried_Total_1.csv', fileEncoding = 'UTF-8-BOM')
head(mydata1)
dim(mydata1)
##read treatment 4 data
mydata4<-read.csv('Dried_Total_4.csv', fileEncoding = 'UTF-8-BOM')
head(mydata4)
dim(mydata4)

data14 <- rbind(mydata1,mydata4)
dim(data14)

data14$cset <- NULL
data14$Block <- NULL
data14$Option1 <- NULL
data14$Option2 <- NULL
data14$Option3 <- NULL

##generate certainty dummy
data14$certain=c()
for (i in 1:9000){
  if (data14$Scale[i]>=7){data14$certain[i]=1}
  else if (data14$Scale[i]<7){data14$certain[i]=0}
} 
##check uncertainty dummy
data14$uncertain=c()
for (i in 1:9000){
  if (data14$Scale[i]<7){data14$uncertain[i]=1}
  else if (data14$Scale[i]>=7){data14$uncertain[i]=0}
} 
##check 
mean(data14$certain)
sum(data14$Scale>=7)/9000

data14$id=c()
for (i in 1:9000){
  data14$id[i]=ceiling(i/3)
}
for (i in 1:9000){
  if (data14$altern[i]==1){data14$altern[i]="Option1"}
  else if (data14$altern[i]==2){data14$altern[i]="Option2"}
  else if (data14$altern[i]==3){data14$altern[i]="Option3"}}
head(data14)
tail(data14)  ##number of obs: 3000, there are 3000 group id

TM_14 <- mlogit.data(data14,choice="Choice",shape="long",alt.var="altern",
                     chid.var="id",id.var="PersonID")
head(TM_14)

rplc_or_14 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                   data=TM_14,
                   ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                   R=500, 
                   #halton=NA,
                   panel=TRUE,
                   correlation=TRUE,
                   model="mixl")
summary(rplc_or_14)
-1*wtp.gmnl(rplc_or_14,wrt="price")
cor.gmnl(rplc_or_14)

wtp_data14<-mlogit.data(data14,choice="Choice",shape="long",alt.var="altern",
                        chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data14)

####Bopt is random
wtp_rplc_or_14<- -1*wtp.gmnl(rplc_or_14,wrt="price")
wtp_rplc_or_14<-data.frame(wtp_rplc_or_14)

start_or_14_1 <- c(1,wtp_rplc_or_14$Estimate[1:4],rep(0.1,11),0.2,0) 
start_or_14_2 <- c(1,rep(0.1,15),0.1,0) 
start_or_14_3 <- c(1,wtp_rplc_or_14$Estimate[1:4],rep(0.1,11),0,0) 
start_or_14_4 <- c(1,rep(1,15),1,0)
start_or_14_5 <- c(1,rep(-1,15),-1,0)

start_14=rbind(start_or_14_1,start_or_14_2,start_or_14_3,
               start_or_14_4,start_or_14_5)

wtpc_or_14 <-list()
ll_14 <- list()

for (i in 1:5){
  wtpc_or_14[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|0|0|0|certain-1,
                          data=wtp_data14,
                          ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                          R=500,
                          #halton = NA,
                          print.init = TRUE, 
                          panel=TRUE,
                          start=start_14[i,],
                          fixed=c(TRUE,rep(FALSE,16),TRUE), 
                          model="gmnl",
                          correlation=TRUE,                           
                          notscale=c(0,0,0,0,1),
                          method="bhhh",
                          iterlim=500)
  ll_14[[i]]=getSummary.gmnl(wtpc_or_14[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_14

summary(wtpc_or_14[[1]])
cor.gmnl(wtpc_or_14[[1]])
getSummary.gmnl(wtpc_or_14[[1]],alpha=0.05)


#---------------------------------------------------------------------------
#-----------------------------  POOLED 2 & 3  ------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")

##read treatment 2 data
mydata2<-read.csv('Dried_Total_2.csv', fileEncoding = 'UTF-8-BOM')
head(mydata2)
dim(mydata2)
##read treatment 3 data
mydata3<-read.csv('Dried_Total_3.csv', fileEncoding = 'UTF-8-BOM')
head(mydata3)
dim(mydata3)

data23 <- rbind(mydata2,mydata3)
dim(data23)

data23$cset <- NULL
data23$Block <- NULL
data23$Option1 <- NULL
data23$Option2 <- NULL
data23$Option3 <- NULL

##generate certainty dummy
data23$certain=c()
for (i in 1:9000){
  if (data23$Scale[i]>=7){data23$certain[i]=1}
  else if (data23$Scale[i]<7){data23$certain[i]=0}
} 
##check uncertainty dummy
data23$uncertain=c()
for (i in 1:9000){
  if (data23$Scale[i]<7){data23$uncertain[i]=1}
  else if (data23$Scale[i]>=7){data23$uncertain[i]=0}
} 
##check 
mean(data23$certain)
sum(data23$Scale>=7)/9000

data23$id=c()
for (i in 1:9000){
  data23$id[i]=ceiling(i/3)
}
for (i in 1:9000){
  if (data23$altern[i]==1){data23$altern[i]="Option1"}
  else if (data23$altern[i]==2){data23$altern[i]="Option2"}
  else if (data23$altern[i]==3){data23$altern[i]="Option3"}}
head(data23)
tail(data23)  ##number of obs: 3000, there are 3000 group id

TM_23 <- mlogit.data(data23,choice="Choice",shape="long",alt.var="altern",
                     chid.var="id",id.var="PersonID")
head(TM_23)

rplc_or_23 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                   data=TM_23,
                   ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                   R=500, 
                   #halton=NA,
                   panel=TRUE,
                   correlation=TRUE,
                   model="mixl")
summary(rplc_or_23)
-1*wtp.gmnl(rplc_or_23,wrt="price")
cor.gmnl(rplc_or_23)

wtp_data23<-mlogit.data(data23,choice="Choice",shape="long",alt.var="altern",
                        chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data23)

####Bopt is random
wtp_rplc_or_23<- -1*wtp.gmnl(rplc_or_23,wrt="price")
wtp_rplc_or_23<-data.frame(wtp_rplc_or_23)

start_or_23_1 <- c(1,wtp_rplc_or_23$Estimate[1:4],rep(0.1,11),0.2,0) 
start_or_23_2 <- c(1,rep(0.1,15),0.1,0) 
start_or_23_3 <- c(1,wtp_rplc_or_23$Estimate[1:4],rep(0.1,11),0,0) 
start_or_23_4 <- c(1,rep(1,15),1,0)
start_or_23_5 <- c(1,rep(-1,15),-1,0)

start_23=rbind(start_or_23_1,start_or_23_2,start_or_23_3,
               start_or_23_4,start_or_23_5)

wtpc_or_23 <-list()
ll_23 <- list()

for (i in 1:5){
  wtpc_or_23[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|0|0|0|certain-1,
                          data=wtp_data23,
                          ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                          R=500,
                          #halton = NA,
                          print.init = TRUE, 
                          panel=TRUE,
                          start=start_23[i,],
                          fixed=c(TRUE,rep(FALSE,16),TRUE), 
                          model="gmnl",
                          correlation=TRUE,
                          notscale=c(0,0,0,0,1),
                          method="bhhh",
                          iterlim=500)
  ll_23[[i]]=getSummary.gmnl(wtpc_or_23[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_23

summary(wtpc_or_23[[3]])
cor.gmnl(wtpc_or_23[[3]])
getSummary.gmnl(wtpc_or_23[[3]],alpha=0.05)


#---------------------------------------------------------------------------
#-----------------------------  POOLED 2 & 4  ------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")

##read treatment 2 data
mydata2<-read.csv('Dried_Total_2.csv', fileEncoding = 'UTF-8-BOM')
head(mydata2)
dim(mydata2)
##read treatment 4 data
mydata4<-read.csv('Dried_Total_4.csv', fileEncoding = 'UTF-8-BOM')
head(mydata4)
dim(mydata4)

data24 <- rbind(mydata2,mydata4)
dim(data24)

data24$cset <- NULL
data24$Block <- NULL
data24$Option1 <- NULL
data24$Option2 <- NULL
data24$Option3 <- NULL

##generate certainty dummy
data24$certain=c()
for (i in 1:9000){
  if (data24$Scale[i]>=7){data24$certain[i]=1}
  else if (data24$Scale[i]<7){data24$certain[i]=0}
} 
##check uncertainty dummy
data24$uncertain=c()
for (i in 1:9000){
  if (data24$Scale[i]<7){data24$uncertain[i]=1}
  else if (data24$Scale[i]>=7){data24$uncertain[i]=0}
} 
##check 
mean(data24$certain)
sum(data24$Scale>=7)/9000

data24$id=c()
for (i in 1:9000){
  data24$id[i]=ceiling(i/3)
}
for (i in 1:9000){
  if (data24$altern[i]==1){data24$altern[i]="Option1"}
  else if (data24$altern[i]==2){data24$altern[i]="Option2"}
  else if (data24$altern[i]==3){data24$altern[i]="Option3"}}
head(data24)
tail(data24)  ##number of obs: 3000, there are 3000 group id

TM_24 <- mlogit.data(data24,choice="Choice",shape="long",alt.var="altern",
                     chid.var="id",id.var="PersonID")
head(TM_24)

rplc_or_24 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                   data=TM_24,
                   ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                   R=500, 
                   #halton=NA,
                   panel=TRUE,
                   correlation=TRUE,
                   model="mixl")
summary(rplc_or_24)
-1*wtp.gmnl(rplc_or_24,wrt="price")
cor.gmnl(rplc_or_24)

wtp_data24<-mlogit.data(data24,choice="Choice",shape="long",alt.var="altern",
                        chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data24)

####Bopt is random
wtp_rplc_or_24<- -1*wtp.gmnl(rplc_or_24,wrt="price")
wtp_rplc_or_24<-data.frame(wtp_rplc_or_24)

start_or_24_1 <- c(1,wtp_rplc_or_24$Estimate[1:4],rep(0.1,11),0.2,0) 
start_or_24_2 <- c(1,rep(0.1,15),0.1,0) 
start_or_24_3 <- c(1,wtp_rplc_or_24$Estimate[1:4],rep(0.1,11),0,0) 
start_or_24_4 <- c(1,rep(1,15),1,0)
start_or_24_5 <- c(1,rep(-1,15),-1,0)

start_24=rbind(start_or_24_1,start_or_24_2,start_or_24_3,
               start_or_24_4,start_or_24_5)

wtpc_or_24 <-list()
ll_24 <- list()

for (i in 1:5){
  wtpc_or_24[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|0|0|0|certain-1,
                          data=wtp_data24,
                          ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                          R=500,
                          #halton = NA,
                          print.init = TRUE, 
                          panel=TRUE,
                          start=start_24[i,],
                          fixed=c(TRUE,rep(FALSE,16),TRUE), 
                          model="gmnl",
                          correlation=TRUE,
                          notscale=c(0,0,0,0,1),
                          method="bhhh",
                          iterlim=500)
  ll_24[[i]]=getSummary.gmnl(wtpc_or_24[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_24

summary(wtpc_or_24[[3]])
cor.gmnl(wtpc_or_24[[3]])
getSummary.gmnl(wtpc_or_24[[3]],alpha=0.05)

#---------------------------------------------------------------------------
#-----------------------------  POOLED 3 & 4  ------------------------------
#---------------------------------------------------------------------------

library(plyr)
library(Formula)
library(gmnl)
library(mlogit)
set.seed(1118)

setwd("C:/Users/m_xue/Desktop/RA/Cranberry-gene/data clean/for R")

##read treatment 3 data
mydata3<-read.csv('Dried_Total_3.csv', fileEncoding = 'UTF-8-BOM')
head(mydata3)
dim(mydata3)
##read treatment 4 data
mydata4<-read.csv('Dried_Total_4.csv', fileEncoding = 'UTF-8-BOM')
head(mydata4)
dim(mydata4)

data34 <- rbind(mydata3,mydata4)
dim(data34)

data34$cset <- NULL
data34$Block <- NULL
data34$Option1 <- NULL
data34$Option2 <- NULL
data34$Option3 <- NULL

##generate certainty dummy
data34$certain=c()
for (i in 1:9000){
  if (data34$Scale[i]>=7){data34$certain[i]=1}
  else if (data34$Scale[i]<7){data34$certain[i]=0}
} 
##check uncertainty dummy
data34$uncertain=c()
for (i in 1:9000){
  if (data34$Scale[i]<7){data34$uncertain[i]=1}
  else if (data34$Scale[i]>=7){data34$uncertain[i]=0}
} 
##check 
mean(data34$certain)
sum(data34$Scale>=7)/9000

data34$id=c()
for (i in 1:9000){
  data34$id[i]=ceiling(i/3)
}
for (i in 1:9000){
  if (data34$altern[i]==1){data34$altern[i]="Option1"}
  else if (data34$altern[i]==2){data34$altern[i]="Option2"}
  else if (data34$altern[i]==3){data34$altern[i]="Option3"}}
head(data34)
tail(data34)  ##number of obs: 3000, there are 3000 group id

TM_34 <- mlogit.data(data34,choice="Choice",shape="long",alt.var="altern",
                     chid.var="id",id.var="PersonID")
head(TM_34)

rplc_or_34 <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|-1,
                   data=TM_34,
                   ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                   R=500, 
                   #halton=NA,
                   panel=TRUE,
                   correlation=TRUE,
                   model="mixl")
summary(rplc_or_34)
-1*wtp.gmnl(rplc_or_34,wrt="price")
cor.gmnl(rplc_or_34)

wtp_data34<-mlogit.data(data34,choice="Choice",shape="long",alt.var="altern",
                        chid.var="id",id.var="PersonID",opposite=c("price"))
head(wtp_data34)

####Bopt is random
wtp_rplc_or_34<- -1*wtp.gmnl(rplc_or_34,wrt="price")
wtp_rplc_or_34<-data.frame(wtp_rplc_or_34)

start_or_34_1 <- c(1,wtp_rplc_or_34$Estimate[1:4],rep(0.1,11),0.2,0) 
start_or_34_2 <- c(1,rep(0.1,15),0.1,0) 
start_or_34_3 <- c(1,wtp_rplc_or_34$Estimate[1:4],rep(0.1,11),0,0) 
start_or_34_4 <- c(1,rep(1,15),1,0)
start_or_34_5 <- c(1,rep(-1,15),-1,0)

start_34=rbind(start_or_34_1,start_or_34_2,start_or_34_3,
               start_or_34_4,start_or_34_5)

wtpc_or_34 <-list()
ll_34 <- list()

for (i in 1:5){
  wtpc_or_34[[i]] <- gmnl(Choice~price+Regular+Gene+BlandFlavor+opt|0|0|0|certain-1,
                          data=wtp_data34,
                          ranp=c(Regular="n",Gene="n",BlandFlavor="n",opt="n"),
                          R=500,
                          #halton = NA,
                          print.init = TRUE, 
                          panel=TRUE,
                          start=start_34[i,],
                          fixed=c(TRUE,rep(FALSE,16),TRUE), 
                          model="gmnl",
                          correlation=TRUE,
                          notscale=c(0,0,0,0,1),
                          method="bhhh",
                          iterlim=500)
  ll_34[[i]]=getSummary.gmnl(wtpc_or_34[[i]],alpha=0.05)$sumstat["logLik"]
}
ll_34

summary(wtpc_or_34[[3]])
cor.gmnl(wtpc_or_34[[3]])
getSummary.gmnl(wtpc_or_34[[3]],alpha=0.05)


####################################plot
par(mfrow=c(2,2))

plot(wtpc_or_tr1[[5]], effect="ce", par = "Regular",col = "grey",
     main="Treatment 1",xlab="conditional WTP")
plot(wtpc_or_tr2[[2]], effect="ce", par = "Regular",col = "grey",
     main="Treatment 2",xlab="conditional WTP")
plot(wtpc_or_tr3[[4]], effect="ce", par = "Regular",col = "grey",
     main="Treatment 3",xlab="conditional WTP")
plot(wtpc_or_tr4[[5]], effect="ce", par = "Regular",col = "grey",
     main="Treatment 4",xlab="conditional WTP")


