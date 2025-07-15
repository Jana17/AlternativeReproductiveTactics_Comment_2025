#########R code to produce the figures in:
#########"Do alternative reproductive tactics facilitate evolutionary rescue? A comment on Knell & Parrett 2024."
#########Jana M. Riederer & Franz J. Weissing
#########j.m.riederer@rug.nl


##########Figure 1###############

##Here, set the working directory which contains the simulation output,
##and in which the figure will be
svg(filename="ART_Fig1.svg", width = 8.80, height = 10.80)

par(mfrow=c(4,2))
l_tr=0.3
writing_size=1.4

Data<-read.csv(paste("ART_Fig1_NoTEvo_initT0.csv",sep = ""))
plot(NULL,xlim=c(0,5000),ylim=c(0,1050),main="No ART",ylab="Population size",xlab="Time",
     cex.lab=writing_size, cex.axis=writing_size, cex.main=writing_size )
for(i in 0:29){
  Data_rep<-subset(Data,Data$Rep==i)
  lines(Data_rep$PopSize~Data_rep$Time,col=rgb(0,0,0,l_tr))
}
text(5000, 960, "A", col = "black", cex = 2, pos=2)


Data<-read.csv(paste("ART_Fig1_NoTEvo_initT0.25.csv",sep = ""))
plot(NULL,xlim=c(0,5000),ylim=c(0,1050),main="T = 0.25",ylab="Population size",xlab="Time", 
     cex.lab=writing_size, cex.axis=writing_size, cex.main=writing_size )
for(i in 0:29){
  Data_rep<-subset(Data,Data$Rep==i)
  lines(Data_rep$PopSize~Data_rep$Time,col=rgb(0,0,0,l_tr))
}
text(5000, 960, "B", col = "black", cex = 2, pos=2)


Data<-read.csv(paste("ART_Fig1_NoTEvo_initT0.2.csv",sep = ""))
plot(NULL,xlim=c(0,5000),ylim=c(0,1050),main="T = 0.2",ylab="Population size",xlab="Time", 
     cex.lab=writing_size, cex.axis=writing_size, cex.main=writing_size )
for(i in 0:29){
  Data_rep<-subset(Data,Data$Rep==i)
  lines(Data_rep$PopSize~Data_rep$Time,col=rgb(0,0,0,l_tr))
}
text(5000, 960, "C", col = "black", cex = 2, pos=2)


Data<-read.csv(paste("ART_Fig1_NoTEvo_initT0.1.csv",sep = ""))
plot(NULL,xlim=c(0,5000),ylim=c(0,1050),main="T = 0.1",ylab="Population size",xlab="Time", 
     cex.lab=writing_size, cex.axis=writing_size, cex.main=writing_size )
for(i in 0:29){
  Data_rep<-subset(Data,Data$Rep==i)
  lines(Data_rep$PopSize~Data_rep$Time,col=rgb(0,0,0,l_tr))
}
text(5000, 960, "D", col = "black", cex = 2, pos=2)


Data<-read.csv(paste("ART_Fig1_TEvo_initT0.25.csv",sep = ""))
plot(NULL,xlim=c(0,5000),ylim=c(0,0.45),main="Evolving T",ylab="Average T",xlab="Time", 
     cex.lab=writing_size, cex.axis=writing_size, cex.main=writing_size )
for(i in 0:29){
  Data_rep<-subset(Data,Data$Rep==i)
  lines(Data_rep$avg_thr~Data_rep$Time,col=rgb(0,0,0,l_tr))
}
text(5000, 0.4114286, "E", col = "black", cex = 2, pos=2)

plot(NULL,xlim=c(0,5000),ylim=c(0,1050),main="Evolving T",ylab="Population size",xlab="Time", 
     cex.lab=writing_size, cex.axis=writing_size, cex.main=writing_size )
for(i in 0:29){
  Data_rep<-subset(Data,Data$Rep==i)
  lines(Data_rep$PopSize~Data_rep$Time,col=rgb(0,0,0,l_tr))
}
text(5000, 960, "F", col = "black", cex = 2, pos=2)


Data<-read.csv(paste("ART_Fig1_NoTEvo_initT0.25_alternativeMatingSystem.csv",sep = ""))
plot(NULL,xlim=c(0,5000),ylim=c(0,1050),main="New mating scenario, 
     T = 0.25",ylab="Population size",xlab="Time", 
     cex.lab=writing_size, cex.axis=writing_size, cex.main=writing_size )
for(i in 0:29){
  Data_rep<-subset(Data,Data$Rep==i)
  lines(Data_rep$PopSize~Data_rep$Time,col=rgb(0,0,0,l_tr))
}
text(5000, 960, "G", col = "black", cex = 2, pos=2)


Data<-read.csv(paste("ART_Fig1_TEvo_initT0.25_alternativeMatingSystem.csv",sep = ""))
plot(NULL,xlim=c(0,5000),ylim=c(0,1050),main="New mating scenario, 
     Evolving T",ylab="Population size",xlab="Time", 
     cex.lab=writing_size, cex.axis=writing_size, cex.main=writing_size )
for(i in 0:29){
  Data_rep<-subset(Data,Data$Rep==i)
  lines(Data_rep$PopSize~Data_rep$Time,col=rgb(0,0,0,l_tr))
}
text(5000, 960, "H", col = "black", cex = 2, pos=2)

dev.off()



###########Figure 2#################

##Here, set the working directory which contains the simulation output,
##and in which the figure will be
AllData1<-data.frame()
for(i in 1:20){
  Data<-read.csv(paste("VarInflux_ART_Fig2_NoTEvo_InheritBlending_many_part",i,".csv",sep = ""))
  MutStDev=Data$MStdDev[1]
  for(j in 1:50){
    MutRate=Data$MRate[100*j]
    dC<-subset(Data,Data$MRate==MutRate)
    if(nrow(dC)==0){print("Error!")
      print(MutStDev)
    }
    dC0<-subset(dC,dC$iniT==0)
    dC1<-subset(dC,dC$iniT==0.25)
    dC0_ex<-subset(dC0,dC0$Extinction==1)
    dC1_ex<-subset(dC1,dC1$Extinction==1)
    nr_ex_0<-nrow(dC0_ex)
    nr_ex_1<-nrow(dC1_ex)
    ART_Window=1
    if(nr_ex_0==nr_ex_1){
      if(nr_ex_0<50){
        ART_Window=0
        if(nr_ex_0!=0){print("Error!")}
        }
      else{
        ART_Window=2
        if(nr_ex_0!=50){print("Error!")}
        }
    }
    if(nrow(dC)!=0){
      newrow<-cbind(MutRate,MutStDev,nr_ex_0,nr_ex_1,ART_Window)
      AllData1<-rbind(AllData1,newrow)
    }
  }
}


AllData2<-data.frame()
for(i in 1:20){
  Data<-read.csv(paste("VarInflux_ART_Fig2_NoTEvo_InheritDipl_many_part",i,".csv",sep = ""))
  MutStDev=Data$MStdDev[1]
  for(j in 1:50){
    MutRate=Data$MRate[100*j]
    dC<-subset(Data,Data$MRate==MutRate)
    if(nrow(dC)==0){print("Error!")
      print(MutStDev)
    }
    dC0<-subset(dC,dC$iniT==0)
    dC1<-subset(dC,dC$iniT==0.25)
    dC0_ex<-subset(dC0,dC0$Extinction==1)
    dC1_ex<-subset(dC1,dC1$Extinction==1)
    nr_ex_0<-nrow(dC0_ex)
    nr_ex_1<-nrow(dC1_ex)
    ART_Window=1
    if(nr_ex_0==nr_ex_1){
      if(nr_ex_0<50){
        ART_Window=0
        if(nr_ex_0!=0){print("Error!")}
        }
      else{
        ART_Window=2
        if(nr_ex_0!=50){print("Error!")}
        }
    }
    if(nrow(dC)!=0){
      newrow<-cbind(MutRate,MutStDev,nr_ex_0,nr_ex_1,ART_Window)
      AllData2<-rbind(AllData2,newrow)
    }
  }
}

AllDataB_W0<-subset(AllData1,AllData1$ART_Window==0)##zero extinctions
AllDataB_W2<-subset(AllData1,AllData1$ART_Window==2)##all extinctions
AllDataB_W1<-subset(AllData1,AllData1$ART_Window==1)

AllDataD_W0<-subset(AllData2,AllData2$ART_Window==0)##zero extinctions
AllDataD_W2<-subset(AllData2,AllData2$ART_Window==2)##all extinctions
AllDataD_W1<-subset(AllData2,AllData2$ART_Window==1)

svg(filename="ART_Fig2.svg",width = 6.30, height = 6.30)
par(mfrow=c(1,2))
plot(NULL,xlim=c(0,0.105),ylim=c(0,1),xaxt='n',
     xlab="Mutational standard deviation", ylab="Mutation rate",main="Blending inheritance")
axis(1, at=c(0,0.05,0.1), labels=TRUE)
points(AllDataB_W0$MutRate~AllDataB_W0$MutStDev,col="aquamarine2",pch=15,cex=1)
points(AllDataB_W2$MutRate~AllDataB_W2$MutStDev,col='dark grey',pch=15,cex=1)
points(AllDataB_W1$MutRate~AllDataB_W1$MutStDev,col=rgb(1,1-abs(AllDataB_W1$nr_ex_0-AllDataB_W1$nr_ex_1)/50,0),
       pch=15,cex=1)
text(0.05, 1, "*", cex = 1.5)

plot(NULL,xlim=c(0,0.105),ylim=c(0,1),xaxt='n',
     xlab="Mutational standard deviation", ylab="Mutation rate",main="Mendelian inheritance")
axis(1, at=c(0,0.05,0.1), labels=TRUE)
points(AllDataD_W0$MutRate~AllDataD_W0$MutStDev,col="aquamarine2",pch=15,cex=1)
points(AllDataD_W2$MutRate~AllDataD_W2$MutStDev,col='dark grey',pch=15,cex=1)
points(AllDataD_W1$MutRate~AllDataD_W1$MutStDev,col=rgb(1,1-abs(AllDataD_W1$nr_ex_0-AllDataD_W1$nr_ex_1)/50,0),
       pch=15,cex=1)

dev.off()



