rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/3eme année/2eme papier/data")
source("Mixing_models_functions.R")
#install.packages("deSolve")
#install.packages("ggplot2")
library("deSolve")
library("ggplot2")
library("ggridges")



#-------------DATA BDV ------------
iter=0.01
Sources_Carbon_BDV<-data.frame(Date=c(0,60,120,180,240,300,360)
                               ,Iso_MPBxU=c(-15.56, -17.60, -15.46, -17.57, -19.85,-20.50,-17.01)
                               ,Iso_Phyto=c(-21.91,-20.47,-22.65,-20.40,-21.31,-22.22,-18.59)
                               ,Iso_TOM=c(-29.35,-26.97,-30.57,-28.50,-28.42,-28.34,-30.00))

Sources_Nitrogen_BDV<-data.frame(Date=c(0,60,120,180,240,300,360)
                                 ,Iso_MPBxU=c(6.14,10.56,7.55,8.42,9.27,9.39,7.17)
                                 ,Iso_Phyto=c(5.46,4.35,6.00,4.25,3.68,3.11,8.15)
                                 ,Iso_TOM=c(7.27,6.43,8.20,3.81,5.16,6.51,9.14))


liste_sources_BDV = list(Sources_Carbon_BDV,Sources_Nitrogen_BDV)


TDF2_Carbon <- data.frame(TDF_MPBxU=0.4,TDF_Phyto=0.4,TDF_TOM=0.4)
TDF2_Nitrogen <- data.frame(TDF_MPBxU=2.2,TDF_Phyto=2.2,TDF_TOM=2.2)
liste_TDF2= list(TDF2_Carbon,TDF2_Nitrogen)

Conso_Carbon_BDV_N <- data.frame(Date=c(0,60,120,180,240,300,360)
                                 ,Iso_conso=c(-19.31,-19.14,-19.11,-19.38,-19.48,-19.97,-19.18)
                                 ,lambda=c(0.027,0.019,0.011,0.004,0.009,0.014,NA))

Conso_Nitrogen_BDV_N <- data.frame(Date=c(0,60,120,180,240,300,360)
                                   ,Iso_conso=c(6.87,8.94,9.94,10.82,10.32,9.07,8.98)
                                   ,lambda=c(0.027,0.019,0.011,0.004,0.009,0.014,NA))
liste_conso_BDV_N = list(Conso_Carbon_BDV_N,Conso_Nitrogen_BDV_N)

conc_Carbon_BDV= data.frame(MPBxU=c(rep(1,7)),Phyto=c(rep(1,7)),TOM=c(rep(1,7)))
conc_Nitrogen_BDV = data.frame(MPBxU=c(rep(1,7)),Phyto=c(rep(1,7)),TOM=c(rep(1,7)))
liste_conc_BDV=list(conc_Carbon_BDV,conc_Nitrogen_BDV)



Sources_C= c(approxfun(Sources_Carbon_BDV$Date,Sources_Carbon_BDV$Iso_MPBxU,method='linear',rule=2),
             approxfun(Sources_Carbon_BDV$Date,Sources_Carbon_BDV$Iso_Phyto,method='linear',rule=2),
             approxfun(Sources_Carbon_BDV$Date,Sources_Carbon_BDV$Iso_TOM,method='linear',rule=2))

Sources_N=c(approxfun(Sources_Nitrogen_BDV$Date,Sources_Nitrogen_BDV$Iso_MPBxU,method="linear",rule=2),
            approxfun(Sources_Nitrogen_BDV$Date,Sources_Nitrogen_BDV$Iso_Phyto,method="linear",rule=2),
            approxfun(Sources_Nitrogen_BDV$Date,Sources_Nitrogen_BDV$Iso_TOM,method="linear",rule=2))
liste_approx=list(Sources_C,Sources_N)
liste_lambda_N=list(data.frame(Conso_Carbon_BDV_N$Date,Conso_Carbon_BDV_N$lambda),data.frame(Conso_Nitrogen_BDV_N$Date,Conso_Nitrogen_BDV_N$lambda))

#----------RUNNING DMM---------------

#fait tourner
res=all_results(sources_list=liste_sources_BDV,data_feat=data_features(liste_sources_BDV),consu_list=liste_conso_BDV_N,conc_list=liste_conc_BDV,iter=iter
                ,list_TEF=liste_TDF2,X=50,sources_approx=liste_approx,lambda_list=liste_lambda_N)


#-----Contribution -----

DMM=res[[3]]
temps=unique(DMM[[1]])

stats=data.frame(row.names=NULL)
 for (i in 1:length(temps)){
  # i=1
  t=temps[i]
  subs=subset(DMM,DMM$time==t)
  stats=rbind(stats,data.frame(time=t,
                               mV1=mean(subs$Var1),Q1V1=quantile(subs$Var1,probs =c(0.25)),
                               Q3V1=quantile(subs$Var1,probs =c(0.75)),
                               mV2=mean(subs$Var2),Q1V2=quantile(subs$Var2,probs =c(0.25)),
                               Q3V2=quantile(subs$Var2,probs =c(0.75)),
                               mV3=mean(subs$Var3),Q1V3=quantile(subs$Var3,probs =c(0.25)),
                               Q3V3=quantile(subs$Var3,probs =c(0.75))))
 }


color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
             rgb(0,0.6,0.45),rgb(0.9,0.6,0),rgb(0.8,0.47,0.65),
             rgb(0.33,0.7,0.9),rgb(0.94,0.89,0.25))
color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),
             rgb(0.9,0.6,0,0.2),rgb(0.8,0.47,0.65,0.2),rgb(0.33,0.7,0.9,0.2),rgb(0.94,0.89,0.25,0.2) )

tiff("Contrib_oysters_DMM.tiff",  width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(NULL,xlim=c(0,tail(stats$time,1)),ylim=c(0,1),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V1,rev(stats$Q3V1))
        ,col=color_poly[1],border=NA)
polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V2,rev(stats$Q3V2))
        ,col=color_poly[2],border=NA)
polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V3,rev(stats$Q3V3))
        ,col=color_poly[3],border=NA)


points(stats$time,stats$mV1,type="l",col=color_main[1],pch=16,lwd=2)
points(stats$time,stats$mV2,type="l",col=color_main[2],pch=16,lwd=2)
points(stats$time,stats$mV3,type="l",col=color_main[3],pch=16,lwd=2)
dev.off()

tiff("legend_oyster.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",legend=c("MPBxU","POM","TOM"),cex=1,col=c(color_main[1:3]),lty=1,box.lty=0,lwd=2)
dev.off()


##----CONTRIB constant interpolation----
stats_bis=data.frame(row.names=NULL)

for (i in 1:(nrow(stats)-1)){
  stats_bis=rbind(stats_bis,data.frame(stats[i,]))
  stats_bis=rbind(stats_bis,data.frame(time=(stats[(i+1),1]-1),stats[i,-1]))
}

stats_bis=rbind(data.frame(time=c(0),stats[1,-1]),stats_bis)

color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
             rgb(0,0.6,0.45),rgb(0.9,0.6,0),rgb(0.8,0.47,0.65),
             rgb(0.33,0.7,0.9),rgb(0.94,0.89,0.25))
color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),
             rgb(0.9,0.6,0,0.2),rgb(0.8,0.47,0.65,0.2),rgb(0.33,0.7,0.9,0.2),rgb(0.94,0.89,0.25,0.2) )

tiff("Contrib_oysters_3_s_DMM_constant.tiff",  width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
par(xpd=TRUE, mar=c(8,4,4,3))
plot(NULL,xlim=c(0,tail(stats_bis$time,1)),ylim=c(0,1),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(stats_bis$time,rev(stats_bis$time))
        ,y=c(stats_bis$Q1V1,rev(stats_bis$Q3V1))
        ,col=color_poly[1],border=NA)
polygon(x=c(stats_bis$time,rev(stats_bis$time))
        ,y=c(stats_bis$Q1V2,rev(stats_bis$Q3V2))
        ,col=color_poly[2],border=NA)
polygon(x=c(stats_bis$time,rev(stats_bis$time))
        ,y=c(stats_bis$Q1V3,rev(stats_bis$Q3V3))
        ,col=color_poly[3],border=NA)

points(stats_bis$time,stats_bis$mV1,type="l",col=color_main[1],pch=16,lwd=2)
points(stats_bis$time,stats_bis$mV2,type="l",col=color_main[2],pch=16,lwd=2)
points(stats_bis$time,stats_bis$mV3,type="l",col=color_main[3],pch=16,lwd=2)
dev.off()


##----FIT----
#carbon
traj_graph(sources_list=liste_sources_BDV,data_feat=data_features(liste_sources_BDV)
           ,sources_approx=liste_approx,consu_list=liste_conso_BDV_N,iso_num=1
           ,lambda_list=liste_lambda_N,conc_list=liste_conc_BDV
           ,list_TEF=liste_TDF2,iter=iter,X=50,time_range=c(0,360),range_y=c(-22,-17),title="C",y_label=expression(paste(" δ"^{13},"C (‰)")))
#nitrogen
traj_graph(sources_list=liste_sources_BDV,data_feat=data_features(liste_sources_BDV)
           ,sources_approx=liste_approx,consu_list=liste_conso_BDV_N,iso_num=2
           ,lambda_list=liste_lambda_N,conc_list=liste_conc_BDV
           ,list_TEF=liste_TDF2,iter=iter,X=50,time_range=c(0,360),range_y=c(6,11),title="N",y_label=expression(paste(" δ"^{15},"N (‰)")))


##----SMM-----

SMM=res[[1]]
temps=unique(SMM[[1]])

stats=data.frame(row.names=NULL)
for (i in 1:length(temps)){
  # i=1
  t=temps[i]
  subs=subset(SMM,SMM$time==t)
  stats=rbind(stats,data.frame(time=t,
                               mV1=mean(subs$Var1),Q1V1=quantile(subs$Var1,probs =c(0.25)),
                               Q3V1=quantile(subs$Var1,probs =c(0.75)),
                               mV2=mean(subs$Var2),Q1V2=quantile(subs$Var2,probs =c(0.25)),
                               Q3V2=quantile(subs$Var2,probs =c(0.75)),
                               mV3=mean(subs$Var3),Q1V3=quantile(subs$Var3,probs =c(0.25)),
                               Q3V3=quantile(subs$Var3,probs =c(0.75))))
}


color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
             rgb(0,0.6,0.45),rgb(0.9,0.6,0),rgb(0.8,0.47,0.65),
             rgb(0.33,0.7,0.9),rgb(0.94,0.89,0.25))
color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),
             rgb(0.9,0.6,0,0.2),rgb(0.8,0.47,0.65,0.2),rgb(0.33,0.7,0.9,0.2),rgb(0.94,0.89,0.25,0.2) )

tiff("Contrib_oysters_SMM.tiff",  width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(NULL,xlim=c(0,tail(stats$time,1)),ylim=c(0,1),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V1,rev(stats$Q3V1))
        ,col=color_poly[1],border=NA)
polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V2,rev(stats$Q3V2))
        ,col=color_poly[2],border=NA)
polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V3,rev(stats$Q3V3))
        ,col=color_poly[3],border=NA)


points(stats$time,stats$mV1,type="l",col=color_main[1],pch=16,lwd=2)
points(stats$time,stats$mV2,type="l",col=color_main[2],pch=16,lwd=2)
points(stats$time,stats$mV3,type="l",col=color_main[3],pch=16,lwd=2)
dev.off()