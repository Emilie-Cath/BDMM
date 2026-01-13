rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/3eme année/2eme papier/data")
source("Mixing_models_functions.R")
#install.packages("deSolve")
#install.packages("ggplot2")
library("deSolve")
library("ggplot2")
library("ggridges")

#---------------DATA---------------------------
iter=0.01
lambda_plasma=0.33
lambda_bloodcell=0.03
Dates=c(0,92,123,396,426,457,488,547)
Sources_Carbon=data.frame(Date=Dates,
                          Zos=rep(-11.17,8),
                          Grass=rep(-30.88,8),
                          Ulvaxe=rep(-12.61,8))
Sources_Nitrogen=data.frame(Date=Dates,
                            Zos=rep(6.49,8),
                            Grass=rep(4.43,8),
                            Ulvaxe=rep(10.50,8))
liste_sources=list(Sources_Carbon,Sources_Nitrogen)

Conc_Carbon=data.frame(Zos=rep(0.35,8),Grass=rep(0.4,8),Ulvaxe=rep(0.197,8))
Conc_Nitrogen=data.frame(Zos=rep(0.03,8),Grass=rep(0.035,8),Ulvaxe=rep(0.0165,8))
liste_Conc=list(Conc_Carbon,Conc_Nitrogen)



TDF_Carbon=data.frame(Zos=1.63,Grass=1.63,Ulvaxe=1.63)
TDF_Nitrogen=data.frame(Zos=3.54,Grass=3.54,Ulvaxe=3.54)
liste_TDF=list(TDF_Carbon,TDF_Nitrogen)



#blood cells

Conso_bloodcell_Carbon=data.frame(Date=Dates,
                                  Iso_C=c(-11.59,-15.97,-21.25,-13.16,-14.83,-20.42,-26.43,-27.57),
                                  lambda=rep(lambda_bloodcell,8))
Conso_bloodcell_Nitrogen=data.frame(Date=Dates,
                                    Iso_C=c(10.30,11.62,10.55,10.75,12.12,10.76,8.86,8.21),
                                    lambda=rep(lambda_bloodcell,8))
liste_Conso_bloodcell=list(Conso_bloodcell_Carbon,Conso_bloodcell_Nitrogen)

Sources_C=c(approxfun(Sources_Carbon$Date,Sources_Carbon$Zos,method='linear',rule=2),
            approxfun(Sources_Carbon$Date,Sources_Carbon$Grass,method='linear',rule=2),
            approxfun(Sources_Carbon$Date,Sources_Carbon$Ulvaxe,method='linear',rule=2))
Sources_N=c(approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Zos,method='linear',rule=2),
            approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Grass,method='linear',rule=2),
            approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Ulvaxe,method='linear',rule=2))
liste_approx=list(Sources_C,Sources_N)

#blood cells
liste_lambda=list(data.frame(Conso_bloodcell_Carbon$Date,Conso_bloodcell_Carbon$lambda),data.frame(Conso_bloodcell_Nitrogen$Date,Conso_bloodcell_Nitrogen$lambda))
#------------fait tourner---------------
#fait tourner
res=all_results(sources_list=liste_sources,data_feat=data_features(liste_sources),consu_list=liste_Conso_bloodcell
                ,conc_list=liste_Conc,iter=iter,list_TEF=liste_TDF,X=50,sources_approx=liste_approx,lambda_list=liste_lambda)

#biplot visuals
sign_biplot(consu_list=liste_Conso_bloodcell,sources_list=liste_sources,list_TEF=liste_TDF,title="",range_C=c(-30,-8),range_N=c(7,15),source_names=c("Zos","grass","Ulva","Entero"))
poly_sources_evol(consu_list=liste_Conso_bloodcell,sources_list=liste_sources,list_TEF=liste_TDF,sources_name=c("Zos","grass","Ulva","Entero"),range_C=c(-30,-8),range_N=c(7,15))


#C
traj_graph(data_feat = data_features(liste_sources),sources_list=liste_sources,sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=1,lambda_list=liste_lambda
           ,conc_list=liste_Conc,list_TEF=liste_TDF,iter=iter,X=50,time_range=c(0,Dates[length(Dates)]),range_y=c(-30,10),title="C",y_label=expression(paste(" δ"^{13},"C (‰)")))
#N
traj_graph(data_feat = data_features(liste_sources),sources_list=liste_sources,sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=2,lambda_list=liste_lambda
           ,conc_list=liste_Conc,list_TEF=liste_TDF,iter=iter,X=50,time_range=c(0,Dates[length(Dates)]),range_y=c(-10,15),title="N",y_label=expression(paste(" δ"^{15},"N (‰)")))


contrib_time(res_model=res[[3]],model_name="DMM",date=Dates
             ,data_feat=data_features(liste_sources),X=50,sources_name=c("Zos","grass","Ulva","Entero"),where_legend="topleft")


##----CONTRIB GRAPH LINEAR----
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

tiff("Contrib_geese_3_s_DMM.tiff",  width =4 , height = 5, units = "in", res = 300)
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

tiff("legend_geese_3_s.tiff",  width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",legend=c("Zos","grass","Ulva x Entero"),cex=1,col=c(color_main[1:3]),lty=1,box.lty=0,lwd=2)
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

tiff("Contrib_geese_3_s_DMM_constant.tiff",  width =4 , height = 5, units = "in", res = 300)
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


#------source poly--------

library(dplyr)
library(tidyr)
library(ggplot2)


# Open the geese data
library(readxl)
geese_dat <- read_excel("geese_data_3_sources.xls")
geese_sources <- read_excel("geese_data_3_sources.xls", sheet = 'Sources')
sources_meanC<- read_excel("geese_data_3_sources.xls", sheet = 'Sources_MeanC')%>% select(-Time) %>% as.matrix()
sources_meanN<- read_excel("geese_data_3_sources.xls", sheet = 'Sources_MeanN')%>% select(-Time) %>% as.matrix()
sources_sdC<- read_excel("geese_data_3_sources.xls", sheet = 'Sources_sdC')%>% select(-Time) %>% as.matrix()
sources_sdN<- read_excel("geese_data_3_sources.xls", sheet = 'Sources_sdN')%>% select(-Time) %>% as.matrix()
geese_TEF <- read_excel("geese_data_3_sources.xls", sheet = 'TEFs')
TEF_meanC<- read_excel("geese_data_3_sources.xls", sheet = 'TEF_MeanC')%>% select(-Time) %>% as.matrix()
TEF_meanN<- read_excel("geese_data_3_sources.xls", sheet = 'TEF_MeanN')%>% select(-Time) %>% as.matrix()
TEF_sdC<- read_excel("geese_data_3_sources.xls", sheet = 'TEF_sdC')%>% select(-Time) %>% as.matrix()
TEF_sdN<- read_excel("geese_data_3_sources.xls", sheet = 'TEF_sdN')%>% select(-Time) %>% as.matrix()
geese_CD <- read_excel("geese_data_3_sources.xls", sheet = 'ConcDep')

lambda_mean=0.03
lambda_sd=0.005
alpha_lambda=lambda_mean/lambda_sd^2
beta_lambda=lambda_mean^2/lambda_sd^2
#t_obs <- sort(jitter(oysters_dat$Time / 100, amount = 0.01))
t_obs <- sort(jitter(geese_dat$Time+0.01 / 1, amount = 0.01))
y_obs <- cbind(geese_dat$d13C_Pl, geese_dat$d15N_Pl) # ADD NOISE HERE
# t_fit <- seq(min(t_obs), max(t_obs), length = 20) # Can't do t_fit on a new grid unless we interpolate the source values
t_fit <- t_obs

bdmm_data <- list(
  n_obs = length(t_obs),
  n_fit = length(t_fit),
  K = ncol(sources_meanC),
  D = 3,
  y_obs = y_obs,
  t_obs = t_obs,
  t_fit = t_fit,
  t0 = 0,
  mu_s1 = sources_meanC + TEF_meanC,
  mu_s2 = sources_meanN + TEF_meanN,
  sigma_s1 = sqrt(sources_sdC^2 + TEF_sdC^2),
  sigma_s2 = sqrt(sources_sdN^2 + TEF_sdN^2),
  alpha_lambda = as.numeric(alpha_lambda), 
  beta_lambda = as.numeric(beta_lambda),
  q1 = geese_CD$d13CPl, 
  q2 = geese_CD$d15NPl, 
  t_mean = mean(t_obs),
  t_sd =  sd(t_obs)
)



####GRAPHS source poly ######
temps=unique(geese_dat$Time)
sources_corr_C=data.frame(cbind(Time=geese_dat$Time,bdmm_data$mu_s1))
sources_corr_N=data.frame(cbind(Time=geese_dat$Time,bdmm_data$mu_s2))
sources_corr_sd_C=data.frame(cbind(Time=geese_dat$Time,bdmm_data$sigma_s1))
sources_corr_sd_N=data.frame(cbind(Time=geese_dat$Time,bdmm_data$sigma_s2))


for(i in 1:length(temps)){
  
  t=temps[i]
  subs_conso=subset(geese_dat,geese_dat$Time==t)
  s_corr_C=subset(sources_corr_C,sources_corr_C$Time==t)[1,]
  s_corr_N=subset(sources_corr_N,sources_corr_N$Time==t)[1,]
  sd_corr_C=subset(sources_corr_sd_C,sources_corr_sd_C$Time==t)[1,]
  sd_corr_N=subset(sources_corr_sd_N,sources_corr_sd_N$Time==t)[1,]
  
  title_graph=paste0("geese_3s_biplot_t_",t,'.tiff')  
  tiff(title_graph,  width =4 , height = 5, units = "in", res = 300)
  plot(NULL,main=paste("t=",t,"d"),xlim=c(-30,-8),ylim=c(5,18),xlab="",ylab="",cex.lab=1.4,font.lab=2,las=1) #plot frame
  # mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,cex=1.5,font=2)
  # mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=3,cex=1.5,font=2)
  
  
  
  for (j in 2:ncol(s_corr_C)){
    # j=2
    arrows(x0=s_corr_C[[j]],y0=s_corr_N[[j]]-sd_corr_N[[j]]
           ,x1=s_corr_C[[j]],y1=s_corr_N[[j]]+sd_corr_N[[j]]
           ,length=0,col=color_main[j-1],lwd=3)
    arrows(x0=s_corr_C[[j]]-sd_corr_C[[j]],y0=s_corr_N[[j]]
           ,x1=s_corr_C[[j]]+sd_corr_C[[j]],y1=s_corr_N[[j]]
           ,length=0,col=color_main[j-1],lwd=3)
  }
  points(subs_conso$d13C_Pl,subs_conso$d15N_Pl,pch=16,cex=1.5)
  dev.off

}



tiff("legend_biplot_geese_3s_.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,main="",xlim=c(-30,-8),ylim=c(5,18),xlab="",ylab="",cex.lab=1.4,font.lab=2,las=1) #plot frame
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=-1,padj=0.8,cex=1.5,font=2)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=-2,cex=1.5,font=2)
legend("center",legend=c("Zos","grass","Ulva x Entero","Consumer"),cex=1.2,col=c(color_main[1:3],"black")
       ,lty=c(1,1,1,NA),box.lty=0,pch=c(NA,NA,NA,16),lwd=c(2,2,2,NA))
dev.off()
