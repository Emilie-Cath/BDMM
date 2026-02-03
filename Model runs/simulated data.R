

rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/3eme année/2eme papier/data")
source("Mixing_models_functions.R")
library(deSolve)
library(philentropy)
library("viridisLite")
library(rstatix)
####FUNCTION###
stat_making<- function(samp_time,lambda,real_mean_diet,MM){
  stats_MM=data.frame(row.names = NULL)
  for (j in 2:nrow(samp_time)){
    # j=2
    times=samp_time[j,]
    lam_T=(samp_time[j,]-samp_time[j-1,])*lambda
    subs_real_diet=subset(real_mean_diet,real_mean_diet$samp_3==times)
    subs_MM=subset(MM,MM$time==times)
    bias=(abs(subs_real_diet$V1-subs_MM$Var1)+abs(subs_real_diet$V2-subs_MM$Var2)+abs(subs_real_diet$V3-subs_MM$Var3))*50
    bias_stats=quantile(bias,probs=c(0.25,0.5,0.75))
    stats_MM=rbind(stats_MM,data.frame(time=times,lambda_T=lam_T
                                       ,B_Q1=bias_stats[1],B_M=bias_stats[2],B_Q3=bias_stats[3]))
    
  }
  return(stats_MM)
}

#----GENERATING THE ENVIRONMENT-----
period=500
t <- seq(0,period+1,1) # output times for the ODE (d)

#CHANGE HERE IF THE SOURCE NUMBER OR ISOTOPE NUMBER CHANGE
TEF=list(data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0)),data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0))) #the TEF=0 for all sources and isotope
conc=list(data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1)),data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1))) #the conc=1 for all sources and isotope
iter=0.05
X=11 #Number of sols kept by the model here 5%

####----LAMBDA VALUE----
#lambda system
lambda=0.02
sd_lambda_syst=0.001
#lambda param
lambda_data=0.2
sd_lambda_data=0.1

liste_lambda=list(data.frame(Date=t,l=rep(lambda,period+2)),data.frame(Date=t,l=rep(lambda,period+2)))
#WARNING IF T IS TOO SHORT COMPARED TO LAMBDA
if((2*log(2))/lambda >period)warning(paste("period should be at least",(2*log(2))/lambda,"days"))




#GENERATING REGULAR SWITCHS in diet
pB_func <- function(t,b,a,k) {a*sin(t*b)+k} 

b=0.025#switch freq
a=0.5 #amplitude
k=0.5 #oscillation axis
combis2=pB_func(t,b,a,k)
combis3=1-combis2

# # windows()
# plot(NULL,xlim=c(0,period+2),ylim=c(0,2),xlab="Time (d)",ylab="Contribution")
# points(t,rep(0,period+2),type="l",lwd=1.4,col="black",lty=2)
# points(t,combis2,type="l",lwd=1.4,col="red")
# points(t,combis3,type="l",lwd=1.4,col="green")
# legend("topright",legend=c("Source1","Source2","Source3"),cex=1,col=c("black","red","green"),lty=c(2,1,1))


#SOURCES VALUES
#all sources constant over time
# sources_list=list(data.frame(Date=t,s1=rep(10,period+2),s2=rep(0,period+2),s3=rep(5,period+2)),
#                   data.frame(Date=t,s1=rep(5,period+2),s2=rep(0,period+2),s3=rep(10,period+2)))

#SOURCES WITH RANDOM VARIATION
#GENERATING
# s3=c(10)
# for (i in 1:(period+1)) {
#   n=rnorm(1,0,0.1)
#   s3=c(s3,s3[i]+n)
# }
# plot(t,s3)
# write.csv2(s3,"s3_N_DMM_degradated.csv")
s2_C=read.csv2("s2_C_DMM_degradated.csv")[-1][[1]]
s3_C=read.csv2("s3_C_DMM_degradated.csv")[-1][[1]]
s2_N=read.csv2("s2_N_DMM_degradated.csv")[-1][[1]]
s3_N=read.csv2("s3_N_DMM_degradated.csv")[-1][[1]]

sources_list=list(data.frame(Date=t,s1=rep(10,period+2)
                             ,s2=s2_C[1:(period+2)]
                             ,s3=s3_C[1:(period+2)]),
                  data.frame(Date=t,s1=rep(-5,period+2),
                             s2=s2_N[1:(period+2)],
                             s3=s3_N[1:(period+2)]))


sources_approxfun=list(c(approxfun(sources_list[[1]]$Date,sources_list[[1]]$s1,method='linear',rule=2)
                         ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s2,method='linear',rule=2)
                         ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s3,method='linear',rule=2)),
                       c(approxfun(sources_list[[2]]$Date,sources_list[[2]]$s1,method='linear',rule=2)
                         ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s2,method='linear',rule=2)
                         ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s3,method='linear',rule=2)))


data_feat=data_features(sources_list )
#CONSUMER INITAL VALUE
conso_ini=c(5,10)

#GENERATING δX and δXinf
res=NULL
for (n in 1:data_feat[1]){ #data_feat[1] is the number of isotopes
  res=c(res,list(RUN_TIMdyn2(t,#time
                             state0=c(X=conso_ini[n]),#initial conditions
                             par=c(#two parameters lambda and Xinf i.e. forcing functions 
                               lambda = approxfun(liste_lambda[[n]]$Date,liste_lambda[[n]]$l,method = 'const',rule=2), 
                               Xinf=approxfun(t, Mixing(t=t,sources_approx=sources_approxfun[[n]],p=data.frame(c(0),approxfun(t,combis2)(t),approxfun(t,combis3)(t))
                                                        ,TEF=TEF[[n]],Conc=conc[[n]])))))) 
}

#----SAMPLING VALUE IN THE ENV------
#CONSUMER VALUE
nb_val=10
# samp_time=c(sample(1:period,nb_val))
# samp_time=as.numeric(samp_time)
# 
# samp_2=order(samp_time,decreasing = FALSE)
# samp_3=samp_time[samp_2]
# samp_time=data.frame(samp_3)

samp_time=data.frame(samp_3=c(8,59,72,84,123,140,291,344,462,471))
#extracting the consumer values at the sampling times
samp_C=data.frame(row.names = NULL)
samp_N=data.frame(row.names = NULL)
for (i in 1:nrow(samp_time)){
  samp_C=rbind(samp_C,subset(res[[1]],res[[1]]$time==samp_time[i,]))
  samp_N=rbind(samp_N,subset(res[[2]],res[[2]]$time==samp_time[i,]))
}

#Source value
samp_source=data.frame(samp_time)
for (i in 1:length(sources_approxfun)){
  for (j in 1:length(sources_approxfun[[i]])){
    l=NULL
    for (k in 1:nrow(samp_time)){
      l=c(l,sources_approxfun[[i]][[j]](samp_time[k,]))
    }
    samp_source=cbind(samp_source,data.frame(l))}}


# # windows()
# dev.off() # ferme les anciennes fenêtes
# par(fig=c(0,0.5,0,1), new=TRUE)
# plot(NA, NA, ylim = c(-10,15), xlim = c(0,period),
#      xlab = c("Time (d) "), ylab ="" , las=1, cex.lab=1.4,font.lab=1, cex.axis=1.4,main="")
# mtext(expression(paste(" δ"^{13},"C (‰)")),side=2,line=4,padj=1.2,cex=1.3,font=2)
# points(res[[1]]$time,sources_list[[1]]$s1,type="l",col="red",lty=2,lwd=2)
# points(res[[1]]$time,res[[1]]$X,col="aquamarine",type="l",lwd=3)
# points(res[[1]]$time,sources_list[[1]]$s2,type="l",col="red",lwd=2)
# points(res[[1]]$time,sources_list[[1]]$s3,type="l",col="red",lwd=2)
# points(samp_C$time,samp_C$X,type="p",pch=16,col="black")
# 
# par(fig=c(0.5,1,0,1), new=TRUE)
# plot(NA, NA, ylim = c(-10,15), xlim = c(0,period),
#      xlab = c("Time (d) "), ylab ="" , las=1, cex.lab=1.4,font.lab=1, cex.axis=1.4,main="")
# mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.2,cex=1.3,font=2)
# # points(res[[2]]$time,res[[2]]$Xinf,lty=3,type="l",col="black"
# points(res[[2]]$time,sources_list[[2]]$s1,type="l",col="red",lty=2,lwd=2)
# points(res[[2]]$time,res[[2]]$X,col="aquamarine",type="l",lwd=3)
# points(res[[2]]$time,sources_list[[2]]$s2,type="l",col="red",lwd=2)
# points(res[[2]]$time,sources_list[[2]]$s3,type="l",col="red",lwd=2)
# points(samp_N$time,samp_N$X,type="p",pch=16,col="black")
# 
# 
# plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
# legend("center",legend=c("Consumer","Source1","Source2-3","Sampled consumer")
#        ,col=c("aquamarine","red","red","black"),lty=c(1,2,1,NA),pch=c(NA,NA,NA,16),cex=1.1,
#        lwd=3, inset = 0.03, bty = "n")

#----ADDING NOISE----


####CHOOSE THE UNCERTAINTY LEVEL####
sd_conso=0.5
sd_sources=0.2
###
####ADD NOISE TO CONSUMER SIGNATURE
noise_consumer_C=data.frame(row.names = NULL)
noise_consumer_N=data.frame(row.names = NULL)

nb_obs=10
for (i in 1:nrow(samp_C)){
  noise_consumer_C=rbind(noise_consumer_C,data.frame(Time=rep(samp_C$time[i],nb_obs)
                                                     ,Conso=rnorm(n=nb_obs,mean=samp_C$X[i],sd=sd_conso)))
  noise_consumer_N=rbind(noise_consumer_N,data.frame(Time=rep(samp_N$time[i],nb_obs)
                                                     ,Conso=rnorm(n=nb_obs,mean=samp_N$X[i],sd=sd_conso)))
}


##ADD NOISE TO SOURCE SIGNATURE
noise_source=data.frame(row.names=NULL)
nb_obs_sources=10
for (i in 1: nrow(samp_source)){
  l=data.frame(rep(NA,nb_obs_sources))
  for (j in 2:ncol(samp_source)){
    l=cbind(l,rnorm(n=nb_obs_sources,mean=samp_source[[i,j]],sd=sd_sources))}
  noise_source=rbind(noise_source,data.frame(Time=rep(samp_source[[i,1]],nb_obs_sources),l[-1]))
}

colnames(noise_source)<- c("Time","s1i1","s2i1","s3i1","s1i2","s2i2","s3i2")


##---- source POLYGON----
color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
             rgb(0,0.6,0.45),rgb(0.9,0.6,0),rgb(0.8,0.47,0.65),
             rgb(0.33,0.7,0.9),rgb(0.94,0.89,0.25))
color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),
             rgb(0.9,0.6,0,0.2),rgb(0.8,0.47,0.65,0.2),rgb(0.33,0.7,0.9,0.2),rgb(0.94,0.89,0.25,0.2) )


for(i in 1:nrow(samp_time)){
   # i=1
  t=samp_time[i,]
  #CONSO
  subs_conso_C=subset(noise_consumer_C,noise_consumer_C$Time==t)
  subs_conso_N=subset(noise_consumer_N,noise_consumer_N$Time==t)
  #source
  
  subs_sources=subset(noise_source,noise_source$Time==t)
  moy_sources=data.frame(subs_sources$Time[1])
  sdev_sources=data.frame(subs_sources$Time[1])
  for (j in 2:ncol(noise_source)){
    moy_sources=cbind(moy_sources,mean(subs_sources[[j]]))
    sdev_sources=cbind(sdev_sources,sd(subs_sources[[j]]))}
  
  
  
  title_graph=paste0("t_",t,"d_l_",lambda_data,"large_distrib",'.tiff')  
  tiff(title_graph,  width =4 , height = 5, units = "in", res = 300)
  plot(NULL,main=paste("t=",t,"d"),xlim=c(-10,15),ylim=c(-10,15),xlab="",ylab="",cex.lab=1.4,font.lab=2,las=1) #plot frame
  # mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,cex=1.5,font=2)
  # mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=3,cex=1.5,font=2)
  
  
  ###4 here is the number 
  for (j in 2:4){
    # j=2
    arrows(x0=moy_sources[[j]],y0=moy_sources[[j+3]]-sdev_sources[[j+3]]
           ,x1=moy_sources[[j]],y1=moy_sources[[j+3]]+sdev_sources[[j+3]]
           ,length=0,col=color_main[j-1],lwd=3)
    arrows(x0=moy_sources[[j]]-sdev_sources[[j]],y0=moy_sources[[j+3]]
           ,x1=moy_sources[[j]]+sdev_sources[[j]],y1=moy_sources[[j+3]]
           ,length=0,col=color_main[j-1],lwd=3)
  }
  points(subs_conso_C$Conso,subs_conso_N$Conso,pch=16,cex=1.5)
  # legend("bottomright",legend=c("Zos","grass","Ulva","Entero"),cex=1,col=c(color_main[1:3],"black")
  #        ,lty=c(1,1,1,NA),box.lty=0,pch=c(NA,NA,NA,16))

  dev.off()
}

# # windows()
# plot(NULL,main="",xlim=c(-10,15),ylim=c(-10,15),xlab="",ylab="",cex.lab=1.4,font.lab=2,las=1) #plot frame
# mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=-1,padj=1.5,cex=1.5,font=2)
# mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=-2,cex=1.5,font=2)
# legend("center",legend=c("Source1","Source2","Source3","Consumer"),cex=1.5,col=c(color_main[1:3],"black")
#        ,lty=c(1,1,1,NA),box.lty=0,pch=c(NA,NA,NA,16))
# 
# 

##legend

# tiff("legend_poly_simu.tiff",  width =4 , height = 5, units = "in", res = 300)
# plot(NULL,main="",xlim=c(-10,15),ylim=c(-10,15),xlab="",ylab="",cex.lab=1.4,font.lab=2,las=1) #plot frame
# mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=-1,padj=1.5,cex=1.2,font=2)
# mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=-2,cex=1.2,font=2)
# legend("center",legend=c("Source1","Source2","Source3","Consumer"),cex=1,col=c(color_main[1:3],"black")
#        ,lty=c(1,1,1,NA),box.lty=0,pch=c(NA,NA,NA,16))
# dev.off()

#----DMM----
########PREPARING THE DATA FOR DMM
##CONSUMER
Conso_C_DMM=data.frame(row.names=NULL)
Conso_N_DMM=data.frame(row.names=NULL)



for (i in 1:nrow(samp_time)){
  t=samp_time[i,]
  subs_C=subset(noise_consumer_C,noise_consumer_C$Time==t)
  Conso_C_DMM=rbind(Conso_C_DMM,data.frame(Dates=t,Iso_C=mean(subs_C$Conso),lambda=lambda_data))
  
  subs_N=subset(noise_consumer_N,noise_consumer_N$Time==t)
  Conso_N_DMM=rbind(Conso_N_DMM,data.frame(Dates=t,Iso_N=mean(subs_N$Conso),lambda=lambda_data))
}

list_conso_DMM<-list(Conso_C_DMM,Conso_N_DMM)
liste_lambda_DMM=list(data.frame(Date=samp_time,l=rep(lambda_data,nrow(samp_time))),data.frame(Date=samp_time,l=rep(lambda_data,nrow(samp_time))))

##source
Source_DMM=data.frame(row.names=NULL)
for (i in 1:nrow(samp_time)){
  t=samp_time[i,]
  subs=subset(noise_source,noise_source$Time==t)
  l=data.frame(c(NA))
  for (j in 2:ncol(noise_source)){l=cbind(l,mean(subs[[j]]))}
  Source_DMM=rbind(Source_DMM,data.frame(Dates=t,l[-1]))
}
colnames(Source_DMM)<- c("Date","s1","s2","s3","s1","s2","s3")
list_sources_DMM<-list(Source_DMM[,1:4],Source_DMM[,c(1,5:7)])
sources_approxfun_bis=list(c(approxfun(list_sources_DMM[[1]]$Date,list_sources_DMM[[1]]$s1,method='linear',rule=2)
                         ,approxfun(list_sources_DMM[[1]]$Date,list_sources_DMM[[1]]$s2,method='linear',rule=2)
                         ,approxfun(list_sources_DMM[[1]]$Date,list_sources_DMM[[1]]$s3,method='linear',rule=2)),
                       c(approxfun(list_sources_DMM[[2]]$Date,list_sources_DMM[[2]]$s1,method='linear',rule=2)
                         ,approxfun(list_sources_DMM[[2]]$Date,list_sources_DMM[[2]]$s2,method='linear',rule=2)
                         ,approxfun(list_sources_DMM[[2]]$Date,list_sources_DMM[[2]]$s3,method='linear',rule=2)))




##RUNNING THE DMM

#concentration
conc_DMM=list(data.frame("s1"=rep(1,nb_obs_sources),"s2"=rep(1,nb_obs_sources),"s3"=rep(1,nb_obs_sources))
          ,data.frame("s1"=rep(1,nb_obs_sources),"s2"=rep(1,nb_obs_sources),"s3"=rep(1,nb_obs_sources))) #the conc=1 for all sources and isotope


res_all=all_results(sources_list=list_sources_DMM,data_feat=data_features(list_sources_DMM)
                    ,consu_list=list_conso_DMM,conc_list=conc_DMM,iter=iter
                    ,list_TEF=TEF,X=X,sources_approx=sources_approxfun_bis,
                    lambda_list=liste_lambda_DMM)
DMM=res_all[[3]]

##----GRAPHS DMM----
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
                               Q3V3=quantile(subs$Var3,probs =c(0.75))
                               # ,
                               # mV4=mean(subs$Var4),Q1V4=quantile(subs$Var4,probs =c(0.25)),
                               # Q3V4=quantile(subs$Var4,probs =c(0.75))
                               ))
}



color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
             rgb(0,0.6,0.45),rgb(0.9,0.6,0),rgb(0.8,0.47,0.65),
             rgb(0.33,0.7,0.9),rgb(0.94,0.89,0.25))
color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),
             rgb(0.9,0.6,0,0.2),rgb(0.8,0.47,0.65,0.2),rgb(0.33,0.7,0.9,0.2),rgb(0.94,0.89,0.25,0.2) )
# dev.off()

title_graph=paste0("Contrib_DMM_l_",lambda_data,'.tiff')  
tiff(title_graph, , width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(NULL,xlim=c(0,tail(stats$time,1)),ylim=c(0,1),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)

# # plot(NULL,xlim=c(0,period+2),ylim=c(0,2),xlab="Time (d)",ylab="Contribution")
# points(t,rep(0,period+2),type="l",lwd=1.4,col="black",lty=2)
# points(t,combis2,type="l",lwd=1.4,col="red")
# points(t,combis3,type="l",lwd=1.4,col="green")



polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V1,rev(stats$Q3V1))
        ,col=color_poly[1],border=NA)
polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V2,rev(stats$Q3V2))
        ,col=color_poly[2],border=NA)
polygon(x=c(stats$time,rev(stats$time))
        ,y=c(stats$Q1V3,rev(stats$Q3V3))
        ,col=color_poly[3],border=NA)
# polygon(x=c(stats$time,rev(stats$time))
#         ,y=c(stats$Q1V4,rev(stats$Q3V4))
#         ,col=color_poly[4],border=NA)

points(stats$time,stats$mV1,type="l",col=color_main[1],pch=16,lwd=2)
points(stats$time,stats$mV2,type="l",col=color_main[2],pch=16,lwd=2)
points(stats$time,stats$mV3,type="l",col=color_main[3],pch=16,lwd=2)
# points(stats$time,stats$mV4,type="l",col=color_main[4],pch=16,lwd=2)
dev.off()




#carbon
titleC=paste0("C_l_",lambda_data)
traj_graph(sources_list=list_sources_DMM,data_feat=data_features(list_sources_DMM)
           ,sources_approx=sources_approxfun_bis,consu_list=list_conso_DMM,iso_num=1
           ,lambda_list=liste_lambda_DMM,conc_list=conc_DMM
           ,list_TEF=TEF,iter=iter,X=X,time_range=c(0,tail(temps,1)),range_y=c(-8,8),title=titleC,y_label=expression(paste(" δ"^{13},"C (‰)")))
#nitrogen
titleN=paste0("N_l_",lambda_data)
traj_graph(sources_list=list_sources_DMM,data_feat=data_features(list_sources_DMM)
           ,sources_approx=sources_approxfun_bis,consu_list=list_conso_DMM,iso_num=2
           ,lambda_list=liste_lambda_DMM,conc_list=conc_DMM
           ,list_TEF=TEF,iter=iter,X=X,time_range=c(0,tail(temps,1)),range_y=c(-2,15),title=titleN,y_label=expression(paste(" δ"^{15},"N (‰)")))

#----BDMM----
#prep lambda
mean_lambda=lambda_data
alpha_lambda_data=mean_lambda^2/sd_lambda_data^2
beta_lambda_data=mean_lambda/sd_lambda_data^2

#prep sources
moy_sources_fin=data.frame(row.names = NULL)
sdev_sources_fin=data.frame(row.names = NULL)
for(i in 1:nrow(samp_time)){
  t=samp_time[i,]
  subs_sources=subset(noise_source,noise_source$Time==t)
  moy_sources=data.frame(subs_sources$Time)
  sdev_sources=data.frame(subs_sources$Time)
  for (j in 2:ncol(noise_source)){
    moy_sources=cbind(moy_sources,mean(subs_sources[[j]]))
    sdev_sources=cbind(sdev_sources,sd(subs_sources[[j]]))}
  moy_sources_fin=rbind(moy_sources_fin,moy_sources)
  sdev_sources_fin=rbind(sdev_sources_fin,sdev_sources)
  
}
##
mean_sources_C=moy_sources_fin[c(2:4)]
mean_sources_N=moy_sources_fin[c(5:7)]
sd_sources_C=sdev_sources_fin[c(2:4)]
sd_sources_N=sdev_sources_fin[c(5:7)]
##



####running the model####
library(dplyr)
library(tidyr)
library(ggplot2)
library(cmdstanr)
options(mc.cores = parallel::detectCores())

t_obs <- sort(jitter(noise_consumer_C$Time+0.01 / 1, amount = 0.01))
y_obs <- cbind(noise_consumer_C$Conso, noise_consumer_N$Conso) # ADD NOISE HERE
t_fit <- t_obs
bdmm_data <- list(
  n_obs = length(t_obs),
  n_fit = length(t_fit),
  K = ncol(mean_sources_C),
  D = 3,
  y_obs = y_obs,
  t_obs = t_obs,
  t_fit = t_fit,
  t0 = samp_time$samp_3[1],
  mu_s1 = mean_sources_C,
  mu_s2 = mean_sources_N,
  sigma_s1 = sd_sources_C,
  sigma_s2 = sd_sources_N,
  alpha_lambda = as.numeric(alpha_lambda_data), 
  beta_lambda = as.numeric(beta_lambda_data),
  q1 = rep(1,ncol(mean_sources_C)), 
  q2 = rep(1,ncol(mean_sources_C)), 
  t_mean = mean(t_obs),
  t_sd =  sd(t_obs)
)

bdmm <- cmdstan_model("bdmm_v12_Emi_version.stan", force_recompile = TRUE)

bdmm_fitted <- bdmm$sample(
  data = bdmm_data,
  seed = 102,
  chains = 4,
  parallel_chains = 4,
  refresh = 100
)
  

# Now generate predictions
bdmm_generated <- bdmm$generate_quantities(
  fitted_params = bdmm_fitted,
  data = bdmm_data,
  seed = 999
)

df <- data.frame(
  t_obs = bdmm_data$t_obs, 
  y_obs1 = bdmm_data$y_obs[,1],
  y_obs2 = bdmm_data$y_obs[,2]
)


bdmm_generated_summary <- bdmm_generated$summary() |>
  dplyr::filter(variable |> stringr::str_detect("^y_fit")) %>% 
  tidyr::separate(variable, c("observation", "isotope"), sep = ",") %>% 
  tidyr::pivot_wider(names_from = "isotope", values_from = c("mean", "median", "sd", "mad", "q5", "q95"))
bdmm_generated_summary$time <- bdmm_data$t_fit

# Save the summary to a text file
generated_title=paste0("generated_summary_simu_l_",lambda_data,"large_distrib",".csv")
write.csv(bdmm_generated_summary, file = generated_title, row.names = FALSE)

##FITTING GRAPHS####
# ggplot(bdmm_generated_summary) + 
#   geom_ribbon(aes(time, ymin = `q5_1]`, ymax = `q95_1]`), fill = "grey70") +
#   geom_line(aes(time, `mean_1]`)) + 
#   geom_point(aes(t_obs, y_obs1), data = df, size = 2) + 
#   labs(x = "time", y = "isotope")
# ggsave('v10_plot1.pdf')
# 
# ggplot(bdmm_generated_summary) + 
#   geom_ribbon(aes(time, ymin = `q5_2]`, ymax = `q95_2]`), fill = "grey70") +
#   geom_line(aes(time, `mean_2]`)) + 
#   geom_point(aes(t_obs, y_obs2), data = df, size = 2) + 
#   labs(x = "time", y = "isotope")
# ggsave('v10_plot2.pdf')

#carbon
# dev.off()
title_graph=paste0("fit_BDMM_C_l_",lambda_data,"large_distrib",'.tiff')  
tiff(title_graph, , width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,tail(bdmm_generated_summary$time,1)),ylim=c(-8,8),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(bdmm_generated_summary$time,rev(bdmm_generated_summary$time))
        ,y=c(bdmm_generated_summary$`q5_1]`,rev(bdmm_generated_summary$`q95_1]`))
        ,col=rgb(0.5,0.5,0.5,0.6),border=NA)
points(bdmm_generated_summary$time,bdmm_generated_summary$`mean_1]`,pch=16,col="grey40",cex=1.5)
points(noise_consumer_C$Time,noise_consumer_C$Conso,pch=16,col="red")
dev.off()
#nitrogen
title_graph=paste0("fit_BDMM_N_l_",lambda_data,"large_distrib",'.tiff')  
tiff(title_graph, , width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,tail(bdmm_generated_summary$time,1)),ylim=c(-2,15),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(bdmm_generated_summary$time,rev(bdmm_generated_summary$time))
        ,y=c(bdmm_generated_summary$`q5_2]`,rev(bdmm_generated_summary$`q95_2]`))
        ,col=rgb(0.5,0.5,0.5,0.6),border=NA)
points(bdmm_generated_summary$time,bdmm_generated_summary$`mean_2]`,pch=16,col="grey40",cex=1.5)
points(noise_consumer_N$Time,noise_consumer_N$Conso,pch=16,col="red")
dev.off()


##CONTRIB #####

p_1 <- bdmm_generated$summary() |>
  # Use a regular expression to filter for ones which start with "p_fit" and end in "1]"
  dplyr::filter(variable |> stringr::str_detect("^p_fit.*1\\]$"))
p_1 $time <- bdmm_data$t_fit

p_2 <- bdmm_generated$summary() |>
  # Use a regular expression to filter for ones which start with "p_fit" and end in "1]"
  dplyr::filter(variable |> stringr::str_detect("^p_fit.*2\\]$"))

p_3 <- bdmm_generated$summary() |>
  # Use a regular expression to filter for ones which start with "p_fit" and end in "1]"
  dplyr::filter(variable |> stringr::str_detect("^p_fit.*3\\]$"))

color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
             rgb(0,0.6,0.45),rgb(0.9,0.6,0),rgb(0.8,0.47,0.65),
             rgb(0.33,0.7,0.9),rgb(0.94,0.89,0.25))
color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),
             rgb(0.9,0.6,0,0.2),rgb(0.8,0.47,0.65,0.2),rgb(0.33,0.7,0.9,0.2),rgb(0.94,0.89,0.25,0.2) )

title_graph=paste0("Contrib_BDMM_l_",lambda_data,"_large_distrib",'.tiff')  
tiff(title_graph, , width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(NULL,xlim=c(0,tail(p_1$time,1)),ylim=c(0,1),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)


# points(t,rep(0,period+2),type="l",lwd=1.4,col=color_main[1])
# points(t,combis2,type="l",lwd=1.4,col=color_main[2])
# points(t,combis3,type="l",lwd=1.4,col=color_main[3])




polygon(x=c(p_1$time,rev(p_1$time))
        ,y=c(p_1$q5,rev(p_1$q95))
        ,col=color_poly[1],border=NA)
polygon(x=c(p_1$time,rev(p_1$time))
        ,y=c(p_2$q5,rev(p_2$q95))
        ,col=color_poly[2],border=NA)
polygon(x=c(p_1$time,rev(p_1$time))
        ,y=c(p_3$q5,rev(p_3$q95))
        ,col=color_poly[3],border=NA)
# polygon(x=c(p_1$time,rev(p_1$time))
#         ,y=c(p_4$q5,rev(p_4$q95))
#         ,col=color_poly[4],border=NA)

points(p_1$time,p_1$mean,type="l",col=color_main[1],pch=16,lwd=2)
points(p_1$time,p_2$mean,type="l",col=color_main[2],pch=16,lwd=2)
points(p_1$time,p_3$mean,type="l",col=color_main[3],pch=16,lwd=2)
# points(p_1$time,p_4$mean,type="l",col=color_main[4],pch=16,lwd=2)
# legend("topright",legend=c("MPBxU","POM","TOM"),cex=1,col=c(color_main[1:3]),lty="solid",box.lty=0)
dev.off()
#----Lambda----

# ##enclave
# lambda=0.02
# sd_lambda_syst=0.001
# lambda_data=0.2
# 
# #prep lambda --> entree BDMM
# mean_lambda=lambda_data
# sd_lambda_data=0.1
# alpha_lambda_data=mean_lambda^2/sd_lambda_data^2
# beta_lambda_data=mean_lambda/sd_lambda_data^2
# 
# 
# prior<- rgamma(n=1000,shape=alpha_lambda_data,rate=beta_lambda_data)
# dens_prior<- density(prior)
# 
# #lambda system
# alpha_real=lambda^2/sd_lambda_syst^2
# beta_real=lambda/sd_lambda_syst^2
# 
# syst= rgamma(n=1000,shape=alpha_real,rate=beta_real)
# dens_syst<-density(syst)
# 
# cols<-c("orange","skyblue","black")
# y_max=max(dens_syst[[2]],dens_prior[[2]])
# 
# plot(dens_syst, col = cols[3], xlab="",ylab="",main="",xlim=c(0,0.05)
#      ,lwd=3, cex.lab=1.3,las=1,ylim=c(0,y_max))
# mtext("Density",side=2,line=4,padj=1.5,cex=1.3,font=1)
# mtext("λ (d-1)",side=1,line=3,cex=1.3,font=1)
# points(dens_prior,type="l",col=cols[1],lwd=3)

####




#posterior
lambda_fitted=data.frame(bdmm_fitted$draws(c("lambda")))
lambda_title=paste0("lambda_fitted_simu_l_",lambda_data,"_large_distrib",".csv")
write.csv2(lambda_fitted, file = lambda_title, row.names = FALSE)

#prior
prior<- rgamma(n=1000,shape=alpha_lambda_data,rate=beta_lambda_data)
dens_prior<- density(prior)

#lambda system
alpha_real=lambda^2/sd_lambda_syst^2
beta_real=lambda/sd_lambda_syst^2

syst= rgamma(n=1000,shape=alpha_real,rate=beta_real)
dens_syst<-density(syst)

cols<-c("orange","skyblue","black")

for (i in 1:ncol(lambda_fitted)){
  # i=1
  title_graph<-paste0("dens_graph_l_",lambda_data,"_large_distrib","_",i,".tiff")
  dens_post<-density(lambda_fitted[[i]])
  # y_max=max(dens_syst[[2]],dens_prior[[2]],dens_post[[2]])
  y_max=max(dens_syst[[2]],dens_post[[2]])
  
  tiff(title_graph,  width =4 , height = 5, units = "in", res = 300)
  plot(dens_syst, col = cols[3], xlab="",ylab="",main="",xlim=c(0,0.2)
       ,lwd=3, cex.lab=1.3,las=1,ylim=c(0,y_max)) 
  mtext("Density",side=2,line=4,padj=1.5,cex=1.3,font=1)
  mtext("λ (d-1)",side=1,line=3,cex=1.3,font=1)
  points(dens_prior,type="l",col=cols[1],lwd=3)
  points(dens_post,type="l",col=cols[2],lwd=3)
  dev.off()

}

#----BIAS FOR THE 2 MODELS----

#extracting real_p_values
b=0.025#switch freq
a=0.5 #amplitude
k=0.5 #oscillation axis

#instant real diet
real_p=data.frame(Time=samp_time$samp_3,
                  p_1=rep(0,nrow(samp_time))
                  ,p_2=pB_func(samp_time$samp_3,b,a,k)
                  ,p_3=(1-pB_func(samp_time$samp_3,b,a,k)))

#mean real diet
mean_real_diet_combi2=NULL
for (i in 1:(nrow(samp_time)-1)){
  mean_real_diet_combi2=c(mean_real_diet_combi2,mean(combis2[samp_time[i,]:samp_time[i+1,]]))
}
#data fram of the real mean diet between the sampling times
real_mean_diet=data.frame(time=samp_time,p_1=rep(0,nrow(samp_time)),p_2=c(NA,mean_real_diet_combi2),
                          p_3=c(NA,1-mean_real_diet_combi2))


BDMM_reassembled=data.frame(time=round(p_1$time),Var1=p_1$mean,Var2=p_2$mean,Var3=p_3$mean)

##Estimating bias for each time


stats_MM=data.frame(row.names = NULL)
for (j in 2:nrow(samp_time)){
       # j=2
  times=samp_time[j,]
  subs_real_p=subset(real_p,real_p$Time==times) #instant diet
  subs_mean_p=subset(real_mean_diet,real_mean_diet$samp_3==times) #mean diet
  subs_DMM=subset(DMM,DMM$time==times)
  subs_BDMM=subset(BDMM_reassembled,BDMM_reassembled$time==times)[1,]
  
  lam_T=(samp_time[j,]-samp_time[j-1,])*lambda
  #DMM on mean diet
  DMM_mean=data.frame(m_p1=mean(subs_DMM$Var1),m_p2=mean(subs_DMM$Var2),m_p3=mean(subs_DMM$Var3))
  bias_DMM=(abs(subs_mean_p$p_1-DMM_mean$m_p1)+abs(subs_mean_p$p_2-DMM_mean$m_p2)+abs(subs_mean_p$p_3-DMM_mean$m_p3))*50
  
  #DMM on instant diet
  bias_instant_DMM=(abs(subs_real_p$p_1-DMM_mean$m_p1)+abs(subs_real_p$p_2-DMM_mean$m_p2)+abs(subs_real_p$p_3-DMM_mean$m_p3))*50
  
  #BDMM
  bias_BDMM=(abs(subs_real_p$p_1-subs_BDMM$Var1)+abs(subs_real_p$p_2-subs_BDMM$Var2)+abs(subs_real_p$p_3-subs_BDMM$Var3))*50
  stats_MM=rbind(stats_MM,data.frame(time=times,L_T=lam_T,
                                     MB_DMM=bias_DMM,MB_instant_DMM=bias_instant_DMM,MB_BDMM=bias_BDMM))
  
}

title_bias=paste0("Bias_lam_",lambda_data,"_large_distrib.csv")
write.csv2(stats_MM,file=title_bias,row.names = F)
##GRAPHS####

# 
# # # windows()
# dev.off() # ferme les anciennes fenêtes
# par(fig=c(0,0.5,0,1) )
# #time plot
# plot(NULL,main="",xlim=c(0,samp_time[nrow(samp_time),]),ylim=c(0,100),xlab="Time(d)",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
# points(stats_MM$time,stats_MM$MB_BDMM,type="p",col="black",pch=16)
# points(stats_MM$time,stats_MM$MB_DMM,type="p",col="#E69F00",pch=16)
# points(stats_MM$time,stats_MM$MB_instant_DMM,type="p",col="#D55E00",pch=16)
# #lambdaT plot
# par(fig=c(0.5,1,0,1), new=TRUE)
# plot(NULL,main="",xlim=c(0,5),ylim=c(0,100),xlab="λT",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
# points(stats_MM$L_T,stats_MM$MB_BDMM,type="p",col="black",pch=16)
# points(stats_MM$L_T,stats_MM$MB_DMM,type="p",col="#E69F00",pch=16)
# points(stats_MM$L_T,stats_MM$MB_instant_DMM,type="p",col="#D55E00",pch=16)
# 
# # plot(NULL,type = 'n', axes = F,xlab = '', ylab = '', main = '')
# # legend("center",legend=c("BDMM","DMM mean","DMM instant"),cex=1,col=c("black","#E69F00","#D55E00" )
# #        ,box.lty=0,pch=16)

colors<- c("orange","skyblue","black")
tiff("legend_density_simu.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",legend=c("prior λ_param","posterior λ_param","distrib λ_system "),cex=1,fill=colors,box.lty=0)
dev.off()

#----COMPARING BIAS incertitude effect----

weak=read.csv2("Bias_faible_inc.csv")
weak$fact<- rep(1,nrow(weak))
moyen=read.csv2("Bias_moyen_inc.csv")
moyen$fact<- rep(2,nrow(moyen))
strong=read.csv2("Bias_fort_inc.csv")
strong$fact<- rep(3,nrow(strong))



#BDMM
# dev.off() # ferme les anciennes fenêtes
# par(fig=c(0,0.5,0,1) )
# plot(NULL,main="",xlim=c(0,samp_time[nrow(samp_time),]),ylim=c(0,100),xlab="Time(d)",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
# points(weak$time,weak$MB_BDMM,type="p",col="grey80",pch=16)
# points(moyen$time,moyen$MB_BDMM,type="p",col="grey50",pch=16)
# points(strong$time,strong$MB_BDMM,type="p",col="grey1",pch=16)
# 
# # par(fig=c(0.5,1,0,1), new=TRUE)
# plot(NULL,main="",xlim=c(0,3.1),ylim=c(0,100),xlab="λT",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
# points(weak$L_T,weak$MB_BDMM,type="p",col="grey80",pch=16)
# points(moyen$L_T,moyen$MB_BDMM,type="p",col="grey50",pch=16)
# points(strong$L_T,strong$MB_BDMM,type="p",col="grey1",pch=16)
# 
# plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
# legend("center",legend=c("low","mid","high"),cex=1,col=c("grey80","grey50","grey1"),pch=16,box.lty=0)
# 
# #DMM
# dev.off() # ferme les anciennes fenêtes
# # par(fig=c(0,0.5,0,1) )
# # plot(NULL,main="",xlim=c(0,samp_time[nrow(samp_time),]),ylim=c(0,100),xlab="Time(d)",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
# # points(weak$time,weak$MB_DMM,type="p",col="grey80",pch=16)
# # points(moyen$time,moyen$MB_DMM,type="p",col="grey50",pch=16)
# # points(strong$time,strong$MB_DMM,type="p",col="grey1",pch=16)
# 
# # par(fig=c(0.5,1,0,1), new=TRUE)
# plot(NULL,main="",xlim=c(0,3.1),ylim=c(0,100),xlab="λT",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
# points(weak$L_T,weak$MB_DMM,type="p",col="grey80",pch=16)
# points(moyen$L_T,moyen$MB_DMM,type="p",col="grey50",pch=16)
# points(strong$L_T,strong$MB_DMM,type="p",col="grey1",pch=16)
cols <- c("grey","grey40")
incert_all=rbind(weak,moyen,strong)
incert_all$fact<- as.factor(incert_all$fact)

tiff("varia_effect_simu_DMM.tiff",  width =5 , height = 5, units = "in", res = 300)
boxplot(incert_all$MB_DMM~incert_all$fact,xlab="Variability ",at=c(1,2,3),names=c("Low","Medium","High")
        ,ylab="β (%)",ylim=c(0,100),las=1, cex.lab=1.3,font.lab=1,col=cols[2])
dev.off()
tiff("varia_effect_simu_BDMM.tiff",  width =5 , height = 5, units = "in", res = 300)
boxplot(incert_all$MB_BDMM~incert_all$fact,xlab="Variability ",at=c(1,2,3),names=c("Low","Medium","High")
        ,ylab="β (%)",ylim=c(0,100),las=1, cex.lab=1.3,font.lab=1,col=cols[1])
dev.off()
##stats####
# install.packages("PMCMRplus") 
# library(PMCMRplus)
Fried_BDMM <- friedman.test(MB_BDMM ~ fact|time, data = incert_all)
Fried_DMM <- friedman.test(MB_DMM ~ fact|time, data = incert_all)

#pour DMM
wilcox_test(MB_DMM ~ fact, paired = TRUE,data=incert_all)




#----COMPARING BIAS LAmbda effect----

Bias_0_02=read.csv2("Bias_lam_0_02.csv")
Bias_0_08=read.csv2("Bias_lam_0.08.csv") 
Bias_0_002=read.csv2("Bias_lam_0.002.csv")
Bias_0_006=read.csv2("Bias_lam_0.006.csv")
Bias_0_2=read.csv2("Bias_lam_0.2.csv")

col_bias=c("deepskyblue","blue","black","red","salmon")
#BDMM
tiff("lambda_effect_BDMM.tiff", , width =4 , height = 5, units = "in", res = 300)
plot(NULL,main="",xlim=c(0,3.1),ylim=c(0,100),xlab="λT",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
points(Bias_0_002$L_T,Bias_0_002$MB_BDMM,type="p",col=col_bias[1],pch=16)
points(Bias_0_006$L_T,Bias_0_006$MB_BDMM,type="p",col=col_bias[2],pch=16)
points(Bias_0_02$L_T,Bias_0_02$MB_BDMM,type="p",col=col_bias[3],pch=16)
points(Bias_0_08$L_T,Bias_0_08$MB_BDMM,type="p",col=col_bias[4],pch=16)
points(Bias_0_2$L_T,Bias_0_2$MB_BDMM,type="p",col=col_bias[5],pch=16)
dev.off()
#DMM
tiff("lambda_effect_DMM.tiff",  width =4 , height = 5, units = "in", res = 300)
# plot(NULL,main="",xlim=c(0,3.1),ylim=c(0,100),xlab="λT",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
plot(NULL,main="",xlim=c(0,3.1),ylim=c(0,100),xlab="λT",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
points(Bias_0_002$L_T,Bias_0_002$MB_DMM,type="p",col=col_bias[1],pch=16)
points(Bias_0_006$L_T,Bias_0_006$MB_DMM,type="p",col=col_bias[2],pch=16)
points(Bias_0_02$L_T,Bias_0_02$MB_DMM,type="p",col=col_bias[3],pch=16)
points(Bias_0_08$L_T,Bias_0_08$MB_DMM,type="p",col=col_bias[4],pch=16)
points(Bias_0_2$L_T,Bias_0_2$MB_DMM,type="p",col=col_bias[5],pch=16)
dev.off()

tiff("lambda_effect_legend.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",title="λ(d-1)",legend=c(0.002,0.006,0.02,0.08,0.2),cex=1,col=col_bias,box.lty=0,pch=16)
dev.off()

#----COMPARING LAMBDA EFFECT ALL----

Bias_0_02=read.csv2("Bias_lam_0.02_run2.csv")
Bias_0_08=read.csv2("Bias_lam_0.08_run2.csv")
Bias_0_002=read.csv2("Bias_lam_0.002_run2.csv")
Bias_0_006=read.csv2("Bias_lam_0.006_run2.csv")
Bias_0_2=read.csv2("Bias_lam_0.2_run2.csv")

Bias_all=rbind(Bias_0_002,Bias_0_006,Bias_0_02,Bias_0_08,Bias_0_2)
Bias_all$Lambda<-as.factor(Bias_all$Lambda)

# plot(NULL,main="",xlim=c(0,0.2),ylim=c(0,100),xlab="λ(d-1)",ylab="β(%)",cex.lab=1.4,font.lab=2,las=1)
# points(Bias_all$Lambda,Bias_all$MB_DMM,col="black",type="p",pch=16)
# points(Bias_all$Lambda,Bias_all$MB_BDMM,col="red",type="p",pch=16)


##for GGplot##
#DMM
Bias_GG_DMM<- Bias_all[,c(3,6)]
Bias_GG_DMM$Bias<- Bias_GG_DMM$MB_DMM
Bias_GG_DMM$Model<- rep("DMM",nrow(Bias_all))

#BDMM
Bias_GG_BDMM<- Bias_all[,c(5,6)]
Bias_GG_BDMM$Bias<- Bias_GG_BDMM$MB_BDMM
Bias_GG_BDMM$Model<- rep("BDMM",nrow(Bias_all))

Bias_GG=rbind(Bias_GG_DMM[,-1],Bias_GG_BDMM[,-1])

library(ggplot2)
ggplot(data = Bias_GG, aes(x = Lambda, y = Bias)) + 
  geom_boxplot(aes(fill = Model), width = 0.8) + theme_bw()


cols <- c("grey","grey40")
# boxplot(Bias ~ Lambda + Model, data = Bias_GG,
#         at = c(1:3, 5:7), col = cols,
#         # names = c("", "A", "", "", "B", ""),
#         xaxs = FALSE)
# legend("topleft", fill = cols, legend = c(1,2,3), horiz = T)

tiff("lambda_effect_simu.tiff",  width =5 , height = 5, units = "in", res = 300)
boxplot(Bias_GG$Bias ~ Bias_GG$Model+Bias_GG$Lambda, col = cols,ylim=c(0,100),
        xlab=" λ_param (d-1)",ylab="",at = c(1:2,3:4,5:6,7:8,9:10)
        ,names = c("0.002","",0.006,"","0.2","","0.8","",0.2,""),las=1, cex.lab=1.3,font.lab=1)
mtext("β (%)",side=2,line=4,padj=1.5,cex=1.3,font=1)
dev.off()

tiff("leg_lambda_effect_simu.tiff",  width =5 , height = 5, units = "in", res = 300)
plot(NULL,main="",,type = 'n', axes = F,xlab = '', ylab = '',xlim=c(0,1),ylim=c(0,1))
legend("center",legend=c("BDMM","DMM"),cex=1.2,fill=cols,box.lty=0)
dev.off()


#boxplot DMM
tiff("lambda_effect_simu_DMM_.tiff",  width =5 , height = 5, units = "in", res = 300)
boxplot(Bias_all$MB_DMM~Bias_all$Lambda,xlab=" λ_param (d-1)"
        ,ylab="",ylim=c(0,100),las=1, cex.lab=1.3,font.lab=1,col=cols[2])
mtext("β (%)",side=2,line=4,padj=1.5,cex=1.3,font=1)
points(Bias_all$Lambda,Bias_all$MB_DMM,pch=16,col=rgb(0,0,0,0.5))
dev.off()
#boxplot BDMM
tiff("lambda_effect_simu_BDMM_.tiff",  width =5 , height = 5, units = "in", res = 300)
boxplot(Bias_all$MB_BDMM~Bias_all$Lambda,xlab=" λ_param (d-1)"
        ,ylab="",ylim=c(0,100),las=1, cex.lab=1.3,font.lab=1,col=cols[1])
mtext("β (%)",side=2,line=4,padj=1.5,cex=1.3,font=1)
points(Bias_all$Lambda,Bias_all$MB_BDMM,pch=16)
dev.off()