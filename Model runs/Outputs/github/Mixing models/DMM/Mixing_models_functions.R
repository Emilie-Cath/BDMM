rm(list=ls())# clear the current environment
graphics.off()# clear the current plots


#---------------------------------------------
#              helping functions
#---------------------------------------------

#-----------------data_features --------------------
#function giving information about the studied system

data_features <- function(sources_list) {
  nb_isotopes = length(sources_list)
  nb_sources= ncol(sources_list[[1]])-1
  nb_obs = nrow(sources_list[[1]])
  return(c(nb_isotopes,nb_sources,nb_obs))
}
#returns a vectors with the number of isotopes, number of sources and number of observations in the data set

#---------------Choose_data-----------------------
#sources_list the list of the sources signatures for all isotopes
#consu_list the list of the consumer signatures for all isotopes
#conc_list the list of source concentration for all isotopes
#iter the precision for creating the solutions
#obs_nb the observation number on which to extract data

Choose_data <- function(sources_list,consu_list,conc_list,iter,obs_nb) {
  infos= data_features(sources_list)
  infos_sources=NULL
  infos_consu=NULL
  infos_consu_t1=NULL
  infos_conc= NULL
  combi <- expand.grid(rep(list(seq(0,1,iter)),infos[2])) #data_feat[2] is the number of food sources
  combi <- subset(combi,rowSums(combi)==1) # list of all possible solution of contribution distribution for a precision=iter
  for (i in 1:infos[1]) { #for each isotope (infos[1])
    infos_sources=c(infos_sources,list(data.frame(sources_list[[i]][obs_nb,-c(1)]))) #extracting the sources signatures
    infos_consu=c(infos_consu,list(data.frame(consu_list[[i]][obs_nb,]))) #extracting the consumer signatures
    infos_conc = c(infos_conc,list(data.frame(conc_list[[i]][obs_nb,])))
    infos_consu_t1=c(infos_consu_t1,list(data.frame(consu_list[[i]][(obs_nb+1),])))} #extracting the consumer singature at the next observation
  print(c("combi","sources","consu","concentration","consu next obs"))
  return(list(combi,infos_sources,infos_consu,infos_conc,infos_consu_t1)) 
}

#return a list with the combinations, the sources signatures, the consumer signatures, the concentrations, the consumer signqture at the next observation

#------------Mean_sources--------------------------------
#t the end of the time window in whiwh to average the sources signatures (d)
#lambda_list the list of lambdas for all isotopes
#sources_approx is the list of function that gives the source signatures

averaged_sources<- function(t,lambda_list,sources_approx){
  V=NULL
  for ( i in 1:length(lambda_list)) { #for all isotopes
    lambda=approxfun(lambda_list[[i]][[1]], lambda_list[[i]][[2]],method = 'const',rule=2) #creating a approxfun of lambdas
    t0=t-(2*log(2)/lambda(t-1)) #  creating the window in which to average the sources signatures
    Sources=NULL
    for (j in 1:length(sources_approx[[i]])){ #for each source
      Sources=cbind(Sources,mean(sources_approx[[i]][[j]](seq(t0,t,1))))} #keeping the mean value of each source
    V=c(V,list(Sources))}
  return(V)
}
#returns a vector containing the averaged sources signatures on the studied time window


#---------------------------------------------
#             Mixing functions
#---------------------------------------------

#--------------Mixing--------------
# t time (d) when to proceed 
#sources_approx is the list of function that gives the source signatures for one isotope
#p: one solution (a combination of proportions) for one isotope
#TEF: trophic enrichment factor for each source for one isotope
#Conc concentrations of one isotope in the sources


Mixing <- function(t,sources_approx,p,TEF,Conc){ 
  Sum1<-NULL
  Sum2<-NULL
  for (i in 1:length(sources_approx)){
    Sum1<-cbind(Sum1,p[[i]]*Conc[[i]]*(sources_approx[[i]](t)+TEF[[i]])) 
    Sum2<-cbind(Sum2,p[[i]]*Conc[[i]])}
  return(rowSums(Sum1)/rowSums(Sum2))} 

#--------------Mixing_SMM_delta--------------
# t time (d) when to proceed 
#averaged_sources is a list of the averaged signatures of the sources on the time window
#p: one solution (a combination of proportions) for one isotope
#TEF: trophic enrichment factor for each source for one isotope
#Conc concentrations of one isotope in the sources


Mixing_SMM_delta <- function(averaged_sources,p,TEF,Conc){ 
  Sum1<-NULL
  Sum2<-NULL
  for (i in 1:length(averaged_sources)){
    Sum1<-cbind(Sum1,p[[i]]*Conc[[i]]*(averaged_sources[[i]]+TEF[[i]]))
    Sum2<-cbind(Sum2,p[[i]]*Conc[[i]])}
  return(rowSums(Sum1)/rowSums(Sum2))} 


#---------------------------------------------
#         computing distance functions
#---------------------------------------------

#--------------diff_obs-pred--------------
#Mix is the data frame with the DX_pred for one isotope and the associated contributions
#DX_obs is the observed consumer signature for one isotope
diff_obs_pred<-function(Mix,Dx_obs){
  diff=data.frame(row.names = NULL)
  for (i in 1:nrow(Mix)){
    a=((Dx_obs-Mix[[i,1]])/Dx_obs)^2 #HERE WE DIVIDE BY DX_OBS and we apply the squaredfucntion /!\
    diff=rbind(diff,a)}
  return(diff)}
#returns a data frame with the weighted difference between the consumer predicted and observed value


#--------------distance--------------
#DIFF a data frame and each column is the difference for one isotope (diff iso1|diff iso 2| etc...)



distance <- function(DIFF){
  V=data.frame(row.names = NULL)
  for (i in  1:nrow(DIFF)){
    dist=sqrt(sum(DIFF[i,]))
    V=rbind(V,dist)}
  return(V)
}

#returns a data frame with and unique column of the distance value for each solution

#--------------sort_and_choose--------------

#Dist is a data frame with the distance and the associated contributions
#X is the number of solution that should be kept
sort_and_choose <- function(dist_combi,X){
  dist_ordo=dist_combi[order(dist_combi[[1]],decreasing = TRUE),]
  return(tail(dist_ordo,X))
}
#returns a data frame with the distance and the associated contributions of the X best solutions
#----------------sol_selection----------------------

#data_feat is the output of the function data_features(
#Mix is a data frame with the results of the mixing model for all isotopes
#combi 
#vect_Dx_obs is a vector containing the observerd consumer signatures for all isotopes
#X the number of solutions to be kept

sol_selection <- function(Mix,vect_Dx_obs,data_feat,combi,X) {
  nb_iso=data_feat[1] 
  diff=data.frame(row.names = NULL)
  for (i in 1:nb_iso) { #for each isotope
    diff=c(diff,list(diff_obs_pred(data.frame(Mix[[i]]),vect_Dx_obs[i]))) #computing the squared weighted difference
  }
  dist=distance(data.frame(diff)) #computing the distance
  dist_combi=cbind(dist,combi)
  selected_sol=sort_and_choose(dist_combi,X) #aranging the solutions and keeping the X best
  colnames(selected_sol)[1]=c("Dist")
  return (selected_sol)
}

#computes the distance and the X best solutions out of the mixing model results

#---------------------------------------------
#             SMM,DMM,SMM_delta
#---------------------------------------------


#----------------SMM_func----------------------
#sources_approx  is the list of the functions approximating the sources for each isotope
#here combi is the datafeame containing all the solutions
#data_feat is the output of the data_feature function
#details_conc is the ouput of choose_data function for the concentration
#details_consu is the ouput of choose_data function for the consumer signatures
#t is the time (in days) when the SMM should be applied
#list_TEF is a list of data frame containing the TEF for each isotope
#source approx is the list of functions expressing the sources signatures along time
#X is the number of solutions that should be kept by the model

SMM_func <- function(t,sources_approx,list_TEF,details_conc,combi,data_feat,details_consu,X){
  Dx_pred=NULL
  vect_DX_obs=NULL
  for (i in 1:data_feat[1]){ #for each isotope
    Mix=data.frame(row.names=NULL)
    for (j in 1:nrow(combi)) {  #for each solution combination 
      Mix=rbind(Mix,Mixing(t=t,sources_approx=sources_approx[[i]],p=combi[j,],TEF=list_TEF[[i]],Conc=details_conc[[i]]))} #creating predictions for the consumer signature for each combination
    Dx_pred=c(Dx_pred,list(Mix))
    vect_DX_obs=c(vect_DX_obs,details_consu[[i]][[2]])
  } 
  Mix_res=data.frame(Dx_pred)
  Prop=sol_selection(Mix=Mix_res,vect_Dx_obs=vect_DX_obs,data_feat=data_feat,combi=combi,X=X)
  return(Prop)
}  
#returns the results of the SMM at a precise time

#----------------SMM_delta_func----------------------
#sources_approx  is the list of the functions approximating the sources for each isotope
#ici combi is the datafeame containing all the solutions
#data_feat is the output of the data_feature function
#details_conc is the ouput of choose_data function for the concentration
#details_consu is the ouput of choose_data function for the consumer signatures
#tis the time (in days) when the SMM should be applied
#list_TEF is a list of data frame containing the TEF for each isotope
#source approx is the list of functions expressing the sources signatures along time
#X is the number of solutions that should be kept by the model
#lambda_list is the list of the lammbda value for each isotope
SMM_delta_func<- function(t,sources_approx,list_TEF,details_conc,combi,data_feat,details_consu,X,lambda_list) {
  Dx_pred=NULL
  vect_DX_obs=NULL
  avr_sources=averaged_sources(t=t,lambda_list=lambda_list,sources_approx =sources_approx )
  for (i in 1:data_feat[1]){ #for each isotope
    Mix=data.frame(row.names=NULL)
    for (j in 1:nrow(combi)) {#pour each solution combination 
      Mix=rbind(Mix,Mixing_SMM_delta(averaged_sources=avr_sources[[i]],p=combi[j,],TEF=list_TEF[[i]],Conc=details_conc[[i]]))}#creating predictions for the consumer signature for each combination
    Dx_pred=c(Dx_pred,list(Mix))
    vect_DX_obs=c(vect_DX_obs,details_consu[[i]][[2]])
  }
  Mix_res=data.frame(Dx_pred)
  Prop=sol_selection(Mix=Mix_res,vect_Dx_obs=vect_DX_obs,data_feat=data_feat,combi=combi,X=X)
  return( Prop)
}
#returns the results of the integrated SMM at a precise time

#----------------DMM_func----------------------
#sources_approx  is the list of the functions approximating the sources for each isotope
#ici combi is the datafeame containing all the solutions
#data_feat is the output of the data_feature function
#details_conc is the ouput of choose_data function for the concentration
#details_consu is the ouput of choose_data function for the consumer signatures at the beginning of the period
#tis the time (in days) when the SMM should be applied
#list_TEF is a list of data frame containing the TEF for each isotope
#source approx is the list of functions expressing the sources signatures along time
#X is the number of solutions that should be kept by the model
#lambda_list is the list of the lammbda value for each isotope
#details_consu_1 is the output of choose_data function for the consumer signatures at the end of the period

DMM_func<- function(details_consu,details_consu_t1,sources_approx,list_TEF,details_conc,Combi,X,lambda_list,data_feat){
  Dx_pred=NULL
  vect_DX_obs=NULL
  for (j in 1:length(sources_approx)){#For each isotope
    Mix=NULL
    for (k in 1:nrow(Combi)) { #on fait tourner le solver pour toutes les combi rt on reccupe la valeur prédite à t+1
      Mix<-rbind(Mix,tail(RUN_TIMdyn2(t=c(details_consu[[j]][[1]],details_consu_t1[[j]][[1]])
                                              ,state0=c(X=details_consu[[j]][[2]]),
                                              par=c(lambda = approxfun(lambda_list[[j]][[1]], lambda_list[[j]][[2]],method = 'const',rule=2), 
                                                    Xinf=approxfun(c(details_consu[[j]][[1]],details_consu_t1[[j]][[1]]),
                                                                   Mixing(t=c(details_consu[[j]][[1]],details_consu_t1[[j]][[1]]),
                                                                          sources_approx=sources_approx[[j]],p=Combi[k,],TEF=list_TEF[[j]]
                                                                          ,Conc=details_conc[[j]]),rule=2))),1))}
    Dx_pred=c(Dx_pred,list(Mix[[2]]))
    vect_DX_obs=c(vect_DX_obs,details_consu_t1[[j]][[2]])
    }
  Mix_res=data.frame(Dx_pred)
  Prop=sol_selection(Mix=Mix_res,vect_Dx_obs=vect_DX_obs,data_feat=data_feat,combi=Combi,X=X)
  return( Prop)
}
#returns the results of the DMM on a time period

#------------ODE solving functions (using the Desolve package)-------------------------------  




#numerical solution of ODE
#with state0 the initial conditions 
#and par composed of lambda and Xinf forcing functions 
RUN_TIMdyn2<-function(t,state0,parms){ 
  out<-as.data.frame(lsoda(state0,t,ODE_TIMdyn2,parms)) # lsoda is an integrator to solve ODE
  return(out)}

# differential equation of the dynamic mixing model
ODE_TIMdyn2<-function(t,state,parameters){ #time incorporation model : ODE
  with(as.list(c(state,parameters)),{  # unpack the state variables, parameters
    #state variable X(t) : the response
    #Forcing variables lambda(t) and Xinf(t)
    dX <- lambda(t)*(Xinf(t)-X)
    # the outputs
    list(c(dX = dX), c(Xinf = Xinf(t))) #the response is saved
  })
}


#------------all_results-------------------------------
#sources_approx  is the list of the functions approximating the sources for each isotope
#data_feat is the output of the data_feature function
#list_TEF is a list of data frame containing the TEF for each isotope
#source approx is the list of functions expressing the sources signatures along time
#X is the number of solutions that should be kept by the model
#lambda_list is the list of the lambda value for each isotope
#source list is list of data frame containing the observed sources signatures for each isotope
#consulist is list of data frame containing the observed consumer signatures for each isotope
#source list is list of data frame containing the concentration for each isotope
#iter is the precision of solution that will be created

all_results <- function(sources_list,data_feat,consu_list,conc_list,iter,list_TEF,X,sources_approx,lambda_list){
  val_SMM=data.frame(row.names = NULL)
  val_SMM_delta=data.frame(row.names=NULL)
  val_DMM=data.frame(row.names = NULL)
  
  #first obs
  A=Choose_data(sources_list=sources_list,consu_list=consu_list,conc_list=conc_list,iter=iter,obs_nb=1)
  SMM=SMM_func(t=A[[3]][[1]][[1]],sources_approx=sources_approx,list_TEF=list_TEF,details_conc=A[[4]]
               ,combi=A[[1]],data_feat=data_feat,details_consu=A[[3]],X=X) #computing the SMM for the first observation
  DMM=DMM_func(details_consu=A[[3]],details_consu_t1=A[[5]],sources_approx=sources_approx,
               list_TEF=list_TEF,details_conc=A[[4]],Combi=A[[1]],X=X,lambda_list=lambda_list,data_feat=data_feat) #computing the DMM for the first time window
  
  val_SMM=rbind(val_SMM,cbind(time=A[[3]][[1]][[1]],SMM)) 
  val_DMM=rbind(val_DMM,cbind(time=A[[5]][[1]][[1]],DMM)) 
  #2nd obs to last obs-1
  
  
  for (i in 2:(data_feat[3]-1)) {
    A=Choose_data(sources_list=sources_list,consu_list=consu_list,conc_list=conc_list,iter=iter,obs_nb=i)
    SMM=SMM_func(t=A[[3]][[1]][[1]],sources_approx=sources_approx,list_TEF=list_TEF,details_conc=A[[4]]
                 ,combi=A[[1]],data_feat=data_feat,details_consu=A[[3]],X=X) #computing the SMM in the middle observations
    DMM=DMM_func(details_consu=A[[3]],details_consu_t1=A[[5]],sources_approx=sources_approx,
                 list_TEF=list_TEF,details_conc=A[[4]],Combi=A[[1]],X=X,lambda_list=lambda_list,data_feat=data_feat)#computing the DMM in the middle observations
    SMM_delta=SMM_delta_func(t=A[[3]][[1]][[1]],sources_approx=sources_approx,list_TEF=list_TEF,details_conc=A[[4]],
                             combi=A[[1]],data_feat=data_feat,details_consu=A[[3]],X=X,lambda_list=lambda_list)#computing the SMM_delta in the middle observations

    window_delta=A[[3]][[1]][[1]]-(2*log(2)/lambda_list[[1]][i-1,2] ) #determining the start of the time window sor the SMM_delta
    val_SMM_delta=rbind(val_SMM_delta,cbind(cbind(t0=window_delta,t1=A[[3]][[1]][[1]]),SMM_delta))
    val_SMM=rbind(val_SMM,cbind(time=A[[3]][[1]][[1]],SMM)) 
    val_DMM=rbind(val_DMM,cbind(time=A[[5]][[1]][[1]],DMM))}
  
  
  #last observation
  A=Choose_data(sources_list=sources_list,consu_list=consu_list,conc_list=conc_list,iter=iter,obs_nb=data_feat[3])
  SMM_delta=SMM_delta_func(t=A[[3]][[1]][[1]],sources_approx=sources_approx,list_TEF=list_TEF,details_conc=A[[4]],
                           combi=A[[1]],data_feat=data_feat,details_consu=A[[3]],X=X,lambda_list=lambda_list)#computing the SMM_delta for the last observation
  SMM=SMM_func(t=A[[3]][[1]][[1]],sources_approx=sources_approx,list_TEF=list_TEF,details_conc=A[[4]]
               ,combi=A[[1]],data_feat=data_feat,details_consu=A[[3]],X=X) #computing the SMM for the last observation
  window_delta=A[[3]][[1]][[1]]-(2*log(2)/lambda_list[[1]][data_feat[3]-1,2] ) #determining the start of the time window sor the SMM_delta
  val_SMM_delta=rbind(val_SMM_delta,cbind(cbind(t0=window_delta,t1=A[[3]][[1]][[1]]),SMM_delta))
  val_SMM=rbind(val_SMM,cbind(time=A[[3]][[1]][[1]],SMM))
  
  return( c(list(val_SMM),list(val_SMM_delta),list(val_DMM)))
  }
#returns the outputs of the three models for all observations

#---------------------------------------------
#             graphs
#---------------------------------------------




#------------traj_graph-------------------------------
# sources_list: is list of data frame containing the observed consumer signatures for each isotope 
# sources_approx: is the list of functions expressing the sources signatures along time
# consu_list:is list of data frame containing the observed consumer signatures for each isotope
# iso_num: specifies the isotope you are working on (usually 1 for Carbon and 2 for Nitrogen)
# lambda_list is the list of data frame containing the lambda value for each isotope
# conc_list is the list of data frames of source concentration for all isotopes
# list_TEF is the list of data frames of TEF for all isotopes
# iter is the precision regarding the generation of solutions
# X is the number of solutions to be kept
# time_range a vector containing the starting and ending time (in days) of the experiment
# range_y a vector containing the isotopic range worked on
# title
# y_label


traj_graph <- function(sources_list,data_feat,sources_approx,consu_list,iso_num,lambda_list,conc_list,list_TEF,iter,X,time_range,range_y,title,y_label) {
  title_file=paste0("fit_DMM_",title,".tiff")
  tiff(title_file,  width =4 , height = 5, units = "in", res = 300)
  plot(NULL,xlim=time_range,ylim=range_y,xlab="Time (d)",ylab="",main=title,cex.lab=1.3,las=1,cex.axis=1.1)
  mtext(y_label,side=2,line=4,padj=1.5,cex=1.3,font=2)
  for (j in 1:(data_feat[3]-1)){ #pour toutes les observations -1 
    A=Choose_data(sources_list=sources_list,consu_list=consu_list,conc_list=conc_list,iter=iter,obs_nb=j)
    DMM=DMM_func(details_consu=A[[3]],details_consu_t1=A[[5]],sources_approx=sources_approx,
                 list_TEF=list_TEF,details_conc=A[[4]],Combi=A[[1]],X=X,lambda_list=lambda_list,data_feat=data_feat)
    for (k in 1:nrow(DMM)) {mix=RUN_TIMdyn2(t=seq(A[[3]][[iso_num]][[1]],A[[5]][[iso_num]][[1]],1)
                                            ,state0=c(X=A[[3]][[iso_num]][[2]]),
                                            par=c(lambda = approxfun(lambda_list[[iso_num]][[1]], lambda_list[[iso_num]][[2]],method = 'const',rule=2), 
                                                  Xinf=approxfun(c(A[[3]][[iso_num]][[1]],A[[5]][[iso_num]][[1]]),
                                                                 Mixing(t=c(A[[3]][[iso_num]][[1]],A[[5]][[iso_num]][[1]]),
                                                                        sources_approx=sources_approx[[iso_num]],p=DMM[-1][k,],TEF=list_TEF[[iso_num]]
                                                                        ,Conc=A[[4]][[iso_num]]),rule=2)))[1:2]
    points(mix[[1]],mix[[2]],type="l",col="grey")}
    }
  points(consu_list[[iso_num]][[1]],consu_list[[iso_num]][[2]],type="p",
         col="red",pch=16) #on ajoute les valeurs observées
  dev.off()
    
}





#------------sign_biplot-----------------------------
# sources_list: is list of data frame containing the observed consumer signatures for each isotope 
# consu_list:is list of data frame containing the observed consumer signatures for each isotope
# lambda_list is the list of data frame containing the lambda value for each isotope
# list_TEF is the list of data frames of TEF for all isotopes
#title
#range_C a vector containing the range of the carbon signatures
#range_N a vector containing the range of the nitrogen signatures
#source_names a vector containing the name of all sources

sign_biplot <-function(consu_list,sources_list,list_TEF,title,range_C,range_N,source_names) {
  data_feat=data_features(sources_list)
  set_couleurs_sources=list(c("deepskyblue4","deepskyblue"),c("aquamarine4","aquamarine"),c("olivedrab3","olivedrab1"),c("seagreen3","seagreen1"),c("turquoise3","turquoise1")) #works for up to 5 sources
  set_legende=c("deepskyblue","aquamarine","olivedrab1","seagreen1","turquoise1")
  set_couleurs_conso=c("firebrick3","firebrick1")
  
  
  sources_corr=NULL #this loop returns a list of data frame with the corrected observed sources signatures
  for (k in 1:data_feat[1]){
    B=NULL
    for (j in 2:ncol(sources_list[[k]])){
      A=NULL
      for (i in 1:nrow(sources_list[[k]])) {A=c(A,sources_list[[k]][[i,j]]+list_TEF[[k]][[j-1]])}
      B=c(B,list(A))}
    B=data.frame(B)
    sources_corr=c(sources_corr,list(B))
  }
  
  
  
  s=consu_list[[1]][[1]] #create a color gradient
  colfunc<-colorRampPalette(set_couleurs_conso)
  colors <- (colfunc(data_feat[3]))
  colors1 <- colors[rank(s)]
  windows()
  plot(NULL,main=title,xlim=range_C,ylim=range_N,xlab=expression(paste(" δ"^{13},"C (‰)")),ylab="",cex.lab=1.4,font.lab=2,las=1) #plotting
  mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=2)
  arrows(x0=consu_list[[1]][[2]][-c(data_feat[3])],x1=consu_list[[1]][[2]][-c(1)],y0=consu_list[[2]][[2]][-c(data_feat[3])],
         y1=consu_list[[2]][[2]][-c(1)],lty=1,length=0.2,col=colors1,cex=1.2) 
  points(consu_list[[1]][[2]],consu_list[[2]][[2]],type="p",col=colors1,pch=16,lwd=1.2,cex=1.6) #plotting the consumer signatures
  
  for (i in 1:data_feat[2]){ #for each source
    s=sources_list[[1]][[1]] 
    colfunc<-colorRampPalette(set_couleurs_sources[[i]]) #create a color gradient for the sources
    colors <- (colfunc(data_feat[3]))
    colors1 <- colors[rank(s)]
    arrows(x0=sources_corr[[1]][[i]][-c(data_feat[3])],x1=sources_corr[[1]][[i]][-c(1)],y0=sources_corr[[2]][[i]][-c(data_feat[3])],
           y1=sources_corr[[2]][[i]][-c(1)],lty=1,length=0.1,col=colors1) 
    points(sources_corr[[1]][[i]],sources_corr[[2]][[i]],type="p",col=colors1,pch=16,cex=1.6)} #plotting the sources
  legend("bottomright",legend=c(source_names,"consumer"),cex=1,col=c(set_legende[1:data_feat[2]],set_couleurs_conso[1]),lty="solid",box.lty=0)
}

#returns a biplot graph containing the corrected signatures of the sources and the consumer signatures


#------------sources_contrib-----------------------------
#mod_results is the output of the all_results function for X=1
#obs_time is a vector containing the date of all observations
#source_num indicates on wi=hich source you are working on
#where_legend can be either "bottomleft"...


sources_contrib <- function(mod_results,obs_time,title,source_num,where_legend) {
  SMM=mod_results[[1]]
  SMM_delta=mod_results[[2]]
  DMM=mod_results[[3]]
  my_col= c("deeppink2", "green","royalblue4")
  windows()
  plot(NULL,xlim=c(obs_time[1],obs_time[nrow(SMM)]),ylim=c(0,1),xlab="Time (d)",ylab="Contribution ",main=title,cex.lab=1.4,font.lab=2,las=1) #plot frame
  for (i in 1:nrow(DMM)) {
    lines(c(obs_time[i],DMM[[i,1]]),rep(DMM[[i,source_num+2]],2),type="l",col=my_col[2],lwd=5) #plotting DMM
    lines(c(SMM_delta[[i,1]],SMM_delta[[i,2]]),rep(SMM_delta[[i,source_num+3]],2),type="l",col=my_col[3],lwd=2)} #plotting SMM_delta
  points(SMM[[1]],SMM[[source_num+2]],type="p",col=my_col[1],pch=16,cex=1.4) #plotting SMM
  legend(where_legend,legend=c("SMM","DMM","SMM delta "),col=my_col,pch=c(19,15,15),bty="n",cex=1.1)

  }
#returns for one source the estimation made by the best solution of each model and shows the intergation time

#------------ridgeline_chart------------------------------
#res_model is a data frame with the results of one of the three models for each observation
#model_name should be either "SMM", "SMM_delta" or "DMM"
#source_names a vector containing the name of all sources
#X the number of kept solutions

ridgeline_chart <- function(res_model,model_name,sources_names,X){
  frame_vide= data.frame(time=rep(0,X))
  nom=c("Var1","Var2","Var3","Var4","Var5","Var6","Var7") #goes up to 7 sources
  if (model_name=="SMM_delta"){modele_corr=res_model[-c(1,3)]
  colnames(modele_corr)[1]<-c("time")
  modele_corr$time=as.numeric(modele_corr$time)
  for (i in 2:ncol(modele_corr)) {
    frame_vide<- cbind(frame_vide,rep(NA,X))
    colnames(frame_vide)[i]<-c(nom[i-1])}
  modele_corr=rbind(frame_vide,modele_corr)}
  else if(model_name=="DMM") {modele_corr=res_model[-c(2)]
  for (i in 2:ncol(modele_corr)) {
    frame_vide<- cbind(frame_vide,rep(NA,X))
    colnames(frame_vide)[i]<-c(nom[i-1])}
  modele_corr=rbind(frame_vide,modele_corr)}
  else{modele_corr=res_model[-c(2)]}
  res=data.frame(row.names = NULL)
  nb_sources=ncol(modele_corr)-1
  for (i in 1:nb_sources){
    Val=data.frame(time=modele_corr[[1]],cbind(modele_corr[[i+1]],sources=rep(i,nrow(modele_corr))))
    res=rbind(res,Val)}
  
  
  my_col=c(rgb(0,0.45,0.7,0.8),rgb(0.83,0.36,0,0.8),rgb(0,0.6,0.45,0.8),rgb(0.94,0.89,0.25,0.8),rgb(0.8,0.47,0.65,0.8),rgb(0.9,0.6,0,0.8),rgb(0.33,0.7,0.9)) #works for up to 7 sources
  
  res$sources=as.factor(res$sources)
  res$time=as.factor(res$time)
  
  
  
   windows()
  ggplot(res, aes(x = V1, y = time, fill = time)) +
    geom_density_ridges(aes(fill=sources))+ 
    theme_ridges()+
    scale_fill_manual(values=c(my_col[c(1:nb_sources)]),name = "Sources", labels = sources_names)+
    labs(x="Contribution ",y="Time (d)")+
    ggtitle(paste("Contribution of the sources over time",model_name ))+
    scale_x_continuous(limits=c(0,1))+
    theme(axis.title=element_text(size=16,face="bold"))
  
}

#returns ridgeline graphs of the contribution along time predicted by one of the three models

#------------evol_sources------------------------------
#temps est un vecteur indiquant le temps de l'experience par ex c(0,360)
#etendue est l'etendue des valeurs que prennent les sources (rajouter un peu en dessous pour la legende)
#par ex c(0,11) et c(-40,-15)
#isotope le nom de l'isotope etudie
#sources la liste d'approxfun de la source voulue
#sources_name le vecteur avec le nom des sources

evol_sources <- function(exp_range,iso_range,isotope,sources,sources_name){
  couleur=c("deepskyblue","aquamarine","olivedrab1","seagreen1","turquoise1") #works with up to 5 sources
  windows()
  plot(NULL,xlim=exp_range, ylim=iso_range,xlab="Time (d)", ylab="isotopic ratio (‰)", main=paste("Sources variation",isotope),cex.lab=1.4,font.lab=2)
  for (i in 1:length(sources)) {
    points(seq(exp_range[1],exp_range[2],1),sources[[i]](seq(exp_range[1],exp_range[2],1)),type="l",col=couleur[[i]],lwd=3)
  }
  legend("bottomright",legend=sources_name,cex=1,col=c(couleur[1:3]),lty="solid",box.lty=0)
  
}


#------------------------poly_sources_evol--------------------
# consu_list
# sources_list
# list_TEF
# sources_name
# range_C
# range_N

poly_sources_evol<-function(consu_list,sources_list,list_TEF,sources_name,range_C,range_N){
  set_couleur_transparent=c(rgb(0,1,0,0.2),rgb(0,0.5,0.5,0.2),rgb(0,0,1,0.2),rgb(0.5,0.5,0,0.2),rgb(0,0.7,0.5,0.2)) #marche jusqu'à 5 sources
  set_couleur=c(rgb(0,1,0,1),rgb(0,0.5,0.5,1),rgb(0,0,1,1),rgb(0.5,0.5,0,1),rgb(0,0.7,0.5,1))
  data_feat=data_features(sources_list)
  sources_corr=NULL #cette boucle renvoie une liste de data frame des sources corrigees
  for (k in 1:data_feat[1]){
    B=NULL
    for (j in 2:ncol(sources_list[[k]])){
      A=NULL
      for (i in 1:nrow(sources_list[[k]])) {A=c(A,sources_list[[k]][[i,j]]+list_TEF[[k]][[j-1]])}
      B=c(B,list(A))}
    B=data.frame(B)
    sources_corr=c(sources_corr,list(B))
  }
  
  
  for (obs in (1:data_feat[3])) {
    windows()
    plot(NULL,main=paste("Period",obs),xlim=range_C,ylim=range_N,xlab="",ylab="",cex.lab=1.4,font.lab=2,las=1) #plot frame
    mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.5,font=2)
    mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=3,cex=1.5,font=2)
    points(consu_list[[1]][[2]],consu_list[[2]][[2]],type="p",pch=16,col=rgb(1,0,0,0.2),cex=1.2) #consumer signatures
    for (i in 1:data_feat[2]) { #sources colored in transparent colors
      points(sources_corr[[1]][[i]],sources_corr[[2]][[i]],type="p",pch=16,col=set_couleur_transparent[i],cex=1.2)
    }
    for (j in 1:(data_feat[2]-1)) { #first two segments of the polygon
      arrows(x0=sources_corr[[1]][[obs,j]],x1=sources_corr[[1]][[obs,j+1]]
             ,y0=sources_corr[[2]][[obs,j]],y1=sources_corr[[2]][[obs,j+1]],length=0)
    }
    arrows(x0=sources_corr[[1]][[obs,data_feat[2]]],x1=sources_corr[[1]][[obs,1]],
           y0=sources_corr[[2]][[obs,data_feat[2]]],y1=sources_corr[[2]][[obs,1]],length=0) #last segment
    for (i in 1:data_feat[2]) { #lcolored source points
      points(sources_corr[[1]][[obs,i]],sources_corr[[2]][[obs,i]],type="o",pch=16,col=set_couleur[i],cex=1.6)
    }
    points(consu_list[[1]][[obs,2]],consu_list[[2]][[obs,2]],type="o",pch=16,col=rgb(1,0,0,1),cex=1.6) #colored consumer point
    legend("bottomright",legend=c(sources_name,"consumer"),cex=1.1,col=c(set_couleur[1:data_feat[2]],rgb(1,0,0,1)),box.lty=0,pch=16)
    
  }
}


#------------dist_graph------------------------------
# 

dist_graph <- function(obs_time,res_model,range_y) {
  colnames(res_model[[2]])[2]<-c("time")
  my_col=c("deeppink2", "royalblue4", "green")
  enlargment=c(1.6,1.4,1.2)
  windows()
  plot(NULL,xlim=c(obs_time[1],obs_time[length(obs_time)]),ylim=range_y,xlab="Time (d)", 
       ylab="Distance",cex.lab=1.4,font.lab=2,las=1)
  for ( i in 1:3){
    res=res_model[[i]]
    points(res$time,res$Dist,type="p",pch=16,col=my_col[[i]],cex=enlargment[i])}
  
  legend("topright",legend=c("SMM","SMM_delta","DMM"),cex=1.1,col=my_col,pch=c(19,19,19),box.lty=0)}




#------------------------contrib_time--------------------
#graph of the contributions over time 
#res_model the output of a mixing model
#model_name either "SMM", "SMM_delta","DMM"
#date a vector of the sampling time
#data_feat the output of the data_features function
#X the number of kept solutions in the mixing model
#source_name a vector with the name of the different sources
#where_legend:either "bottomleft","topleft","bottomright","topright"

contrib_time <- function(res_model,model_name,date,data_feat,X,sources_name,where_legend) {
  nb_sources=data_feat[2]
  if (model_name!="SMM") {date=date[-1]}
  if (model_name=="SMM_delta") {
    res_model=res_model[-1]
    colnames(res_model)[1] ="time"}
  aranged_data=data.frame(row.names = NULL)
  for (j in 1:length(date)) {
    sub=subset(res_model,res_model$time==date[j])
    df=data.frame(time=date[j])
    for (i in 1:nb_sources){
      df=cbind(df,sub[X,i+2],min(sub[i+2]),max(sub[i+2]))
    }
    aranged_data=rbind(aranged_data,df)
  }
  ind=seq(from=2,by=3,length=nb_sources)
  color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
               rgb(0,0.6,0.45),rgb(0.94,0.89,0.25)
               ,rgb(0.8,0.47,0.65),rgb(0.9,0.6,0),
               rgb(0.33,0.7,0.9))
  color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),rgb(0.94,0.89,0.25,0.2)
               ,rgb(0.8,0.47,0.65,0.2),rgb(0.9,0.6,0,0.2),rgb(0.33,0.7,0.9,0.2) )
  windows()
  par(xpd=TRUE, mar=c(8,4,4,3))
  plot(NULL,xlim=c(0,date[length(date)]),ylim=c(0,1),xlab="time(day)",ylab="",main=model_name,
       , las=1, cex.lab=1.3,font=1)
  mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)
  for (i in 1:nb_sources) {polygon(x=c(aranged_data[[1]],rev(aranged_data[[1]]))
                                   ,y=c(aranged_data[[ind[i]+1]],rev(aranged_data[[ind[i]+2]]))
                                   ,col=color_poly[i],border=NA)}
  for (j in 1:nb_sources){points(aranged_data[[1]],aranged_data[[ind[j]]],type="l"
                                 ,col=color_main[j],
                                 pch=16,lwd=2)}
  
  legend(where_legend,legend=sources_name,cex=1,col=c(color_main[1:nb_sources]),lty="solid",box.lty=0)
}
