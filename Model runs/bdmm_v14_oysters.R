# Version 14 has different sources for different time points
# FIT TO THE OYSTERS DATA
# And proper source values
# This version contains concentration dependence
rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/3eme année/2eme papier/data")


# Load in packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(cmdstanr)
options(mc.cores = parallel::detectCores())

# Open teh oysters data
library(readxl)
oysters_dat <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'Targets')
oysters_sources <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'Sources')
oysters_sources_MeanC <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'Sources_MeanC') %>% select(-Time) %>% as.matrix()
oysters_sources_MeanN <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'Sources_MeanN') %>% select(-Time) %>% as.matrix()
oysters_sources_SdC <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'Sources_SdC') %>% select(-Time) %>% as.matrix()
oysters_sources_SdN <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'Sources_SdN') %>% select(-Time) %>% as.matrix()
oysters_TEF <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'TEFs')
oysters_TEF_MeanC <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'TEFs_MeanC') %>% select(-Time) %>% as.matrix()
oysters_TEF_MeanN <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'TEFs_MeanN') %>% select(-Time) %>% as.matrix()
oysters_TEF_SdC <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'TEFs_SdC') %>% select(-Time) %>% as.matrix()
oysters_TEF_SdN <- read_excel("Donnees_Marin_Leal_use_3sources.xlsx", sheet = 'TEFs_SdN') %>% select(-Time) %>% as.matrix()

# Now use these data to fit the model -------------------------------------

#t_obs <- sort(jitter(oysters_dat$Time / 100, amount = 0.01))
t_obs <- sort(jitter(oysters_dat$Time+0.01 / 1, amount = 0.01))
y_obs <- cbind(oysters_dat$d13C_Pl, oysters_dat$d15N_Pl) # ADD NOISE HERE
# t_fit <- seq(min(t_obs), max(t_obs), length = 20) # Can't do t_fit on a new grid unless we interpolate the source values
t_fit <- t_obs
bdmm_data <- list(
  n_obs = length(t_obs),
  n_fit = length(t_fit),
  K = length(unique(oysters_TEF$Source)),
  D = 3,
  y_obs = y_obs,
  t_obs = t_obs,
  t_fit = t_fit,
  t0 = 0,
  mu_s1 = oysters_sources_MeanC + oysters_TEF_MeanC,
  mu_s2 = oysters_sources_MeanN + oysters_TEF_MeanN,
  sigma_s1 = sqrt(oysters_sources_SdC^2 + oysters_TEF_SdC^2),
  sigma_s2 = sqrt(oysters_sources_SdN^2 + oysters_TEF_SdN^2),
  lambda_mean=as.numeric(oysters_dat$lambda),
  alpha_lambda = as.numeric(oysters_dat$alpha), # Units are inverse days! Used 5 when Time was divided by 100
  beta_lambda = as.numeric(oysters_dat$beta),
  q1 = rep(1, length(unique(oysters_TEF$Source))), # No concentration dependence
  q2 = rep(1, length(unique(oysters_TEF$Source))), # No concentration dependence
  t_mean = mean(t_obs),
  t_sd =  sd(t_obs)
)

bdmm <- cmdstan_model("bdmm_distribs_lambda.stan", force_recompile = TRUE)
bdmm_fitted <- bdmm$sample(
  data = bdmm_data,
  seed = 102,
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  iter_sampling	= 4000,
  iter_warmup = 2000,
  thin = 2
)

# bdmm_fitted <- bdmm$sample(
#   data = bdmm_data,
#   seed = 102,
#   chains = 4,
#   parallel_chains = 4,
#   refresh = 100
# )

write.csv(bdmm_fitted$summary(), file = "oysters_fitted_run2.csv", row.names = FALSE)

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
write.csv(bdmm_generated_summary, file = "oysters_generated_run2.csv", row.names = FALSE)



##----FIT GRPAHS----

#carbon
dev.off() 
tiff("fit_BDMM_C_oysters_run2.tiff", , width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,tail(bdmm_generated_summary$time,1)),ylim=c(-22,-17),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(bdmm_generated_summary$time,rev(bdmm_generated_summary$time))
        ,y=c(bdmm_generated_summary$`q5_1]`,rev(bdmm_generated_summary$`q95_1]`))
        ,col=rgb(0.5,0.5,0.5,0.6),border=NA)
points(bdmm_generated_summary$time,bdmm_generated_summary$`mean_1]`,pch=16,col="grey40",cex=1.5)
points(oysters_dat$Time,oysters_dat$d13C_Pl,pch=16,col="red")
dev.off()
#nitrogen
tiff("fit_BDMM_N_oysters_run2.tiff", , width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,tail(bdmm_generated_summary$time,1)),ylim=c(6,11),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(bdmm_generated_summary$time,rev(bdmm_generated_summary$time))
        ,y=c(bdmm_generated_summary$`q5_2]`,rev(bdmm_generated_summary$`q95_2]`))
        ,col=rgb(0.5,0.5,0.5,0.6),border=NA)
points(bdmm_generated_summary$time,bdmm_generated_summary$`mean_2]`,pch=16,col="grey40",cex=1.5)
points(oysters_dat$Time,oysters_dat$d15N_Pl,pch=16,col="red")
dev.off()

##----CONTRIBUTION GRAPHS----


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


tiff("Contrib_BDMM_oysters_run2.tiff", , width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(NULL,xlim=c(0,tail(p_1$time,1)),ylim=c(0,1),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)


polygon(x=c(p_1$time,rev(p_1$time))
        ,y=c(p_1$q5,rev(p_1$q95))
        ,col=color_poly[1],border=NA)
polygon(x=c(p_1$time,rev(p_1$time))
        ,y=c(p_2$q5,rev(p_2$q95))
        ,col=color_poly[2],border=NA)
polygon(x=c(p_1$time,rev(p_1$time))
        ,y=c(p_3$q5,rev(p_3$q95))
        ,col=color_poly[3],border=NA)

points(p_1$time,p_1$mean,type="l",col=color_main[1],pch=16,lwd=2)
points(p_1$time,p_2$mean,type="l",col=color_main[2],pch=16,lwd=2)
points(p_1$time,p_3$mean,type="l",col=color_main[3],pch=16,lwd=2)
dev.off()


##----LAMBDA GRAPHS----
#posterior
# lambda_fitted=data.frame(bdmm_fitted$draws(c("lambda")))
# write.csv2(lambda_fitted, file = "lambda_fitted_oysters.csv", row.names = FALSE)

lambda_fitted= read.csv2("lambda_fitted_oysters.csv")
rep_dat_oyst<- oysters_dat[rep(1:nrow(oysters_dat), each = 4), ]

cols<-c("orange","skyblue","black")

for (i in 1:ncol(lambda_fitted)){
   # i=1
  #posterior
  post= lambda_fitted[[i]]
  dens_post<- density(post)
  #prior
  prior= rgamma(n=nrow(lambda_fitted),shape=rep_dat_oyst$alpha[i],rate=rep_dat_oyst$beta[i])
  dens_prior<-density(prior)
  
  y_max=max(dens_prior[[2]],dens_post[[2]])
  title_graph=paste0("lambda_oysters_t_",rep_dat_oyst$Time[i],"_num",i,".tiff")
  tiff(title_graph,  width =4 , height = 5, units = "in", res = 300)
  plot(dens_prior,col=cols[1], xlab="",ylab="",main="",xlim=c(0,0.03)
       ,lwd=3, cex.lab=1.3,las=1,ylim=c(0,y_max)) 
  mtext("Density",side=2,line=4,padj=1.5,cex=1.3,font=1)
  mtext("λ (d-1)",side=1,line=3,cex=1.3,font=1)
  points(dens_post,type="l",col=cols[2],lwd=3)
  dev.off()

}



plot(NULL)

##----GRAPHS source poly----
temps=unique(oysters_dat$Time)
sources_corr_C=data.frame(cbind(Time=oysters_dat$Time,bdmm_data$mu_s1))
sources_corr_N=data.frame(cbind(Time=oysters_dat$Time,bdmm_data$mu_s2))
sources_corr_sd_C=data.frame(cbind(Time=oysters_dat$Time,bdmm_data$sigma_s1))
sources_corr_sd_N=data.frame(cbind(Time=oysters_dat$Time,bdmm_data$sigma_s2))


for(i in 1:length(temps)){

t=temps[i]
subs_conso=subset(oysters_dat,oysters_dat$Time==t)
s_corr_C=subset(sources_corr_C,sources_corr_C$Time==t)[1,]
s_corr_N=subset(sources_corr_N,sources_corr_N$Time==t)[1,]
sd_corr_C=subset(sources_corr_sd_C,sources_corr_sd_C$Time==t)[1,]
sd_corr_N=subset(sources_corr_sd_N,sources_corr_sd_N$Time==t)[1,]

title_graph=paste0("oysters_biplot_t_",t,'.tiff')  
tiff(title_graph,  width =4 , height = 5, units = "in", res = 300)
plot(NULL,main=paste("t=",t,"d"),xlim=c(-30,-14),ylim=c(5,14),xlab="",ylab="",cex.lab=1.4,font.lab=2,las=1) #plot frame
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
dev.off()

}


tiff("legend_biplot_oyster.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,main="",xlim=c(-30,-14),ylim=c(5,14),xlab="",ylab="",cex.lab=1.4,font.lab=2,las=1) #plot frame
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=-1,padj=1.5,cex=1.5,font=2)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=-2,cex=1.5,font=2)
legend("center",legend=c("MPBxU","POM","TOM","Consumer"),cex=1.5,col=c(color_main[1:3],"black")
       ,lty=c(1,1,1,NA),box.lty=0,pch=c(NA,NA,NA,16))
dev.off()


