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







# Now use these data to fit the model -------------------------------------
lambda_mean=0.03
lambda_sd=0.01
alpha_lambda=lambda_mean^2/(lambda_sd)^2
beta_lambda=lambda_mean/(lambda_sd)^2
# alpha_lambda=1
# beta_lambda=20
#t_obs <- sort(jitter(oysters_dat$Time / 100, amount = 0.01))
t_obs <- sort(jitter(geese_dat$Time+0.01 / 1, amount = 0.01))
y_obs <- cbind(geese_dat$d13C_Pl, geese_dat$d15N_Pl) #
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

bdmm <- cmdstan_model("bdmm_one_distrib_lambda.stan", force_recompile = TRUE)
# bdmm_fitted <- bdmm$sample(
#   data = bdmm_data,
#   seed = 102,
#   chains = 4,
#   parallel_chains = 4,
#   refresh = 100,
#   iter_sampling	= 4000,
#   iter_warmup = 2000,
#   thin = 2
# )

bdmm_fitted <- bdmm$sample(
  data = bdmm_data,
  seed = 102,
  chains = 4,
  parallel_chains = 4,
  refresh = 100
)




write.csv(bdmm_fitted$summary(), file = "v12_GEESE_3_sources.csv", row.names = FALSE)

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
# write.csv(bdmm_generated_summary, file = "v14_ver_longue_bdmm_generated_summary_3sources.csv", row.names = FALSE)






----#CONTRIB ----

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


write.csv(p_1, file = "GEESE_p1_3s.csv", row.names = FALSE)
write.csv(p_2, file = "GEESE_p2_3s.csv", row.names = FALSE)
write.csv(p_3, file = "GEESE_p3_3s.csv", row.names = FALSE)



color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
             rgb(0,0.6,0.45),rgb(0.9,0.6,0),rgb(0.8,0.47,0.65),
             rgb(0.33,0.7,0.9),rgb(0.94,0.89,0.25))
color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),
             rgb(0.9,0.6,0,0.2),rgb(0.8,0.47,0.65,0.2),rgb(0.33,0.7,0.9,0.2),rgb(0.94,0.89,0.25,0.2) )

tiff("contrib_BDMM_geese_3.tiff",  width =4 , height = 5, units = "in", res = 300)
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
# polygon(x=c(p_1$time,rev(p_1$time))
#         ,y=c(p_4$q5,rev(p_4$q95))
#         ,col=color_poly[4],border=NA)

points(p_1$time,p_1$mean,type="l",col=color_main[1],pch=16,lwd=2)
points(p_1$time,p_2$mean,type="l",col=color_main[2],pch=16,lwd=2)
points(p_1$time,p_3$mean,type="l",col=color_main[3],pch=16,lwd=2)
# points(p_1$time,p_4$mean,type="l",col=color_main[4],pch=16,lwd=2)
dev.off()

----#Lambda####

lambda_fitted<-data.frame(bdmm_fitted$draws(c("lambda")))
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(lambda_fitted)

write.csv(lambda_fitted, file = "GEESE_3s_lambda_fitted.csv", row.names = FALSE)

for (i in 1:ncol(lambda_fitted)){
   # i=4
  val=data.frame(lambda=c(lambda_fitted[[i]],rgamma(n=1000,shape=alpha_lambda,rate=beta_lambda))
                                                     ,type=as.factor(c(rep("posterior",1000),rep("prior",1000))))
  title=paste0('lambda_geese_3',i,'.tiff')
  ggplot(val, aes(x=lambda, color=type)) +
    xlab("λ (d-1)") + ylab("Density")+
    geom_density()+
    ggtitle(i)+
    theme_classic()
  ggsave(title)
  }

----#fit to data####

consumer_summary <- bdmm_generated$summary() |>
  dplyr::filter(variable |> stringr::str_detect("^y_fit")) %>% 
  tidyr::separate(variable, c("observation", "isotope"), sep = ",") %>% 
  tidyr::pivot_wider(names_from = "isotope", values_from = c("mean", "median", "sd", "mad", "q5", "q95"))
consumer_summary$time <- bdmm_data$t_fit
write.csv(consumer_summary, file = "GEESE_fit_consumer.csv", row.names = FALSE)

#carbon
 
tiff("fit_geese_C.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,tail(consumer_summary$time,1)),ylim=c(-30,10),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(consumer_summary$time,rev(consumer_summary$time))
        ,y=c(consumer_summary$`q5_1]`,rev(consumer_summary$`q95_1]`))
        ,col=rgb(0.5,0.5,0.5,0.6),border=NA)
points(consumer_summary$time,consumer_summary$`mean_1]`,pch=16,col="grey40",cex=1.5)
points(geese_dat$Time,geese_dat$d13C_Pl,pch=16,col="red")
dev.off()
 

#nitrogen

tiff("fit_geese_N.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,tail(consumer_summary$time,1)),ylim=c(-10,15),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(consumer_summary$time,rev(consumer_summary$time))
        ,y=c(consumer_summary$`q5_2]`,rev(consumer_summary$`q95_2]`))
        ,col=rgb(0.5,0.5,0.5,0.6),border=NA)
points(consumer_summary$time,consumer_summary$`mean_2]`,pch=16,col="grey40",cex=1.5)
points(geese_dat$Time,geese_dat$d15N_Pl,pch=16,col="red")
dev.off()

#----FIT v2----

y_fitted=data.frame(bdmm_fitted$draws(c("y_fit")))
#y fit for isotope 1 and chain 1
chain_1=seq(1,1004,by=4)
y_1_ch_1=y_fitted[,chain_1]
stats_y1_ch1=data.frame(row.names = NULL)
for (i in 1:ncol(y_1_ch_1)) {
 y_val= y_1_ch_1[,i]
 probs=quantile(y_val,probs=c(0.05,0.95))
 stats_y1_ch1=rbind(stats_y1_ch1,data.frame(M=mean(y_val),q5=probs[1],q95=probs[2]))
}



###

consumer_fit <- bdmm_fitted$summary() |>
  dplyr::filter(variable |> stringr::str_detect("^y_fit")) %>% 
  tidyr::separate(variable, c("observation", "isotope"), sep = ",") %>% 
  tidyr::pivot_wider(names_from = "isotope", values_from = c("mean", "median", "sd", "mad", "q5", "q95"))
consumer_fit$time <- rep(bdmm_data$t_obs,2)
consumer_fit$iso <- c(rep(1,251),rep(2,251))



#carbon
conso_c=subset(consumer_fit,consumer_fit$iso==1)
tiff("fit_geese_C_3_s.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,tail(conso_c$time,1)),ylim=c(-30,10),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(conso_c$time,rev(conso_c$time))
        ,y=c(conso_c$`q5_1]`,rev(conso_c$`q95_1]`))
        ,col=rgb(0.5,0.5,0.5,0.6),border=NA)
points(conso_c$time,conso_c$`mean_1]`,pch=16,col="grey40",cex=1.5)
points(geese_dat$Time,geese_dat$d13C_Pl,pch=16,col="red")
dev.off()

#nitrogen

conso_n=subset(consumer_fit,consumer_fit$iso==2)
tiff("fit_geese_N_3_s.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,tail(conso_n$time,1)),ylim=c(-10,15),xlab="Time(d)",ylab="",main="",
     , las=1, cex.lab=1.3,font=1)
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=1)
polygon(x=c(conso_n$time,rev(conso_n$time))
        ,y=c(conso_n$`q5_2]`,rev(conso_n$`q95_2]`))
        ,col=rgb(0.5,0.5,0.5,0.6),border=NA)
points(conso_n$time,conso_n$`mean_2]`,pch=16,col="grey40",cex=1.5)
points(geese_dat$Time,geese_dat$d15N_Pl,pch=16,col="red")

dev.off()

