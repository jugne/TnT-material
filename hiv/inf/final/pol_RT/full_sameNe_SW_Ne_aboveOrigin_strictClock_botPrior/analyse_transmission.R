rm(list = ls())
library(stringr)
library(ggplot2)
library(PRROC)
library(pROC)

'%notin%' <- Negate('%in%')


inference_file <- ""
output_file <- ""
tnt_dir <- "~/Documents/Source/TnT-material/hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/combined/"
setwd(tnt_dir)

### Tools are part of TnT repository: https://github.com/jugne/TnT/
#https://github.com/jugne/TnT/blob/master/src/tnt/tntAnnotator/TransmissionAnalyser.java
transmission_analyser_path <- "../TransmissionAnalyser.jar"
# https://github.com/jugne/TnT/blob/master/scripts/summarize_transmission.py
summariser_path<- "../summarize_transmission.py"



cmd<-paste0("python3 ", summariser_path, " -i env_pol_saConstrained_s1_sameNe_combined.transmission.log -o ./")
system(cmd)

cmd<-paste0("python3 ", summariser_path, " -i ../../full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/combined/env_pol_saConstrained_s1_sameNe_noBot_combined.transmission.log -o ../../full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/combined/")
system(cmd)

trace_tnt<-read.table("env_pol_saConstrained_s1_sameNe_combined.log", header = T, check.names=T)
log_tnt <- read.csv("inferred_transmission.csv", header = T)
log_tnt_noBot <- read.csv("../../full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/combined/inferred_transmission.csv", header = T)
log_tnt_times <- read.table("env_pol_saConstrained_s1_sameNe_combined.transmission_infectionTimes.txt", header = T, check.names=FALSE)
log_tnt_times_no_bot <- read.table("../../full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/combined/env_pol_saConstrained_s1_sameNe_noBot_combined.transmission_infectionTimes.txt", header = T, check.names=FALSE)
ll<- read.table("env_pol_saConstrained_s1_sameNe_combined.transmission.log", header = T, check.names=FALSE)
ll_noBot <- read.table("../../full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/combined/env_pol_saConstrained_s1_sameNe_noBot_combined.transmission.log", header = T, check.names=FALSE)
true_times<- read.table("../../../../../transmission_dates.csv", header = T, sep = ",")
true_times<-true_times[1:9,]

e_low <- c()
e_high<-c()
e_median <-c()

e_no_bot_low <- c()
e_no_bot_high<-c()
e_no_bot_median <-c()
for (h in true_times$to){
  hpd<-HPDinterval(as.mcmc(log_tnt_times[[h]]))
  e_low<-c(e_low,hpd[1])
  e_high<-c(e_high,hpd[2])
  e_median<-c(e_median,median(log_tnt_times[[h]], na.rm=T))
  
  hpd_not_bot<-HPDinterval(as.mcmc(log_tnt_times_no_bot[[h]]))
  e_no_bot_low<-c(e_no_bot_low,hpd_not_bot[1])
  e_no_bot_high<-c(e_no_bot_high,hpd_not_bot[2])
  e_no_bot_median <- c(e_no_bot_median, median(log_tnt_times_no_bot[[h]], na.rm=T))
}



times<-data.frame(status=c(rep("true", length(true_times$to)),
                           rep("est", length(true_times$to)),
                           rep("est_no_bot", length(true_times$to))),
                  host=c(true_times$to, true_times$to, true_times$to),
                  low=c(true_times$low_date_distance_from_0, e_low, e_no_bot_low),
                  high=c(true_times$high_date_distance_from_0, e_high, e_no_bot_high),
                  median=c(rep(0, length(true_times$to)), e_median, e_no_bot_median))
times$status <- factor(times$status, levels = c("est_no_bot", "est", "true"))


pd <- position_dodge(width=0.4)

p<-ggplot(times, aes(color=status)) +
  geom_point(aes(x=2006.16-median, y=host ,shape=status),size=2, position=pd) +
  theme_minimal() +
  theme(axis.title = element_text(size = 14)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "Year", y = "Host") +
  geom_errorbar(aes(y=host, xmin=2006.16-low, xmax=2006.16-high, color=status),position=pd, width=0.2, size=1)+
  scale_color_manual(name="status",values=c("#FFBAA0", "coral", "red"),labels=c('TnT, no bottleneck', 'TnT', 'True')) +
  scale_shape_manual(name="status",values=c(19,19, -1), labels=c('TnT, no bottleneck', 'TnT', 'True'))
ggsave(plot=p,paste("infectionTimes2", ".pdf",
                                  sep=""),width=6, height=5)

hpd <- HPDinterval(as.mcmc(trace_tnt$pairwiseCoalProbBottleneck.env_pol))
ids <- which(trace_tnt$pairwiseCoalProbBottleneck.env_pol>=hpd[1] &
               trace_tnt$pairwiseCoalProbBottleneck.env_pol<=hpd[2] )
dd<- data.frame(value=c(trace_tnt$pairwiseCoalProbBottleneck.env_pol,
                        trace_tnt$pairwiseCoalProbBottleneck.env_pol[ids]),
                hpd=c(rep("no", nrow(trace_tnt)), rep("yes", length(ids))))
pp<-ggplot(dd)+aes(x=0,y=value, fill = hpd, alpha=hpd)+
  geom_violin(scale="width", position="identity", color = NA)+
  scale_fill_manual(values=c("#69b3a2", "#69b3a2"))+
  stat_summary(fun = "median", geom = "crossbar", width = 0.25,
               position = position_dodge(width = .25), color = "#5A5A5A", alpha=0.8)+
  scale_alpha_manual(values=c(0.3, 0.8))+
  theme_minimal() + theme(axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          legend.position = "None")+
  labs(x = "Pairwise Coalescent Probability", y="")

ggsave(plot=pp,paste("bottleneck", ".pdf",
                     sep=""),width=4, height=6)

tt <- data.frame(mean=unlist(lapply(trace_tnt[6:26], mean)),
                 median=unlist(lapply(trace_tnt[6:26], median)),
                 hpd_l=unlist(lapply(trace_tnt[6:26],
                              function(x) HPDinterval(as.mcmc(x))[1])),
                 hpd_h=unlist(lapply(trace_tnt[6:26],
                              function(x) HPDinterval(as.mcmc(x))[2])))
tt[,1:ncol(tt)]<-round(tt[,1:ncol(tt)],3)

write.csv(tt, "param_summary.csv")

