######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)


# clear workspace
rm(list = ls())

test_vals <- function(prob, true=0, estimated=c()){
  lower <- HPDinterval(as.mcmc(estimated), prob=prob)[1]
  upper <- HPDinterval(as.mcmc(estimated), prob=prob)[2]
  test <- c(as.numeric(true >= lower & true <= upper))
  return(test)
}



# Set the directory to the directory of the file
this.dir<-"~/Documents/Source/TnT-material/validation_3009/runs/"
setwd(this.dir)

# read in the true rates
true_rates_file <- "../true_rates.csv"
true.rates <- read.table(true_rates_file, header=TRUE, sep=",")

# read in the true rates
sampling_file <- "../summarySamples.csv"
true.sampling <- read.table(sampling_file, header=TRUE, sep=",")


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

first = T
n<-200
remove_rows<-c()
m<-c()

pairwiseCoalProbBottleneck = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n), lower=numeric(n), rel.err.meadian=numeric(n),
                                        rel.err.mean=numeric(n), rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n), hpd.rel.width=numeric(n), test=numeric(n),
                                        test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                        test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                        test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                        test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                        test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                        test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                        test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

R0 = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n), lower=numeric(n), rel.err.meadian=numeric(n),
                       rel.err.mean=numeric(n), rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n), hpd.rel.width=numeric(n), test=numeric(n),
                test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

delta = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n), lower=numeric(n), rel.err.meadian=numeric(n),
                rel.err.mean=numeric(n), rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n), hpd.rel.width=numeric(n), test=numeric(n),
                test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))


popSize.1 = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n), lower=numeric(n), rel.err.meadian=numeric(n),
                       rel.err.mean=numeric(n), rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n), hpd.rel.width=numeric(n), test=numeric(n),
                       test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                       test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                       test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                       test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                       test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                       test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                       test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

popSize.2 = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n), lower=numeric(n), rel.err.meadian=numeric(n),
                       rel.err.mean=numeric(n), rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n), hpd.rel.width=numeric(n), test=numeric(n),
                       test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                       test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                       test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                       test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                       test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                       test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                       test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

removalProb = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n), lower=numeric(n), rel.err.meadian=numeric(n),
                       rel.err.mean=numeric(n), rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n), hpd.rel.width=numeric(n), test=numeric(n),
                       test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                       test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                       test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                       test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                       test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                       test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                       test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

rsamplingAtPresentProb = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n), lower=numeric(n), rel.err.meadian=numeric(n),
                                    rel.err.mean=numeric(n), rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n), hpd.rel.width=numeric(n), test=numeric(n),
                                    test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                    test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                    test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                    test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                    test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                    test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                    test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))


post_ess<- numeric(n)
for (i in 1:200){
  
  log <- list.files(path=paste0("run_",i,"/inf"), pattern="inference.log", full.names = TRUE)
  
  # read in the log file
  if (length(log)==0){
    remove_rows = cbind(remove_rows, i)
    next
  }
  
  used.rates = true.rates[which(true.rates$X==i),];
  
  t <- read.table(log[[1]], header=TRUE, sep="\t")
  
  # if (length(t[,1])<4001){
  #   m = cbind(m, i)
  #   next
  # }
  
  
  if (length(t$Sample)<2){
    remove_rows = cbind(remove_rows, i)
    next
  }
  
  # take a 10% burnin
  t <- t[-seq(1,ceiling(length(t$multiCoalescent)/10)), ]

  
  ess<-  effectiveSize(t)
  post_ess[i] = as.numeric(ess["posterior"])

  # if(post_ess[i]<20){
  #   remove_rows = cbind(remove_rows, i)
  #   next
  # }

  pairwiseCoalProbBottleneck$true[i] <- used.rates$pairwiseCoalProbBottleneck
  pairwiseCoalProbBottleneck$estimated[i] <- median(t$pairwiseCoalProbBottleneck)
  pairwiseCoalProbBottleneck$upper[i] <- quantile(t$pairwiseCoalProbBottleneck,0.975)
  pairwiseCoalProbBottleneck$lower[i] <- quantile(t$pairwiseCoalProbBottleneck,0.025)


  pairwiseCoalProbBottleneck$rel.err.meadian[i] <- abs(median(t$pairwiseCoalProbBottleneck)-used.rates$pairwiseCoalProbBottleneck)/used.rates$pairwiseCoalProbBottleneck
  pairwiseCoalProbBottleneck$rel.err.mean[i] <-abs(mean(t$pairwiseCoalProbBottleneck)-used.rates$pairwiseCoalProbBottleneck)/used.rates$pairwiseCoalProbBottleneck
  pairwiseCoalProbBottleneck$rel.err.mode[i] <- abs(Mode(t$pairwiseCoalProbBottleneck)-used.rates$pairwiseCoalProbBottleneck)/used.rates$pairwiseCoalProbBottleneck
  pairwiseCoalProbBottleneck$cv[i] <- sqrt(exp(sd(log(t$pairwiseCoalProbBottleneck))**2)-1)
  pairwiseCoalProbBottleneck$hpd.lower[i] <- HPDinterval(as.mcmc(t$pairwiseCoalProbBottleneck))[1]
  pairwiseCoalProbBottleneck$hpd.upper[i] <- HPDinterval(as.mcmc(t$pairwiseCoalProbBottleneck))[2]
  pairwiseCoalProbBottleneck$test[i] <- c(as.numeric(used.rates$pairwiseCoalProbBottleneck >= pairwiseCoalProbBottleneck$hpd.lower[i] & used.rates$pairwiseCoalProbBottleneck <= pairwiseCoalProbBottleneck$hpd.upper[i]))
  pairwiseCoalProbBottleneck$hpd.rel.width[i] <- (pairwiseCoalProbBottleneck$hpd.upper[i]-pairwiseCoalProbBottleneck$hpd.lower[i])/used.rates$pairwiseCoalProbBottleneck
  pairwiseCoalProbBottleneck[i,(ncol(pairwiseCoalProbBottleneck)-21+1):ncol(pairwiseCoalProbBottleneck)] <- lapply(seq(0.,1,0.05), test_vals, pairwiseCoalProbBottleneck$true[i], t$pairwiseCoalProbBottleneck)

  popSize.1$true[i] <- used.rates$Ne
  popSize.1$estimated[i] <- median(t$popSize.1)
  popSize.1$upper[i] <- quantile(t$popSize.1,0.975)
  popSize.1$lower[i] <- quantile(t$popSize.1,0.025)


  popSize.1$rel.err.meadian[i] <- abs(median(t$popSize.1)-used.rates$Ne)/used.rates$Ne
  popSize.1$rel.err.mean[i] <-abs(mean(t$popSize.1)-used.rates$Ne)/used.rates$Ne
  popSize.1$rel.err.mode[i] <- abs(Mode(t$popSize.1)-used.rates$Ne)/used.rates$Ne
  popSize.1$cv[i] <- sqrt(exp(sd(log(t$popSize.1))**2)-1)
  popSize.1$hpd.lower[i] <- HPDinterval(as.mcmc(t$popSize.1))[1]
  popSize.1$hpd.upper[i] <- HPDinterval(as.mcmc(t$popSize.1))[2]
  popSize.1$test[i] <- c(as.numeric(used.rates$Ne >= popSize.1$hpd.lower[i] & used.rates$Ne <= popSize.1$hpd.upper[i]))
  popSize.1$hpd.rel.width[i] <- (popSize.1$hpd.upper[i]-popSize.1$hpd.lower[i])/used.rates$Ne
  popSize.1[i,(ncol(popSize.1)-21+1):ncol(popSize.1)] <- lapply(seq(0.,1,0.05), test_vals, popSize.1$true[i], t$popSize.1)
  

  popSize.2$true[i] <- used.rates$Ne_orig
  popSize.2$estimated[i] <- median(t$popSize.2)
  popSize.2$upper[i] <- quantile(t$popSize.2,0.975)
  popSize.2$lower[i] <- quantile(t$popSize.2,0.025)


  popSize.2$rel.err.meadian[i] <- abs(median(t$popSize.2)-used.rates$Ne_orig)/used.rates$Ne_orig
  popSize.2$rel.err.mean[i] <-abs(mean(t$popSize.2)-used.rates$Ne_orig)/used.rates$Ne_orig
  popSize.2$rel.err.mode[i] <- abs(Mode(t$popSize.2)-used.rates$Ne_orig)/used.rates$Ne_orig
  popSize.2$cv[i] <- sqrt(exp(sd(log(t$popSize.2))**2)-1)
  popSize.2$hpd.lower[i] <- HPDinterval(as.mcmc(t$popSize.2))[1]
  popSize.2$hpd.upper[i] <- HPDinterval(as.mcmc(t$popSize.2))[2]
  popSize.2$test[i] <- c(as.numeric(used.rates$Ne_orig >= popSize.2$hpd.lower[i] & used.rates$Ne_orig <= popSize.2$hpd.upper[i]))
  popSize.2$hpd.rel.width[i] <- (popSize.2$hpd.upper[i]-popSize.2$hpd.lower[i])/used.rates$Ne_orig
  popSize.2[i,(ncol(popSize.2)-21+1):ncol(popSize.2)] <- lapply(seq(0.,1,0.05), test_vals, popSize.2$true[i], t$popSize.2)
  

  R0$true[i] <- used.rates$R0
  R0$estimated[i] <- median(t$R0)
  R0$upper[i] <- quantile(t$R0,0.975)
  R0$lower[i] <- quantile(t$R0,0.025)


  R0$rel.err.meadian[i] <- abs(median(t$R0)-used.rates$R0)/used.rates$R0
  R0$rel.err.mean[i] <-abs(mean(t$R0)-used.rates$R0)/used.rates$R0
  R0$rel.err.mode[i] <- abs(Mode(t$R0)-used.rates$R0)/used.rates$R0
  R0$cv[i] <- sqrt(exp(sd(log(t$R0))**2)-1)
  R0$hpd.lower[i] <- HPDinterval(as.mcmc(t$R0))[1]
  R0$hpd.upper[i] <- HPDinterval(as.mcmc(t$R0))[2]
  R0$test[i] <- c(as.numeric(used.rates$R0 >= R0$hpd.lower[i] & used.rates$R0 <= R0$hpd.upper[i]))
  R0$hpd.rel.width[i] <- (R0$hpd.upper[i]-R0$hpd.lower[i])/used.rates$R0
  R0[i,(ncol(R0)-21+1):ncol(R0)] <- lapply(seq(0.,1,0.05), test_vals, R0$true[i], t$R0)
  
  
  delta$true[i] <- used.rates$delta
  delta$estimated[i] <- median(t$becominUninfectiousRate)
  delta$upper[i] <- quantile(t$becominUninfectiousRate,0.975)
  delta$lower[i] <- quantile(t$becominUninfectiousRate,0.025)


  delta$rel.err.meadian[i] <- abs(median(t$becominUninfectiousRate)-used.rates$delta)/used.rates$delta
  delta$rel.err.mean[i] <-abs(mean(t$becominUninfectiousRate)-used.rates$delta)/used.rates$delta
  delta$rel.err.mode[i] <- abs(Mode(t$becominUninfectiousRate)-used.rates$delta)/used.rates$delta
  delta$cv[i] <- sqrt(exp(sd(log(t$becominUninfectiousRate))**2)-1)
  delta$hpd.lower[i] <- HPDinterval(as.mcmc(t$becominUninfectiousRate))[1]
  delta$hpd.upper[i] <- HPDinterval(as.mcmc(t$becominUninfectiousRate))[2]
  delta$test[i] <- c(as.numeric(used.rates$delta >= delta$hpd.lower[i] & used.rates$delta <= delta$hpd.upper[i]))
  delta$hpd.rel.width[i] <- (delta$hpd.upper[i]-delta$hpd.lower[i])/used.rates$delta
  delta[i,(ncol(delta)-21+1):ncol(delta)] <- lapply(seq(0.,1,0.05), test_vals, delta$true[i], t$becominUninfectiousRate)
  

  removalProb$true[i] <- used.rates$r
  removalProb$estimated[i] <- median(t$removalProb)
  removalProb$upper[i] <- quantile(t$removalProb,0.975)
  removalProb$lower[i] <- quantile(t$removalProb,0.025)


  removalProb$rel.err.meadian[i] <- abs(median(t$removalProb)-used.rates$r)/used.rates$r
  removalProb$rel.err.mean[i] <-abs(mean(t$removalProb)-used.rates$r)/used.rates$r
  removalProb$rel.err.mode[i] <- abs(Mode(t$removalProb)-used.rates$r)/used.rates$r
  removalProb$cv[i] <- sqrt(exp(sd(log(t$removalProb))**2)-1)
  removalProb$hpd.lower[i] <- HPDinterval(as.mcmc(t$removalProb))[1]
  removalProb$hpd.upper[i] <- HPDinterval(as.mcmc(t$removalProb))[2]
  removalProb$test[i] <- c(as.numeric(used.rates$r >= removalProb$hpd.lower[i] & used.rates$r <= removalProb$hpd.upper[i]))
  removalProb$hpd.rel.width[i] <- (removalProb$hpd.upper[i]-removalProb$hpd.lower[i])/used.rates$r
  removalProb[i,(ncol(removalProb)-21+1):ncol(removalProb)] <- lapply(seq(0.,1,0.05), test_vals, removalProb$true[i], t$removalProb)
  


  rsamplingAtPresentProb$true[i] <- used.rates$rho
  rsamplingAtPresentProb$estimated[i] <- median(t$rho)
  rsamplingAtPresentProb$upper[i] <- quantile(t$rho,0.975)
  rsamplingAtPresentProb$lower[i] <- quantile(t$rho,0.025)


  rsamplingAtPresentProb$rhoel.err.meadian[i] <- abs(median(t$rho)-used.rates$rho)/used.rates$rho
  rsamplingAtPresentProb$rhoel.err.mean[i] <-abs(mean(t$rho)-used.rates$rho)/used.rates$rho
  rsamplingAtPresentProb$rhoel.err.mode[i] <- abs(Mode(t$rho)-used.rates$rho)/used.rates$rho
  rsamplingAtPresentProb$cv[i] <- sqrt(exp(sd(log(t$rho))**2)-1)
  rsamplingAtPresentProb$hpd.lower[i] <- HPDinterval(as.mcmc(t$rho))[1]
  rsamplingAtPresentProb$hpd.upper[i] <- HPDinterval(as.mcmc(t$rho))[2]
  rsamplingAtPresentProb$test[i] <- c(as.numeric(used.rates$rho >= rsamplingAtPresentProb$hpd.lower[i] & used.rates$rho <= rsamplingAtPresentProb$hpd.upper[i]))
  rsamplingAtPresentProb$hpd.rel.width[i] <- (rsamplingAtPresentProb$hpd.upper[i]-rsamplingAtPresentProb$hpd.lower[i])/used.rates$rho
  rsamplingAtPresentProb[i,(ncol(rsamplingAtPresentProb)-21+1):ncol(rsamplingAtPresentProb)] <- lapply(seq(0.,1,0.05), test_vals, rsamplingAtPresentProb$true[i], t$rho)
  

}

if (!is.null(remove_rows)){
  pairwiseCoalProbBottleneck <- pairwiseCoalProbBottleneck[-remove_rows, ]
  R0 <- R0[-remove_rows, ]
  delta <- delta[-remove_rows, ]
  popSize.1 <- popSize.1[-remove_rows, ]
  popSize.2 <- popSize.2[-remove_rows, ]
  removalProb <- removalProb[-remove_rows, ]
  rsamplingAtPresentProb <- removalProb[-remove_rows, ]
}


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(pairwiseCoalProbBottleneck [,(ncol(pairwiseCoalProbBottleneck)-21+1):ncol(pairwiseCoalProbBottleneck)], sum)/length(pairwiseCoalProbBottleneck$test_0))

l <- length(pairwiseCoalProbBottleneck$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(pairwiseCoalProbBottleneck[,(ncol(pairwiseCoalProbBottleneck)-21+1):ncol(pairwiseCoalProbBottleneck)]))

p.qq_div_rate <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()+ theme(text = element_text(size=20), plot.title=element_text(size=19, hjust=0))
  
ggsave(plot=p.qq_div_rate,paste("../figures/qq_pairwiseCoalProbBottleneck", ".pdf", sep=""),width=5, height=5)

p.qq_div_rate <- p.qq_div_rate+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Pairwise coal. prob. at bottleneck")

# p.qq_div_rate <- p.qq_div_rate+
#   stat_summary(data=dd, aes(x=p, y=vals), fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
#   ggtitle("Diversification rate")

ggsave(plot=p.qq_div_rate,paste("../figures/qq_pairwiseCoalProbBottleneck_boot", ".pdf", sep=""),width=5, height=5)

x.min_pairwiseCoalProbBottleneck = min(pairwiseCoalProbBottleneck$true)
x.max_pairwiseCoalProbBottleneck = max(pairwiseCoalProbBottleneck$true)

y.min_pairwiseCoalProbBottleneck = min(pairwiseCoalProbBottleneck$lower)
y.max_pairwiseCoalProbBottleneck = max(pairwiseCoalProbBottleneck$upper)


lim.min = min(x.min_pairwiseCoalProbBottleneck, y.min_pairwiseCoalProbBottleneck)
lim.max = max(x.max_pairwiseCoalProbBottleneck, y.max_pairwiseCoalProbBottleneck)

p.pairwiseCoalProbBottleneck <- ggplot(pairwiseCoalProbBottleneck)+
    geom_abline(intercept = 0, color="red", linetype="dashed")+
    geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
    geom_point(aes(x=true, y=estimated), size=2) +
    theme_minimal()
    # scale_y_log10(limits=c(lim.min, lim.max)) +
    # scale_x_log10(limits=c(lim.min, lim.max))



ggsave(plot=p.pairwiseCoalProbBottleneck,paste("../figures/pairwiseCoalProbBottleneck", ".pdf", sep=""),width=6, height=5)

p.pairwiseCoalProbBottleneck_log <- ggplot(pairwiseCoalProbBottleneck)+
    geom_abline(intercept = 0, color="red", linetype="dashed")+
    geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
    geom_point(aes(x=true, y=estimated), size=2) +
    theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

  ggsave(plot=p.pairwiseCoalProbBottleneck_log,paste("../figures/pairwiseCoalProbBottleneck_logScale", ".pdf", sep=""),width=6, height=5)

  hpd.test = data.frame("Parameter"=c('pairwiseCoalProbBottleneck'),
                        "HPD coverage"=c(mean(pairwiseCoalProbBottleneck$test)))
  write.csv(hpd.test, file = paste("../figures/pairwiseCoalProbBottleneck_HPD_test",".csv", sep="" ))

  hpd.width = data.frame("Parameter"=c('pairwiseCoalProbBottleneck'),
                         "Average relative HPD width"=c(mean(pairwiseCoalProbBottleneck$hpd.rel.width)))
  write.csv(hpd.width, file = paste("../figures/pairwiseCoalProbBottleneck_HPD_rel_width_score_",".csv", sep="" ))

  cv = data.frame("Parameter"=c('pairwiseCoalProbBottleneck'),
                  "Average Coefficient of Variation"=c(mean(pairwiseCoalProbBottleneck$cv)))
  write.csv(cv, file = paste("../figures/pairwiseCoalProbBottleneck_CV_score",".csv", sep="" ))

  medians = data.frame("Parameter"=c('pairwiseCoalProbBottleneck'),
                       "Medians"=c(mean(pairwiseCoalProbBottleneck$rel.err.meadian)))
  write.csv(medians, file = paste("../figures/pairwiseCoalProbBottleneck_rel_error_median",".csv", sep="" ))

  means = data.frame("Parameter=c('pairwiseCoalProbBottleneck')", "Means"=c(mean(pairwiseCoalProbBottleneck$rel.err.mean)))
  write.csv(means, file = paste("../figures/pairwiseCoalProbBottleneck_rel_error_mean",".csv", sep=""))


dt_rel.err.median.pairwiseCoalProbBottleneck = data.frame(param='pairwiseCoalProbBottleneck', value=pairwiseCoalProbBottleneck$rel.err.meadian)
p_rel_err_pairwiseCoalProbBottleneck<- ggplot(dt_rel.err.median.pairwiseCoalProbBottleneck, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.pairwiseCoalProbBottleneck$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_pairwiseCoalProbBottleneck,paste("../figures/pairwiseCoalProbBottleneck_rel_err_median", ".pdf", sep=""),width=6, height=5)

dt_rel.err.mean.pairwiseCoalProbBottleneck = data.frame(param='pairwiseCoalProbBottleneck', value=pairwiseCoalProbBottleneck$rel.err.mean)
p_rel_err_mean_pairwiseCoalProbBottleneck<- ggplot(dt_rel.err.mean.pairwiseCoalProbBottleneck, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.pairwiseCoalProbBottleneck$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_mean_pairwiseCoalProbBottleneck,paste("../figures/pairwiseCoalProbBottleneck_rel_err_mean", ".pdf", sep=""),width=6, height=5)

dt_hpd.rel.width.pairwiseCoalProbBottleneck = data.frame(param='pairwiseCoalProbBottleneck', value=pairwiseCoalProbBottleneck$hpd.rel.width)
p_hpd.rel.width.pairwiseCoalProbBottleneck <- ggplot(dt_hpd.rel.width.pairwiseCoalProbBottleneck, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.pairwiseCoalProbBottleneck$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_hpd.rel.width.pairwiseCoalProbBottleneck,paste("../figures/pairwiseCoalProbBottleneck_rel_hpd_width", ".pdf", sep=""),width=6, height=5)




########## Ne ################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(popSize.1[,(ncol(popSize.1)-21+1):ncol(popSize.1)], sum)/length(popSize.1$test_0))

l <- length(popSize.1$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(popSize.1[,(ncol(popSize.1)-21+1):ncol(popSize.1)]))

p.qq_div_rate <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()+ theme(text = element_text(size=20), plot.title=element_text(size=19, hjust=0))
ggsave(plot=p.qq_div_rate,paste("../figures/qq_popSize.1", ".pdf", sep=""),width=5, height=5)

p.qq_div_rate <- p.qq_div_rate+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Eff. pop. size, below origin")

# p.qq_div_rate <- p.qq_div_rate+
#   stat_summary(data=dd, aes(x=p, y=vals), fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
#   ggtitle("Diversification rate")

ggsave(plot=p.qq_div_rate,paste("../figures/qq_popSize.1_boot", ".pdf", sep=""),width=5, height=5)


x.min_popSize.1 = min(popSize.1$true)
x.max_popSize.1 = max(popSize.1$true)

y.min_popSize.1 = min(popSize.1$lower)
y.max_popSize.1 = max(popSize.1$upper)


lim.min = min(x.min_popSize.1, y.min_popSize.1)
lim.max = max(x.max_popSize.1, y.max_popSize.1)

p.popSize.1 <- ggplot(popSize.1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))



ggsave(plot=p.popSize.1,paste("../figures/popSize.1", ".pdf", sep=""),width=6, height=5)

p.popSize.1_log <- ggplot(popSize.1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.popSize.1_log,paste("../figures/popSize.1_logScale", ".pdf", sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('popSize.1'),
                      "HPD coverage"=c(mean(popSize.1$test)))
write.csv(hpd.test, file = paste("../figures/popSize.1_HPD_test",".csv", sep="" ))

hpd.width = data.frame("Parameter"=c('popSize.1'),
                       "Average relative HPD width"=c(mean(popSize.1$hpd.rel.width)))
write.csv(hpd.width, file = paste("../figures/popSize.1_HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('popSize.1'),
                "Average Coefficient of Variation"=c(mean(popSize.1$cv)))
write.csv(cv, file = paste("../figures/popSize.1_CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('popSize.1'),
                     "Medians"=c(mean(popSize.1$rel.err.meadian)))
write.csv(medians, file = paste("../figures/popSize.1_rel_error_median",".csv", sep="" ))

means = data.frame("Parameter=c('popSize.1')", "Means"=c(mean(popSize.1$rel.err.mean)))
write.csv(means, file = paste("../figures/popSize.1_rel_error_mean",".csv", sep=""))


dt_rel.err.median.popSize.1 = data.frame(param='popSize.1', value=popSize.1$rel.err.meadian)
p_rel_err_popSize.1<- ggplot(dt_rel.err.median.popSize.1, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.popSize.1$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_popSize.1,paste("../figures/popSize.1_rel_err_median", ".pdf", sep=""),width=6, height=5)

dt_rel.err.mean.popSize.1 = data.frame(param='popSize.1', value=popSize.1$rel.err.mean)
p_rel_err_mean_popSize.1<- ggplot(dt_rel.err.mean.popSize.1, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.popSize.1$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_mean_popSize.1,paste("../figures/popSize.1_rel_err_mean", ".pdf", sep=""),width=6, height=5)

dt_hpd.rel.width.popSize.1 = data.frame(param='popSize.1', value=popSize.1$hpd.rel.width)
p_hpd.rel.width.popSize.1 <- ggplot(dt_hpd.rel.width.popSize.1, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.popSize.1$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_hpd.rel.width.popSize.1,paste("../figures/popSize.1_rel_hpd_width", ".pdf", sep=""),width=6, height=5)


########## Ne_orig ################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(popSize.2[,(ncol(popSize.2)-21+1):ncol(popSize.2)], sum)/length(popSize.2$test_0))

l <- length(popSize.2$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(popSize.2[,(ncol(popSize.2)-21+1):ncol(popSize.2)]))

p.qq_div_rate <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()+ theme(text = element_text(size=20), plot.title=element_text(size=19, hjust=0))
ggsave(plot=p.qq_div_rate,paste("../figures/qq_popSize.2", ".pdf", sep=""),width=5, height=5)

p.qq_div_rate <- p.qq_div_rate+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Eff. pop. size, above origin")

ggsave(plot=p.qq_div_rate,paste("../figures/qq_popSize.2_boot", ".pdf", sep=""),width=5, height=5)

x.min_popSize.2 = min(popSize.2$true)
x.max_popSize.2 = max(popSize.2$true)

y.min_popSize.2 = min(popSize.2$lower)
y.max_popSize.2 = max(popSize.2$upper)


lim.min = min(x.min_popSize.2, y.min_popSize.2)
lim.max = max(x.max_popSize.2, y.max_popSize.2)

p.popSize.2 <- ggplot(popSize.2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))



ggsave(plot=p.popSize.2,paste("../figures/popSize.2", ".pdf", sep=""),width=6, height=5)

p.popSize.2_log <- ggplot(popSize.2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.popSize.2_log,paste("../figures/popSize.2_logScale", ".pdf", sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('popSize.2'),
                      "HPD coverage"=c(mean(popSize.2$test)))
write.csv(hpd.test, file = paste("../figures/popSize.2_HPD_test",".csv", sep="" ))

hpd.width = data.frame("Parameter"=c('popSize.2'),
                       "Average relative HPD width"=c(mean(popSize.2$hpd.rel.width)))
write.csv(hpd.width, file = paste("../figures/popSize.2_HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('popSize.2'),
                "Average Coefficient of Variation"=c(mean(popSize.2$cv)))
write.csv(cv, file = paste("../figures/popSize.2_CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('popSize.2'),
                     "Medians"=c(mean(popSize.2$rel.err.meadian)))
write.csv(medians, file = paste("../figures/popSize.2_rel_error_median",".csv", sep="" ))

means = data.frame("Parameter=c('popSize.2')", "Means"=c(mean(popSize.2$rel.err.mean)))
write.csv(means, file = paste("../figures/popSize.2_rel_error_mean",".csv", sep=""))


dt_rel.err.median.popSize.2 = data.frame(param='popSize.2', value=popSize.2$rel.err.meadian)
p_rel_err_popSize.2<- ggplot(dt_rel.err.median.popSize.2, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.popSize.2$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_popSize.2,paste("../figures/popSize.2_rel_err_median", ".pdf", sep=""),width=6, height=5)

dt_rel.err.mean.popSize.2 = data.frame(param='popSize.2', value=popSize.2$rel.err.mean)
p_rel_err_mean_popSize.2<- ggplot(dt_rel.err.mean.popSize.2, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.popSize.2$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_mean_popSize.2,paste("../figures/popSize.2_rel_err_mean", ".pdf", sep=""),width=6, height=5)

dt_hpd.rel.width.popSize.2 = data.frame(param='popSize.2', value=popSize.2$hpd.rel.width)
p_hpd.rel.width.popSize.2 <- ggplot(dt_hpd.rel.width.popSize.2, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.popSize.2$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_hpd.rel.width.popSize.2,paste("../figures/popSize.2_rel_hpd_width", ".pdf", sep=""),width=6, height=5)

########## reproductive number plots ####################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(R0[,(ncol(R0)-21+1):ncol(R0)], sum)/length(R0$test_0))

l <- length(R0$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(R0[,(ncol(R0)-21+1):ncol(R0)]))

p.qq_div_rate <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()+ theme(text = element_text(size=20), plot.title=element_text(size=19, hjust=0))
ggsave(plot=p.qq_div_rate,paste("../figures/qq_R0", ".pdf", sep=""),width=5, height=5)

p.qq_div_rate <- p.qq_div_rate+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Reproductive number")

ggsave(plot=p.qq_div_rate,paste("../figures/qq_R0_boot", ".pdf", sep=""),width=5, height=5)

x.min_R0 = min(R0$true)
x.max_R0 = max(R0$true)

y.min_R0 = min(R0$lower)
y.max_R0 = max(R0$upper)


lim.min = min(x.min_R0, y.min_R0)
lim.max = max(x.max_R0, y.max_R0)

p.R0 <- ggplot(R0)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))



ggsave(plot=p.R0,paste("../figures/R0", ".pdf", sep=""),width=6, height=5)

p.R0_log <- ggplot(R0)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.R0_log,paste("../figures/R0_logScale", ".pdf", sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('R0'),
                      "HPD coverage"=c(mean(R0$test)))
write.csv(hpd.test, file = paste("../figures/R0_HPD_test",".csv", sep="" ))

hpd.width = data.frame("Parameter"=c('R0'),
                       "Average relative HPD width"=c(mean(R0$hpd.rel.width)))
write.csv(hpd.width, file = paste("../figures/R0_HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('R0'),
                "Average Coefficient of Variation"=c(mean(R0$cv)))
write.csv(cv, file = paste("../figures/R0_CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('R0'),
                     "Medians"=c(mean(R0$rel.err.meadian)))
write.csv(medians, file = paste("../figures/R0_rel_error_median",".csv", sep="" ))

means = data.frame("Parameter=c('R0')", "Means"=c(mean(R0$rel.err.mean)))
write.csv(means, file = paste("../figures/R0_rel_error_mean",".csv", sep=""))


dt_rel.err.median.R0 = data.frame(param='R0', value=R0$rel.err.meadian)
p_rel_err_R0<- ggplot(dt_rel.err.median.R0, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.R0$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_R0,paste("../figures/R0_rel_err_median", ".pdf", sep=""),width=6, height=5)

dt_rel.err.mean.R0 = data.frame(param='R0', value=R0$rel.err.mean)
p_rel_err_mean_R0<- ggplot(dt_rel.err.mean.R0, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.R0$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_mean_R0,paste("../figures/R0_rel_err_mean", ".pdf", sep=""),width=6, height=5)

dt_hpd.rel.width.R0 = data.frame(param='R0', value=R0$hpd.rel.width)
p_hpd.rel.width.R0 <- ggplot(dt_hpd.rel.width.R0, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.R0$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_hpd.rel.width.R0,paste("../figures/R0_rel_hpd_width", ".pdf", sep=""),width=6, height=5)

########## become uninfectious rate plots ####################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(delta[,(ncol(delta)-21+1):ncol(delta)], sum)/length(delta$test_0))

l <- length(delta$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(delta[,(ncol(delta)-21+1):ncol(delta)]))

p.qq_div_rate <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()+ theme(text = element_text(size=20), plot.title=element_text(size=19, hjust=0))
ggsave(plot=p.qq_div_rate,paste("../figures/qq_delta", ".pdf", sep=""),width=5, height=5)

p.qq_div_rate <- p.qq_div_rate+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Become uninfectious rate")

ggsave(plot=p.qq_div_rate,paste("../figures/qq_delta_boot", ".pdf", sep=""),width=5, height=5)

x.min_delta = min(delta$true)
x.max_delta = max(delta$true)

y.min_delta = min(delta$lower)
y.max_delta = max(delta$upper)


lim.min = min(x.min_delta, y.min_delta)
lim.max = max(x.max_delta, y.max_delta)

p.delta <- ggplot(delta)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))



ggsave(plot=p.delta,paste("../figures/delta", ".pdf", sep=""),width=6, height=5)

p.delta_log <- ggplot(delta)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.delta_log,paste("../figures/delta_logScale", ".pdf", sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('delta'),
                      "HPD coverage"=c(mean(delta$test)))
write.csv(hpd.test, file = paste("../figures/delta_HPD_test",".csv", sep="" ))

hpd.width = data.frame("Parameter"=c('delta'),
                       "Average relative HPD width"=c(mean(delta$hpd.rel.width)))
write.csv(hpd.width, file = paste("../figures/delta_HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('delta'),
                "Average Coefficient of Variation"=c(mean(delta$cv)))
write.csv(cv, file = paste("../figures/delta_CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('delta'),
                     "Medians"=c(mean(delta$rel.err.meadian)))
write.csv(medians, file = paste("../figures/delta_rel_error_median",".csv", sep="" ))

means = data.frame("Parameter=c('delta')", "Means"=c(mean(delta$rel.err.mean)))
write.csv(means, file = paste("../figures/delta_rel_error_mean",".csv", sep=""))


dt_rel.err.median.delta = data.frame(param='delta', value=delta$rel.err.meadian)
p_rel_err_delta<- ggplot(dt_rel.err.median.delta, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.delta$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_delta,paste("../figures/delta_rel_err_median", ".pdf", sep=""),width=6, height=5)

dt_rel.err.mean.delta = data.frame(param='delta', value=delta$rel.err.mean)
p_rel_err_mean_delta<- ggplot(dt_rel.err.mean.delta, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.delta$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_mean_delta,paste("../figures/delta_rel_err_mean", ".pdf", sep=""),width=6, height=5)

dt_hpd.rel.width.delta = data.frame(param='delta', value=delta$hpd.rel.width)
p_hpd.rel.width.delta <- ggplot(dt_hpd.rel.width.delta, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.delta$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_hpd.rel.width.delta,paste("../figures/delta_rel_hpd_width", ".pdf", sep=""),width=6, height=5)


########## removal prob plots ####################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(removalProb[,(ncol(removalProb)-21+1):ncol(removalProb)], sum)/length(removalProb$test_0))

l <- length(removalProb$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(removalProb[,(ncol(removalProb)-21+1):ncol(removalProb)]))

p.qq_div_rate <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()+ theme(text = element_text(size=20), plot.title=element_text(size=19, hjust=0))
ggsave(plot=p.qq_div_rate,paste("../figures/qq_r", ".pdf", sep=""),width=5, height=5)

p.qq_div_rate <- p.qq_div_rate+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Removal probability")

ggsave(plot=p.qq_div_rate,paste("../figures/qq_r_boot", ".pdf", sep=""),width=5, height=5)

x.min_removalProb = min(removalProb$true)
x.max_removalProb = max(removalProb$true)

y.min_removalProb = min(removalProb$lower)
y.max_removalProb = max(removalProb$upper)


lim.min = min(x.min_removalProb, y.min_removalProb)
lim.max = max(x.max_removalProb, y.max_removalProb)

p.removalProb <- ggplot(removalProb)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))



ggsave(plot=p.removalProb,paste("../figures/removalProb", ".pdf", sep=""),width=6, height=5)

p.removalProb_log <- ggplot(removalProb)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.removalProb_log,paste("../figures/removalProb_logScale", ".pdf", sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('removalProb'),
                      "HPD coverage"=c(mean(removalProb$test)))
write.csv(hpd.test, file = paste("../figures/removalProb_HPD_test",".csv", sep="" ))

hpd.width = data.frame("Parameter"=c('removalProb'),
                       "Average relative HPD width"=c(mean(removalProb$hpd.rel.width)))
write.csv(hpd.width, file = paste("../figures/removalProb_HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('removalProb'),
                "Average Coefficient of Variation"=c(mean(removalProb$cv)))
write.csv(cv, file = paste("../figures/removalProb_CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('removalProb'),
                     "Medians"=c(mean(removalProb$rel.err.meadian)))
write.csv(medians, file = paste("../figures/removalProb_rel_error_median",".csv", sep="" ))

means = data.frame("Parameter=c('removalProb')", "Means"=c(mean(removalProb$rel.err.mean)))
write.csv(means, file = paste("../figures/removalProb_rel_error_mean",".csv", sep=""))


dt_rel.err.median.removalProb = data.frame(param='removalProb', value=removalProb$rel.err.meadian)
p_rel_err_removalProb<- ggplot(dt_rel.err.median.removalProb, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.removalProb$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_removalProb,paste("../figures/removalProb_rel_err_median", ".pdf", sep=""),width=6, height=5)

dt_rel.err.mean.removalProb = data.frame(param='removalProb', value=removalProb$rel.err.mean)
p_rel_err_mean_removalProb<- ggplot(dt_rel.err.mean.removalProb, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.removalProb$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_mean_removalProb,paste("../figures/removalProb_rel_err_mean", ".pdf", sep=""),width=6, height=5)

dt_hpd.rel.width.removalProb = data.frame(param='removalProb', value=removalProb$hpd.rel.width)
p_hpd.rel.width.removalProb <- ggplot(dt_hpd.rel.width.removalProb, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.removalProb$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_hpd.rel.width.removalProb,paste("../figures/removalProb_rel_hpd_width", ".pdf", sep=""),width=6, height=5)


########## sampling prob at present plots ####################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(rsamplingAtPresentProb[,(ncol(rsamplingAtPresentProb)-21+1):ncol(rsamplingAtPresentProb)], sum)/length(rsamplingAtPresentProb$test_0))

l <- length(rsamplingAtPresentProb$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(rsamplingAtPresentProb[,(ncol(rsamplingAtPresentProb)-21+1):ncol(rsamplingAtPresentProb)]))

p.qq_div_rate <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()+ theme(text = element_text(size=20), plot.title=element_text(size=19, hjust=0))
ggsave(plot=p.qq_div_rate,paste("../figures/qq_rho", ".pdf", sep=""),width=5, height=5)

p.qq_div_rate <- p.qq_div_rate+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Sampling at present prob.")

ggsave(plot=p.qq_div_rate,paste("../figures/qq_rho_boot", ".pdf", sep=""),width=5, height=5)

x.min_rsamplingAtPresentProb = min(rsamplingAtPresentProb$true)
x.max_rsamplingAtPresentProb = max(rsamplingAtPresentProb$true)

y.min_rsamplingAtPresentProb = min(rsamplingAtPresentProb$lower)
y.max_rsamplingAtPresentProb = max(rsamplingAtPresentProb$upper)


lim.min = min(x.min_rsamplingAtPresentProb, y.min_rsamplingAtPresentProb)
lim.max = max(x.max_rsamplingAtPresentProb, y.max_rsamplingAtPresentProb)

p.rsamplingAtPresentProb <- ggplot(rsamplingAtPresentProb)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))



ggsave(plot=p.rsamplingAtPresentProb,paste("../figures/rsamplingAtPresentProb", ".pdf", sep=""),width=6, height=5)

p.rsamplingAtPresentProb_log <- ggplot(rsamplingAtPresentProb)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.rsamplingAtPresentProb_log,paste("../figures/rsamplingAtPresentProb_logScale", ".pdf", sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('rsamplingAtPresentProb'),
                      "HPD coverage"=c(mean(rsamplingAtPresentProb$test)))
write.csv(hpd.test, file = paste("../figures/rsamplingAtPresentProb_HPD_test",".csv", sep="" ))

hpd.width = data.frame("Parameter"=c('rsamplingAtPresentProb'),
                       "Average relative HPD width"=c(mean(rsamplingAtPresentProb$hpd.rel.width)))
write.csv(hpd.width, file = paste("../figures/rsamplingAtPresentProb_HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('rsamplingAtPresentProb'),
                "Average Coefficient of Variation"=c(mean(rsamplingAtPresentProb$cv)))
write.csv(cv, file = paste("../figures/rsamplingAtPresentProb_CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('rsamplingAtPresentProb'),
                     "Medians"=c(mean(rsamplingAtPresentProb$rel.err.meadian)))
write.csv(medians, file = paste("../figures/rsamplingAtPresentProb_rel_error_median",".csv", sep="" ))

means = data.frame("Parameter=c('rsamplingAtPresentProb')", "Means"=c(mean(rsamplingAtPresentProb$rel.err.mean)))
write.csv(means, file = paste("../figures/rsamplingAtPresentProb_rel_error_mean",".csv", sep=""))


dt_rel.err.median.rsamplingAtPresentProb = data.frame(param='rsamplingAtPresentProb', value=rsamplingAtPresentProb$rel.err.meadian)
p_rel_err_rsamplingAtPresentProb<- ggplot(dt_rel.err.median.rsamplingAtPresentProb, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.rsamplingAtPresentProb$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_rsamplingAtPresentProb,paste("../figures/rsamplingAtPresentProb_rel_err_median", ".pdf", sep=""),width=6, height=5)

dt_rel.err.mean.rsamplingAtPresentProb = data.frame(param='rsamplingAtPresentProb', value=rsamplingAtPresentProb$rel.err.mean)
p_rel_err_mean_rsamplingAtPresentProb<- ggplot(dt_rel.err.mean.rsamplingAtPresentProb, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.rsamplingAtPresentProb$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_rel_err_mean_rsamplingAtPresentProb,paste("../figures/rsamplingAtPresentProb_rel_err_mean", ".pdf", sep=""),width=6, height=5)

dt_hpd.rel.width.rsamplingAtPresentProb = data.frame(param='rsamplingAtPresentProb', value=rsamplingAtPresentProb$hpd.rel.width)
p_hpd.rel.width.rsamplingAtPresentProb <- ggplot(dt_hpd.rel.width.rsamplingAtPresentProb, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.rsamplingAtPresentProb$value)$stats[c(1, 5)]*1.05)
ggsave(plot=p_hpd.rel.width.rsamplingAtPresentProb,paste("../figures/rsamplingAtPresentProb_rel_hpd_width", ".pdf", sep=""),width=6, height=5)

save.image(file="plots.RData")








