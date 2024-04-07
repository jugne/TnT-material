library(stringr)
library(coda)
library(PRROC)
library(pROC)
library(superheat)
library(ggplot2)



#Path where this file is located:
comparisons_dir <- "~/Documents/Source/TnT-material/comparisons_180822/"

'%notin%' <- Negate('%in%')
nruns=100

wd<-comparisons_dir
setwd(wd)

tnt_path <- paste0(comparisons_dir, "tnt/")
scotti_path <- paste0(comparisons_dir, "SCOTTI/")

inference_file <- ""
output_file <- ""

### Tools are part of TnT repository: https://github.com/jugne/TnT/

#https://github.com/jugne/TnT/blob/master/src/tnt/mapping/ReMapTool.java
mapper_path <- "~/Documents/Source/TnT/out/artifacts/ReMap_jar/ReMapTool.jar"
#https://github.com/jugne/TnT/blob/master/src/tnt/tntAnnotator/TransmissionAnalyser.java
annotator_path <- "~/Documents/Source/TnT/out/artifacts/TransmissionAnalyser_jar/TransmissionAnalyser.jar"
# https://github.com/jugne/TnT/blob/master/scripts/summarize_transmission.py
analyser_path<- "~/Documents/Source/TnT/scripts/summarize_transmission.py"
#This is from scotti repository
scotti_analyser_path <- "~/Documents/Source/SCOTTI/scripts/Make_transmission_tree_alternative.py"

kk <- 1
used.rates<- c()
post_ess <- numeric(nruns)
post_ess2 <- numeric(nruns)
remove_rows<-c()
m <- matrix(ncol=300,nrow=200)
n <- nruns*100
n_<-1
transmissionTimes =data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n), lower=numeric(n), rel.err.median=numeric(n),
                                        hpd.rel.width=numeric(n), test=numeric(n), unobserved=numeric(n))
nHiddenEvents =data.frame(true=numeric(nruns), estimated=numeric(nruns), upper=numeric(nruns), lower=numeric(nruns), rel.err.median=numeric(nruns),
                              hpd.rel.width=numeric(nruns), test=numeric(nruns))
d_<-1
data <- data.frame(truth=numeric(10*n), tnt=numeric(10*n), scotti=numeric(10*n))
data_indirect <- data.frame(truth=numeric(10*n), tnt=numeric(10*n), scotti=numeric(10*n))
nd_<-1
data_no_direction <- data.frame(truth=numeric(10*n), tnt=numeric(10*n), scotti=numeric(10*n))
data_no_direction_indirect <- data.frame(truth=numeric(10*n), tnt=numeric(10*n), scotti=numeric(10*n))
for (i in 1:nruns){
  log <- list.files(path=paste0(tnt_path,"run_",i,"/inf"), pattern="inference.log", full.names = TRUE)
  log2 <- list.files(path=paste0(scotti_path,"run_",i), pattern=paste0("scotti_run_",i,".log"), full.names = TRUE)
  # read in the log file
  if (length(log)==0 | length(log2)==0){
    remove_rows = cbind(remove_rows, i)
    next
  }
  t <- read.table(log[[1]], header=TRUE, sep="\t")
  t2 <- read.table(log2[[1]], header=TRUE, sep="\t")
  # take a 20% burnin
  t <- t[-seq(1,ceiling(length(t$posterior)/30)), ]
  t2 <- t2[-seq(1,ceiling(length(t2$posterior)/30)), ]
  
  ess<- coda::effectiveSize(t$posterior)
  post_ess[i] = as.numeric(ess)#["posterior"])
  ess2<-  coda::effectiveSize(t2$posterior)
  post_ess2[i] = as.numeric(ess2)#["posterior"])
  
  # if(post_ess[i]<200 | post_ess2[i]<200){
  #   remove_rows = cbind(remove_rows, i)
  #   next
  # }
rm(t)
rm(t2)

dir.create(paste0(tnt_path,"run_", i, "/out"))
setwd(paste0(tnt_path,"run_", i, "/out"))

if (!file.exists("inference_map.xml")){
  cmd<- paste0("java -jar ", mapper_path, " -xml ../inf/inference.xml -tree ../inf/inference.transmission.trees -log ../inf/inference.log -out inference_map.xml")
  system(cmd)
}
if (!file.exists("inference_map_infectionTimes.txt")){
  cmd <- paste0("java -jar ", annotator_path, " -tree inference_map.trees -out inf_direct_tranmission.log -burnIn 30")
  system(cmd)
}
if (!file.exists("true_direct_sampled_transmission.csv")){
  cmd<-paste0("python3 ", analyser_path, " -i inf_direct_tranmission.log -o ./ -t ../sim/trees_and_pars.txt -tt inference_map_infectionTimes.txt")
  system(cmd)
}

df_truth_direct <- read.csv("true_direct_sampled_transmission.csv", header=T)
df_truth_indirect <- read.csv("true_indirect_sampled_transmission.csv", header=T)
df_truth_tr_times <- read.csv("true_sampled_transmission_times.csv", header=T)
df_truth <- data.frame(from=c(df_truth_direct$from, df_truth_indirect$from), to=c(df_truth_direct$to,df_truth_indirect$to))

df_tnt_tr_direct <- read.csv("inferred_direct_transmission.csv", header=T)
df_tnt_tr <- read.csv("inferred_transmission.csv", header=T)
df_tnt_tr_hosts <- read.csv("inferred_host_summary.csv", header=T)
df_tnt_tr_times <- read.table("inference_map_infectionTimes.txt", header = T)
df_tnt_hiden_events <- read.table("inference_map.stats.log", header=T)
df_tnt_hiden_events <- df_tnt_hiden_events[-seq(1,ceiling(length(df_tnt_hiden_events$Sample)/20)), ]

setwd(paste0(scotti_path, "run_", i))
if (!file.exists("tranmission_network.txt")){
  cmd<- paste0("python3 ", scotti_analyser_path, " --input scotti_run_",i,".trees --outputF tranmission")
  system(cmd)
}

log_scotti <- readLines("tranmission_network.txt")


scotti_probs_in_true_direct <- numeric(nrow(df_truth_direct))
tnt_probs_in_true_direct <- numeric(nrow(df_truth_direct))
for (j in 1:nrow(df_truth_direct)){
    d<- df_truth_direct$from[j]
    r<- df_truth_direct$to[j]

    direct<- log_scotti[which(log_scotti==paste0("To host ",r," from : "))[1]+1]
    direct <- strsplit(direct, ", ")

    root<- log_scotti[grep(paste0("Probabilities of being root: *"),log_scotti, ignore.case=TRUE)]
    root<- substring(root, 30)
    root<- strsplit(root, ", ")

    direct_idx <-grep(paste0("^",d," *"),direct[[1]], ignore.case=TRUE)
    root_idx <- grep(paste0("^",r," *"),root[[1]], ignore.case=TRUE)
    # unsampled_idx<-grep(paste0("^Unsampled *"),direct[[1]])

    if(length(direct_idx)==0){
      scotti_probs_in_true_direct[j] <- 0
    }else{
      if (d=="unsampled"){
        if (as.numeric(strsplit(root[[1]][root_idx]," ")[[1]][2]) > as.numeric(strsplit(direct[[1]][direct_idx]," ")[[1]][2])){
          scotti_probs_in_true_direct[j] <- 0
        } else {
          scotti_probs_in_true_direct[j] <- as.numeric(strsplit(direct[[1]][direct_idx]," ")[[1]][2]) - as.numeric(strsplit(root[[1]][root_idx]," ")[[1]][2])
        }
      } else{
        scotti_probs_in_true_direct[j] <- as.numeric(strsplit(direct[[1]][direct_idx]," ")[[1]][2])
      }
    }

    ## tnt
    idx <- which(df_tnt_tr_direct$from==d & df_tnt_tr_direct$to==r)
    if (length(idx)>0){
      tnt_probs_in_true_direct[j] <- as.numeric(df_tnt_tr_direct$probability[idx])
    } else{
      tnt_probs_in_true_direct[j] <- 0
    }
}

for (h in unique(df_truth_tr_times$host)){
    row_true <- df_truth_tr_times[which(df_truth_tr_times$host==h),]
    row_est <- df_tnt_tr_hosts[which(df_tnt_tr_hosts$id==h),]
    transmissionTimes$true[n_]<- row_true[,"time"]
    transmissionTimes$estimated[n_] <- row_est[,"median_infection_times"]
    transmissionTimes$lower[n_] <- row_est[,"hpd_lower_infection_time"]
    transmissionTimes$upper[n_] <- row_est[,"hpd_upper_infection_time"]
    transmissionTimes$rel.err.median[n_]<- abs(row_est[,"median_infection_times"]-row_true[,"time"])/row_true[,"time"]
    transmissionTimes$test[n_] <- c(as.numeric(transmissionTimes$true[n_] >=  transmissionTimes$lower[n_] & transmissionTimes$true[n_] <= transmissionTimes$upper[n_]))
    transmissionTimes$hpd.rel.width[n_] <-  (transmissionTimes$upper[n_]- transmissionTimes$lower[n_]) / transmissionTimes$true[n_]
    transmissionTimes$unobserved[n_] <- as.numeric((h %in% df_truth_indirect$to))
    n_=n_+1
}

nHiddenEvents$true[i] <- length(unique(unlist(str_split(df_truth_indirect$intermediate, " "))))
nHiddenEvents$estimated[i] <- median(df_tnt_hiden_events$nHiddenEvents)
nHiddenEvents$lower[i] <- HPDinterval(mcmc(df_tnt_hiden_events$nHiddenEvents))[1]
nHiddenEvents$upper[i] <- HPDinterval(mcmc(df_tnt_hiden_events$nHiddenEvents))[2]
nHiddenEvents$rel.err.median[i] <- abs(nHiddenEvents$estimated[i]-nHiddenEvents$true[i])/nHiddenEvents$true[i]
nHiddenEvents$test[i] <- c(as.numeric(nHiddenEvents$true[i] >=  nHiddenEvents$lower[i] & nHiddenEvents$true[i] <= nHiddenEvents$upper[i]))
nHiddenEvents$hpd.rel.width[i] <-  (nHiddenEvents$upper[i]- nHiddenEvents$lower[i]) / nHiddenEvents$true[i]



if (nrow(df_truth_direct) > 0){

  m[kk, 1:length(tnt_probs_in_true_direct)]<-tnt_probs_in_true_direct
  m[kk, 101:(100+length(scotti_probs_in_true_direct))]<-scotti_probs_in_true_direct
  # # m[kk, 201:(200+length(tmp_vect_scotti_indirect))]<-tmp_vect_scotti_indirect
  kk=kk+1

for (ii in 2:length(df_tnt_tr_hosts$id)){
  ## no direction
  d <- df_tnt_tr_hosts$id[ii-1]
  for (jj in ii:length(df_tnt_tr_hosts$id)){
    r_ <-  df_tnt_tr_hosts$id[jj]
    ## truth
    data_no_direction$truth[nd_] <- 0.0
    for (l in 1:length(df_truth_direct$from)){
      if ((df_truth_direct$from[l]==d & df_truth_direct$to[l]==r_) | (df_truth_direct$from[l]==r_ & df_truth_direct$to[l]==d)){
        data_no_direction$truth[nd_] <- 1.0
        break
      }
    }

    data_no_direction_indirect$truth[nd_] <- 0.0
    for (l in 1:length(df_truth$from)){
      if ((df_truth$from[l]==d & df_truth$to[l]==r_) | (df_truth$from[l]==r_ & df_truth$from[l]==d)){
        data_no_direction_indirect$truth[nd_]<-1.0
        break
      }
    }

    ## tnt
    idx <- which((df_tnt_tr_direct$from==d & df_tnt_tr_direct$to==r_) | (df_tnt_tr_direct$from==r_ & df_tnt_tr_direct$to==d))
    if (length(idx)>0){
      data_no_direction$tnt[nd_] <- sum(as.numeric(df_tnt_tr_direct$probability[idx]))
    } else{
      data_no_direction$tnt[nd_] <- 0
    }

    idx_indirect <- which((df_tnt_tr$from==d & df_tnt_tr$to==r_) | (df_tnt_tr$from==r_ & df_tnt_tr$to==d))
    if (length(idx)>0){
      data_no_direction_indirect$tnt[nd_] <- sum(as.numeric(df_tnt_tr$probability[idx_indirect]))
    } else{
      data_no_direction_indirect$tnt[nd_] <- 0
    }

    ## scotti
    direct_1<- log_scotti[which(log_scotti==paste0("To host ",r_," from : "))[1]+1]
    direct_1 <- strsplit(direct_1, ", ")

    direct_2<- log_scotti[which(log_scotti==paste0("To host ",d," from : "))[1]+1]
    direct_2 <- strsplit(direct_2, ", ")

    root<- log_scotti[grep(paste0("Probabilities of being root: *"),log_scotti, ignore.case=TRUE)]
    root<- substring(root, 30)
    root<- strsplit(root, ", ")

    direct_1_idx <- grep(paste0("^",d," *"),direct_1[[1]], ignore.case=TRUE)
    direct_2_idx <- grep(paste0("^",r_," *"),direct_2[[1]], ignore.case=TRUE)
    root_1_idx <- grep(paste0("^",r_," *"),root[[1]], ignore.case=TRUE)
    root_2_idx <- grep(paste0("^",d," *"),root[[1]], ignore.case=TRUE)

    if(length(direct_1_idx)!=0 | length(direct_2_idx)!=0){
      if (d=="unsampled"){
        if(length(direct_2_idx)!=0){
          data_no_direction$scotti[nd_] <- data_no_direction$scotti[nd_] + as.numeric(strsplit(direct_2[[1]][direct_2_idx]," ")[[1]][2])
        }

        tmp_1 <- as.numeric(strsplit(direct_1[[1]][direct_1_idx]," ")[[1]][2]) - as.numeric(strsplit(root[[1]][root_1_idx]," ")[[1]][2])
        if (tmp_1>0){
          data_no_direction$scotti[nd_] <- data_no_direction$scotti[nd_] + tmp_1
        }
      } else if (r_=="unsampled"){
        if(length(direct_1_idx)!=0){
          data_no_direction$scotti[nd_] <- data_no_direction$scotti[nd_] + as.numeric(strsplit(direct_1[[1]][direct_1_idx]," ")[[1]][2])
        }

        tmp_2 <- as.numeric(strsplit(direct_2[[1]][direct_2_idx]," ")[[1]][2]) - as.numeric(strsplit(root[[1]][root_2_idx]," ")[[1]][2])
        if (tmp_2>0){
          data_no_direction$scotti[nd_] =  data_no_direction$scotti[nd_] + tmp_2
        }
      } else {
        if(length(direct_1_idx)!=0){
          data_no_direction$scotti[nd_] <- data_no_direction$scotti[nd_] + as.numeric(strsplit(direct_1[[1]][direct_1_idx]," ")[[1]][2])
        }

        if(length(direct_2_idx)!=0){
          data_no_direction$scotti[nd_] <- data_no_direction$scotti[nd_] + as.numeric(strsplit(direct_2[[1]][direct_2_idx]," ")[[1]][2])
        }
      }
    }
    nd_=nd_+1
  }




  ## with direction
  for (r in df_tnt_tr_hosts$id){
    if (d==r){
      next
    }

    ## truth
    data$truth[d_] <- 0.0
    for (l in 1:length(df_truth_direct$from)){
      if (df_truth_direct$from[l]==d & df_truth_direct$to[l]==r){
        data$truth[d_] <- 1.0
        break
      }
    }

    data_indirect$truth[d_] <- 0.0
    for (l in 1:length(df_truth$from)){
      if (df_truth$from[l]==d & df_truth$to[l]==r){
        data_indirect$truth[d_] <- 1.0
        break
      }
    }

    ## tnt
    idx <- which(df_tnt_tr_direct$from==d & df_tnt_tr_direct$to==r)
    if (length(idx)>0){
      data$tnt[d_] <- as.numeric(df_tnt_tr_direct$probability[idx])
    } else{
      data$tnt[d_] <- 0
    }

    idx_indirect <- which(df_tnt_tr$from==d & df_tnt_tr$to==r)
    if (length(idx)>0){
      data_indirect$tnt[d_] <- as.numeric(df_tnt_tr$probability[idx_indirect])
    } else{
      data_indirect$tnt[d_] <- 0
    }

    ## scotti
    direct<- log_scotti[which(log_scotti==paste0("To host ",r," from : "))[1]+1]
    direct <- strsplit(direct, ", ")

    root<- log_scotti[grep(paste0("Probabilities of being root: *"),log_scotti, ignore.case=TRUE)]
    root<- substring(root, 30)
    root<- strsplit(root, ", ")

    direct_idx <-grep(paste0("^",d," *"),direct[[1]], ignore.case=TRUE)
    root_idx <- grep(paste0("^",r," *"),root[[1]], ignore.case=TRUE)

    if(length(direct_idx)==0){
      data$scotti[d_] <- 0
    }else{
      if (d=="unsampled"){
        if (as.numeric(strsplit(root[[1]][root_idx]," ")[[1]][2]) > as.numeric(strsplit(direct[[1]][direct_idx]," ")[[1]][2])){
          data$scotti[d_] <- 0
        } else {
          data$scotti[d_] <- as.numeric(strsplit(direct[[1]][direct_idx]," ")[[1]][2]) - as.numeric(strsplit(root[[1]][root_idx]," ")[[1]][2])
        }
      } else{
        data$scotti[d_] <- as.numeric(strsplit(direct[[1]][direct_idx]," ")[[1]][2])
      }
    }
    d_=d_+1
  }
}
}
}
save(m, file="/Users/jugne/Documents/TnT-material/comparisons_180822/tnt_scotti_matrix.Rdata")
save(data, file="/Users/jugne/Documents/TnT-material/comparisons_180822/tnt_scotti_data.Rdata")
save(data_no_direction, file="/Users/jugne/Documents/TnT-material/comparisons_180822/tnt_scotti_data_no_direction.Rdata")
# load("/Users/jugne/Documents/TnT-material/comparisons/tnt_scotti_matrix.Rdata")

nHiddenEvents <- nHiddenEvents[-remove_rows, ]

x.min_nHiddenEvents = min(nHiddenEvents$true)
x.max_nHiddenEvents = max(nHiddenEvents$true)
y.min_nHiddenEvents = min(nHiddenEvents$lower)
y.max_nHiddenEvents = max(nHiddenEvents$upper)
lim.min = min(x.min_nHiddenEvents, y.min_nHiddenEvents)
lim.max = max(x.max_nHiddenEvents, y.max_nHiddenEvents)

p.nHiddenEvents <- ggplot(nHiddenEvents)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper),
                colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() 
ggsave(plot=p.nHiddenEvents,paste("../../figures/nHiddenEvents", ".pdf",
                                  sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('nHiddenEvents'),
                      "HPD coverage"=c(mean(nHiddenEvents$test)))
write.csv(hpd.test, file = paste("../../figures/nHiddenEvents_HPD_test",".csv",
                                 sep="" ))



transmissionTimes <- transmissionTimes[1:(n_-1),]

x.min_transmissionTimes = min(transmissionTimes$true)
x.max_transmissionTimes = max(transmissionTimes$true)
y.min_transmissionTimes = min(transmissionTimes$lower)
y.max_transmissionTimes = max(transmissionTimes$upper)
lim.min = min(x.min_transmissionTimes, y.min_transmissionTimes)
lim.max = max(x.max_transmissionTimes, y.max_transmissionTimes)

p.transmissionTimes <- ggplot(transmissionTimes)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated,colour=unobserved), size=2, alpha=0.6) +
  theme_minimal() +
  xlim(c(0,3.1)) +
  ylim(c(0,3.1))

ggsave(plot=p.transmissionTimes,paste0(comparisons_dir,
                                       "figures/transmissionTimes_legend.pdf"),width=6, height=5)
p.transmissionTimes <- p.transmissionTimes+theme(legend.position = "none")
ggsave(plot=p.transmissionTimes,paste0(comparisons_dir,
                                       "figures/transmissionTimes.pdf"),width=6, height=5)

id <- which(transmissionTimes$unobserved==1)
tr_un <- transmissionTimes[id,]
tr_ob <- transmissionTimes[-id,]

hpd.test = data.frame("State"=c('full', 'unobserved', 'observed'),
                      "HPD coverage"=c(mean(transmissionTimes$test, na.rm=T), 
                                       mean(tr_un$test, na.rm=T), 
                                       mean(tr_ob$test, na.rm=T)),
                      "Relative HPD width"=c(mean((transmissionTimes$upper-transmissionTimes$lower)/transmissionTimes$true, na.rm = T),
                                             mean((tr_un$upper-tr_un$lower)/tr_un$true, na.rm = T),
                                             mean((tr_ob$upper-tr_ob$lower)/tr_ob$true, na.rm = T)))
write.csv(hpd.test, file = paste0(comparisons_dir,
                                  "figures/transmissionTimes_HPD_test.csv"))


## plots: heatmap

m<-m[rowSums(is.na(m)) != ncol(m),]
m<-m[,colSums(is.na(m)) != nrow(m)]
m<-m[sort(rowSums(is.na(m)), decreasing = T, index.return=T)$ix,]
#
m_tnt<-m[,1:(dim(m)[2]/2)]
m_scotti<-m[,((dim(m)[2]/2)+1):(dim(m)[2])]


png(paste0(comparisons_dir,"diff_heatmap.png"), width = 1000, height = 1000)
superheat(m_tnt-m_scotti, heat.lim = c(-1,1), heat.pal = c("#b35806", "white", "#542788"),
          heat.pal.values=c(0,0.5,1),
          grid.hline.col = "white", #X.text = round(as.matrix(m), 1), X.text.col = "black",
          #X.text.size=1,
          left.label = "none",
          bottom.label = "none",
          grid.vline.col = "white",
          heat.na.col = "white",

          yr = apply(m_tnt-m_scotti,1,sum,na.rm=TRUE)/apply(m_tnt-m_scotti,1, function(x) length(na.omit(x))),
          yr.axis.name = "normalised row sum",
          yr.plot.type = "bar")
dev.off()

## plots: with direction

data<-data[1:d_,]
roc_tnt <- roc.curve(scores.class0 = data$tnt, weights.class0 = data$truth, curve=TRUE)
roc_scotti <- roc.curve(scores.class0 = data$scotti, weights.class0 = data$truth, curve=TRUE)

roc<-ggplot() + geom_line(data=data.frame(roc_tnt$curve), aes(x=X1,y=X2, colour="TnT"),lwd=1.5, show.legend = T) +
  geom_line(data=data.frame(roc_scotti$curve), aes(x=X1,y=X2, colour="SCOTTI"),lwd=1.5,show.legend = T) +
  geom_abline(linetype="dashed",color="grey")+
  theme_minimal() +
  theme(legend.position=c(0.70,0.35), legend.text=element_text(size=15),
        legend.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  coord_cartesian(ylim=c(0, 1), xlim = c(0,1))+
  xlab("Specificity") +
  ylab("Sensitivity") +
  scale_x_continuous(labels=c("0.00" = "1.00", "0.25" = "0.75", ## because we label specificity not FPR
                              "0.50" = "0.50", "0.75" = "0.25", "1.00"="0.00"))+
  scale_colour_manual(name = 'Method',
                      values =c('TnT'='salmon','SCOTTI'='goldenrod'),
                      labels=c('TnT'=paste0("TnT, AUC = ", round(roc_tnt$auc,3)),
                               'SCOTTI'=paste0("SCOTTI, AUC = ", round(roc_scotti$auc,3))))

dir.create(paste0(comparisons_dir,"figures/pr_roc/"))
ggsave(plot=roc,paste0(comparisons_dir,
                       "figures/pr_roc/roc_direct.svg"),width=6, height=6)


pr_tnt <- pr.curve(scores.class0 = data$tnt, weights.class0 = data$truth, curve=TRUE)
pr_scotti <- pr.curve(scores.class0 = data$scotti, weights.class0 = data$truth, curve=TRUE)

prevalence = sum(data$truth)/length(data$truth)

pr<-ggplot() + geom_line(data=data.frame(pr_tnt$curve), aes(x=X1,y=X2, colour="TnT"),lwd=1.5, show.legend = T) +
  geom_line(data=data.frame(pr_scotti$curve), aes(x=X1,y=X2, colour="SCOTTI"),lwd=1.5,show.legend = T) +
  geom_hline(yintercept=prevalence, linetype="dashed",color="grey")+
  theme_minimal() +
  theme(legend.position=c(0.70,0.75), legend.text=element_text(size=15),
        legend.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  coord_cartesian(ylim=c(0, 1), xlim = c(0,1))+
  xlab("Recall") +
  ylab("Precision") +
  scale_colour_manual(name = 'Method',
                      values =c('TnT'='salmon','SCOTTI'='goldenrod'),
                      labels=c('TnT'=paste0("TnT, AUC = ", round(pr_tnt$auc.integral,3)),
                               'SCOTTI'=paste0("SCOTTI, AUC = ", round(pr_scotti$auc.integral,3))))

ggsave(plot=pr,paste0(comparisons_dir, "figures/pr_roc/pr_direct.svg"),
       width=6, height=6)


  roc_tnt<- roc(data$truth,data$tnt)
roc_scotti<- roc(data$truth,data$scotti)
c_tnt<-coords(roc_tnt, .5, "threshold", ret=c("sensitivity","specificity","ppv","npv"))
c_scotti<-coords(roc_scotti, .5, "threshold", ret=c("sensitivity","specificity","ppv","npv"))


## plots: directed transmission, with indirect events

data_indirect<-data_indirect[1:d_,]
roc_indirect_tnt <- roc.curve(scores.class0 = data_indirect$tnt, weights.class0 = data_indirect$truth, curve=TRUE)

roc_indirect<-ggplot() + geom_line(data=data.frame(roc_indirect_tnt$curve), aes(x=X1,y=X2, colour="TnT"),lwd=1.5, show.legend = T) +
  geom_abline(linetype="dashed",color="grey")+
  theme_minimal() +
  theme(legend.position=c(0.70,0.35), legend.text=element_text(size=15),
        legend.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  coord_cartesian(ylim=c(0, 1), xlim = c(0,1))+
  xlab("Specificity") +
  ylab("Sensitivity") +
  scale_x_continuous(labels=c("0.00" = "1.00", "0.25" = "0.75", ## because we label specificity not FPR
                              "0.50" = "0.50", "0.75" = "0.25", "1.00"="0.00"))+
  scale_colour_manual(name = 'Method',
                      values =c('TnT'='salmon'),
                      labels=c(paste0("TnT, AUC = ", round(roc_indirect_tnt$auc,3))))
ggsave(plot=roc_indirect,paste0(comparisons_dir, "figures/pr_roc/roc_indirect.pdf"),
       width=6, height=6)

pr_indirect_tnt <- pr.curve(scores.class0 = data_indirect$tnt, weights.class0 = data_indirect$truth, curve=TRUE)

prevalence_indirect = sum(data_indirect$truth)/length(data_indirect$truth)

pr_indirect<-ggplot() + geom_line(data=data.frame(pr_indirect_tnt$curve), aes(x=X1,y=X2, colour="TnT"),lwd=1.5, show.legend = T) +
  geom_hline(yintercept=prevalence_indirect, linetype="dashed",color="grey")+
  theme_minimal() +
  theme(legend.position=c(0.70,0.75), legend.text=element_text(size=15),
        legend.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  coord_cartesian(ylim=c(0, 1), xlim = c(0,1))+
  xlab("Recall") +
  ylab("Precision") +
  scale_colour_manual(name = 'Method',
                      values =c('TnT'='salmon'),
                      labels=c(paste0("TnT, AUC = ", round(pr_indirect_tnt$auc.integral,3))))
ggsave(plot=pr_indirect,paste0(comparisons_dir, "figures/pr_roc/pr_indirect.pdf"),
       width=6, height=6)


roc_indirect_tnt<- roc(data_indirect$truth,data_indirect$tnt)
c_tnt_indirect<-coords(roc_indirect_tnt, .5, "threshold", ret=c("sensitivity","specificity","ppv","npv"))

## plots: no direction

data_no_direction<-data_no_direction[1:nd_,]

roc_nd_tnt <- roc.curve(scores.class0 = data_no_direction$tnt, weights.class0 = data_no_direction$truth, curve=TRUE)
roc_nd_scotti <- roc.curve(scores.class0 = data_no_direction$scotti, weights.class0 = data_no_direction$truth, curve=TRUE)


roc_nd<-ggplot() + geom_line(data=data.frame(roc_nd_tnt$curve), aes(x=X1,y=X2, colour="TnT"),lwd=1.5, show.legend = T) +
  geom_line(data=data.frame(roc_nd_scotti$curve), aes(x=X1,y=X2, colour="SCOTTI"),lwd=1.5,show.legend = T) +
  geom_abline(linetype="dashed",color="grey")+
  theme_minimal() +
  theme(legend.position=c(0.70,0.35), legend.text=element_text(size=15),
        legend.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  coord_cartesian(ylim=c(0, 1), xlim = c(0,1))+
  xlab("Specificity") +
  ylab("Sensitivity") +
  scale_x_continuous(labels=c("0.00" = "1.00", "0.25" = "0.75", ## because we label specificity not FPR
                              "0.50" = "0.50", "0.75" = "0.25", "1.00"="0.00"))+
  scale_colour_manual(name = 'Method',
                      values =c('TnT'='salmon','SCOTTI'='goldenrod'),
                      labels=c(paste0("TnT, AUC = ", round(roc_nd_tnt$auc,3)),
                               paste0("SCOTTI, AUC = ", round(roc_nd_scotti$auc,3))))
ggsave(plot=roc_nd,paste0(comparisons_dir, "figures/pr_roc/roc_nd_direct.pdf"),
       width=6, height=6)


pr_nd_tnt <- pr.curve(scores.class0 = data_no_direction$tnt, weights.class0 = data_no_direction$truth, curve=TRUE)
pr_nd_scotti <- pr.curve(scores.class0 = data_no_direction$scotti, weights.class0 = data_no_direction$truth, curve=TRUE)
prevalence_nd = sum(data_no_direction$truth)/length(data_no_direction$truth)

pr_nd<-ggplot() + geom_line(data=data.frame(pr_nd_tnt$curve), aes(x=X1,y=X2, colour="TnT"),lwd=1.5, show.legend = T) +
  geom_line(data=data.frame(pr_nd_scotti$curve), aes(x=X1,y=X2, colour="SCOTTI"),lwd=1.5,show.legend = T) +
  geom_hline(yintercept=prevalence_nd, linetype="dashed",color="grey")+
  theme_minimal() +
  theme(legend.position=c(0.70,0.75), legend.text=element_text(size=15),
        legend.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  coord_cartesian(ylim=c(0, 1), xlim = c(0,1))+
  xlab("Recall") +
  ylab("Precision") +
  scale_colour_manual(name = 'Method',
                      values =c('TnT'='salmon','SCOTTI'='goldenrod'),
                      labels=c(paste0("TnT, AUC = ", round(pr_nd_tnt$auc.integral,3)),
                               paste0("SCOTTI, AUC = ", round(pr_nd_scotti$auc.integral,3))))
ggsave(plot=pr_nd,paste0(comparisons_dir, "figures/pr_roc/pr_nd_direct.pdf"),
       width=6, height=6)


roc_nd_tnt<- roc(data_no_direction$truth,data_no_direction$tnt)
roc_nd_scotti<- roc(data_no_direction$truth,data_no_direction$scotti)


c_nd_tnt<-coords(roc_nd_tnt, .5, "threshold", ret=c("sensitivity","specificity","ppv","npv"))
c_nd_scotti<-coords(roc_nd_scotti, .5, "threshold", ret=c("sensitivity","specificity","ppv","npv"))

## plots: no direction. with indirect events
data_no_direction_indirect<-data_no_direction_indirect[1:nd_,]

roc_nd_indirect_tnt <- roc.curve(scores.class0 = data_no_direction_indirect$tnt, weights.class0 = data_no_direction_indirect$truth, curve=TRUE)

roc_nd_indirect<-ggplot() + geom_line(data=data.frame(roc_nd_indirect_tnt$curve), aes(x=X1,y=X2, colour="TnT"),lwd=1.5, show.legend = T) +
  geom_abline(linetype="dashed",color="grey")+
  theme_minimal() +
  theme(legend.position=c(0.70,0.35), legend.text=element_text(size=15),
        legend.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  coord_cartesian(ylim=c(0, 1), xlim = c(0,1))+
  xlab("Specificity") +
  ylab("Sensitivity") +
  scale_x_continuous(labels=c("0.00" = "1.00", "0.25" = "0.75", ## because we label specificity not FPR
                              "0.50" = "0.50", "0.75" = "0.25", "1.00"="0.00"))+
  scale_colour_manual(name = 'Method',
                      values =c('TnT'='salmon'),
                      labels=c(paste0("TnT, AUC = ", round(roc_nd_indirect_tnt$auc,3))))

ggsave(plot=roc_nd_indirect,paste0(comparisons_dir, "figures/pr_roc/roc_nd_indirect.pdf"),
       width=6, height=6)


pr_nd_indirect_tnt <- pr.curve(scores.class0 = data_no_direction_indirect$tnt, weights.class0 = data_no_direction_indirect$truth, curve=TRUE)
prevalence_nd_indirect = sum(data_no_direction_indirect$truth)/length(data_no_direction_indirect$truth)

pr_nd_indirect<-ggplot() + geom_line(data=data.frame(pr_nd_indirect_tnt$curve), aes(x=X1,y=X2, colour="TnT"),lwd=1.5, show.legend = T) +
  geom_hline(yintercept=prevalence_nd_indirect, linetype="dashed",color="grey")+
  theme_minimal() +
  theme(legend.position=c(0.70,0.75), legend.text=element_text(size=15),
        legend.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  coord_cartesian(ylim=c(0, 1), xlim = c(0,1))+
  xlab("Recall") +
  ylab("Precision") +
  scale_colour_manual(name = 'Method',
                      values =c('TnT'='salmon'),
                      labels=c(paste0("TnT, AUC = ", round(pr_nd_indirect_tnt$auc.integral,3))))
ggsave(plot=pr_nd_indirect,paste0(comparisons_dir, "figures/pr_roc/pr_nd_indirect.pdf"),
       width=6, height=6)


roc_nd_indirect_tnt<- roc(data_no_direction_indirect$truth,data_no_direction_indirect$tnt)
c_nd_indirect_tnt<-coords(roc_nd_indirect_tnt, .5, "threshold", ret=c("sensitivity","specificity","ppv","npv"))



summary<-data.frame(Sensitivity=numeric(6), Specificity=numeric(6), PPV=numeric(6), NPV=numeric(6))

summary$Sensitivity <- c(c_tnt$sensitivity, c_scotti$sensitivity,
                         c_nd_tnt$sensitivity, c_nd_scotti$sensitivity,
                         c_tnt_indirect$sensitivity, c_nd_indirect_tnt$sensitivity)
summary$Specificity <- c(c_tnt$specificity, c_scotti$specificity,
                         c_nd_tnt$specificity, c_nd_scotti$specificity,
                         c_tnt_indirect$specificity, c_nd_indirect_tnt$specificity)
summary$PPV <- c(c_tnt$ppv, c_scotti$ppv, c_nd_tnt$ppv, c_nd_scotti$ppv,
                 c_tnt_indirect$ppv, c_nd_indirect_tnt$ppv)
summary$NPV <- c(c_tnt$npv, c_scotti$npv, c_nd_tnt$npv, c_nd_scotti$npv,
                 c_tnt_indirect$npv, c_nd_indirect_tnt$npv)

summary <- round(summary, 3)
rownames(summary) <- c("TnT", "SCOTTI", "TnT_no_direction", "SCOTTI_no_direction",
                       "TnT_with_indirect", "TnT_nd_with_indirect")

write.csv(summary, file = paste0(comparisons_dir,
                                  "figures/pr_roc/summary.csv"))


