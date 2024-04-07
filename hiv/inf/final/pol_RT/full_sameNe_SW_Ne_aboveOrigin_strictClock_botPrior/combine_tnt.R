rm(list = ls())
library(coda)

# You need to install beast before runnng this and change the appropriate paths
log_combiner_path <- "/Applications/BEAST\\ 2.6.7/bin/logcombiner"

# The tool can be found https://github.com/jugne/TnT/blob/master/src/tnt/tntAnnotator/TransmissionAnalyser.java
transmission_analyser_path <- "~/Documents/Source/TnT/out/artifacts/TransmissionAnalyser_jar/TransmissionAnalyser.jar"
wd <- "~/Documents/Source/TnT-material/hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/combined/"

dir.create(wd)
setwd(wd)
burnIn <- 30
logs <- list.files(path=paste0("../"), pattern="env_pol_saConstrained_s1_sameNe_run[0-9]*\\.log$", full.names = TRUE)
trees<- list.files(path=paste0("../"), pattern="env_pol_saConstrained_s1_sameNe_run[0-9]*\\.transmission.trees$", full.names = TRUE)


in_command <- paste0(" -b ",burnIn)
posteriors <- list()
i=1
for (log in logs){
  in_command <- paste0(in_command, " -log ", log)
  posteriors[[i]] <-as.mcmc(read.table(log, header=T)["posterior"])
  i=i+1
}

out_command = gsub("run0", "combined", logs[1])
out_command = gsub("..//", "", out_command)

in_command_trees <- gsub("\\.log", ".transmission.trees", in_command)
out_command_trees = gsub("run0*", "combined.transmission.trees", logs[1])
out_command_trees = gsub("\\.log", "", out_command_trees)
out_command_trees = gsub("..//", "", out_command_trees)


system(paste(log_combiner_path, in_command," -o ",out_command, sep=" "))
system(paste(log_combiner_path, in_command_trees," -o ",out_command_trees, sep=" "))
 
# system(paste("java -jar",remap_tool_path,"-xml ..//env_pol_saConstrained_s1_sameNe_run0.xml -tree",out_command_trees,"-log",out_command,"-out combined_remap.xml"))

system(paste("java -jar",transmission_analyser_path,"-tree env_pol_saConstrained_s1_sameNe_combined.transmission.trees",
             "-out env_pol_saConstrained_s1_sameNe_combined.transmission.log"))



comb_log <- read.table(out_command, header = T)
comb_post_ess <-effectiveSize(as.mcmc(comb_log$posterior))
pdf(file="comb_log_trace.pdf", width=15, height=8)
traceplot(as.mcmc(comb_log$posterior), smooth = T,
          type = "l", xlab = "logged iterations", ylab = "posterior", col="lightblue",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
legend(
  'topright', ncol = 2L, cex = 1.2,
  legend = c(
    'Chain', "chain 1","chain 2","chain 3", "combined",
    'Ess', round(effectiveSize(as.mcmc(posteriors[[1]][-c(1:round(length(posteriors[[1]])*0.30))]))), 
    round(effectiveSize(as.mcmc(posteriors[[2]][-c(1:round(length(posteriors[[2]])*0.30))]))),
    round(effectiveSize(as.mcmc(posteriors[[3]][-c(1:round(length(posteriors[[3]])*0.30))]))),
    round(comb_post_ess)
  )
)

dev.off()
