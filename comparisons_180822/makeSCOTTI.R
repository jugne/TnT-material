rm(list = ls())

library(ape)
library(seqinr)
library(beastio)
library(phytools)

'%notin%' <- Negate('%in%')

wd<-"~/Documents/Source/TnT-material/comparisons_180822/SCOTTI"
setwd(wd)


tnt_path<-"../tnt/"
# xml generator from SCOTTI repository, version 2.0.2
SCOTTI_xml_generator_path<-"~/Documents/Source/SCOTTI/scripts/SCOTTI_generate_xml.py"

n_gene_samples<-5

rates<-read.csv(paste0(tnt_path,"true_rates.csv"))

for (ii in rates$X){
  setwd(wd)
  rundir<-paste0("run_",ii)
  dir.create(rundir)
  
  alignment<-read.nexus.data(paste0(tnt_path,rundir,"/sim/sampledGeneTreeSim.alignment.nexus"))
  write.fasta(alignment, names(alignment), file.out = paste0(rundir,"/alignment.fasta"))
  
  tree_f <- list.files(path=paste0(tnt_path,rundir, "/sim/"), pattern="trees_and_pars.txt", full.names = TRUE)
  tp <- readLines(tree_f)
  tree_full <- tp[2]
  tree_sampled <- paste0(tp[4], ";")
  tree<- read.tree(text=tree_sampled)
  offset <- as.numeric(tp[which(tp=="offset")+1])
  line_first_trait <- which(tp=="sampled tree traits")+1
  line_last_trait <- which(tp=="same patient samples")-1
  n_hosts<-as.numeric(tp[which(tp=="number of hosts")+1])
  
  n_host_samples <- length(tree$tip.label)
  
  dates<-data.frame(tip=character(n_host_samples*n_gene_samples), date=numeric(n_host_samples*n_gene_samples),stringsAsFactors = FALSE)
  hosts<-data.frame(sample=character(n_host_samples*n_gene_samples), host=character(n_host_samples*n_gene_samples), stringsAsFactors = FALSE)
  hostTimes<-data.frame(host=character(n_host_samples), startTime=numeric(n_host_samples), endTime=numeric(n_host_samples), stringsAsFactors = FALSE)
  i=1
  h=1
  for(k in line_first_trait:line_last_trait){
    for (j in 1:n_gene_samples){
      dates$tip[(i-1)*n_gene_samples+j]=paste0(str_split(tp[k], "=")[[1]][1],j)
      dates$date[(i-1)*n_gene_samples+j]=3.0-as.numeric(str_split(str_split(tp[k], "=")[[1]][2], ",")[[1]][1])+offset
      
      hosts$sample[(i-1)*n_gene_samples+j]=paste0(str_split(tp[k], "=")[[1]][1],j)
      hosts$host[(i-1)*n_gene_samples+j]=str_split(str_split(tp[k], "=")[[1]][1], "_")[[1]][1]
    }
    if(str_split(str_split(tp[k], "=")[[1]][1], "_")[[1]][1] %notin% hostTimes$host){
      hostTimes$host[h]=str_split(str_split(tp[k], "=")[[1]][1], "_")[[1]][1]
      hostTimes$startTime[h]=-3.
      hostTimes$endTime[h]=3.0+offset+0.00001#-as.numeric(str_split(str_split(tp[k], "=")[[1]][2], ",")[[1]][1])+offset+0.00001
      h=h+1
    } 

    i=i+1
  }
  
  if (length(which(hostTimes$host==""))!=0){
    hostTimes <- hostTimes[-which(hostTimes$host==""),]
  }
  setwd(paste0(wd,"/",rundir))
  write.table(dates, file="dates.csv", sep=",", row.names = FALSE, col.names = FALSE, quote=FALSE)
  write.table(hosts, file="hosts.csv", sep=",", row.names = FALSE, col.names = FALSE, quote=FALSE)
  write.table(hostTimes, file="host_exposure_times.csv", sep=",", row.names = FALSE, col.names = FALSE, quote=FALSE)
  
  
  
  cmd<-paste0('python2 ',SCOTTI_xml_generator_path,' --overwrite --fasta alignment.fasta --dates dates.csv --hosts ', 
              'hosts.csv --hostTimes host_exposure_times.csv --output scotti_',rundir,' --mutationModel "JC" --numIter 10000000 --maxHosts=1000')
  system(cmd)
  
  
  xml <- readLines(paste0('scotti_',rundir,'.xml'))
  xml  <- gsub(pattern = '<substModel spec="jc69">',
               replace = '<substModel spec="JukesCantor">', x = xml)
  xml  <- gsub(pattern = "<frequencies estimate=\"false\" spec='Frequencies'>",
               replace = '', x = xml)
  xml  <- gsub(pattern = '<frequencies spec=\'RealParameter\' id="JC.freq" value="0.25 0.25 0.25 0.25"/>',
               replace = '', x = xml)
  xml  <- gsub(pattern = '</frequencies>',
               replace = '', x = xml)
  xml  <- gsub(pattern = '<stateNode idref="JC.freq"/>',
               replace = '', x = xml)
  xml  <- gsub(pattern = '<operator spec="DeltaExchangeOperator" id="freqExchanger" parameter="@JC.freq" delta="0.01" weight="0.1"/>',
               replace = '', x = xml)
  xml  <- gsub(pattern = '<log idref="JC.freq"/>',
               replace = '', x = xml)
  
  xml  <- gsub(pattern = '<mutationRate spec=\'RealParameter\' id="mutationRate" value="1.0"/>',
               replace = '<mutationRate spec=\'RealParameter\' id="mutationRate" value="0.005"/>', x = xml)
  xml  <- gsub(pattern = '<operator spec="ScaleOperator" id="muRateScaler" parameter="@mutationRate" scaleFactor="0.8" weight="1"/>',
               replace = '', x = xml)

  
  writeLines(xml, con=paste0('scotti_',rundir,'.xml'))
}


