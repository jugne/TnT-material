rm(list = ls())

library(beastio)
library(ape)
library(phytools)
library(ggplot2)
library(stringr)

'%notin%' <- Negate('%in%')

f_lambda <- function(R0, delta){
  return(R0*delta)
}

f_mu <- function(delta, s, r){
  return(delta*(1-s)/(1-(1-r)*s))
}

f_psi <- function(s, delta, r){
  return(s*delta/(1-(1-r)*s))
}


wd<- "~/Documents/Source/TnT-material/comparisons_180822/tnt/"
setwd(wd)

# source: https://github.com/jugne/sampled-ancestors/blob/v2.6/src/beast/app/simulators/Simulator.java
SaSimPath <- '../SAsimulator_01112021.jar'
# source https://github.com/jugne/TnT
TnTPath<- '../TnT.jar'

n_runs = 100
n_gene_samples = 5

# the below seed draw code proved problematic when switching between machines
# now seed is just the run number
# seed = 64418
# set.seed(seed)
# seeds <- sample(0:2147483647, n_runs)


# define params of the lognormal distribution of the Ne
mean_ne <- 2;
sigma_ne <- 0.5;
mu_ne <- log(mean_ne) - sigma_ne^2/2

origin_true = 3.0;
min_samples=5;
max_samples=50;


true_rates<-data.frame(seed=numeric(n_runs), t_or=numeric(n_runs), R0=numeric(n_runs), delta=numeric(n_runs), s=numeric(n_runs), lambda=numeric(n_runs), mu=numeric(n_runs), psi=numeric(n_runs), r=numeric(n_runs), rho=numeric(n_runs), Ne=numeric(n_runs),
                       Ne_orig=numeric(n_runs), tau=numeric(n_runs), botDuration=numeric(n_runs), pairwiseCoalProbBottleneck=numeric(n_runs))

for (j in 1:n_runs){
  set.seed(j)

  
  run_dir <- paste0(wd, "run_",j)
  unlink(run_dir, recursive = TRUE)
  dir.create(run_dir)
  sim_dir <- paste0(run_dir,"/sim")
  dir.create(sim_dir)
  inf_dir <- paste0(run_dir,"/inf")
  dir.create(inf_dir)
  
    
  # MSC params
  Ne <- rlnorm(1,mu_ne, sigma_ne)
  Ne_orig <- rlnorm(1,mu_ne, sigma_ne)
  probCoalBottleneck <- runif(1, 0.6, 1.0)
  tau <- -log(1-probCoalBottleneck)
  botDuration <- tau*Ne
  
  setwd(sim_dir)
  result = 124;
  while(result==124){
    # SABD params

    l_R0<-1.
    h_R0<-4.
    l_delta<-0.1
    h_delta<-2.
    l_s<-.3
    h_s<-1.
    l_r<-0.2
    h_r<-1.
    l_rho<-0.
    h_rho<-1.

    R0<-runif(1, l_R0,h_R0)
    delta<-runif(1, l_delta, h_delta)
    s<- runif(1.0, l_s, h_s)
    r<- runif(1.0, l_r, h_r)

    lambda <- f_lambda(R0, delta)
    mu <- f_mu(delta, s, r)
    psi<- f_psi(s, delta, r)
    # rho <- runif(1, l_rho, h_rho)
    rho <- 0.


    cmd <- paste('java -jar',SaSimPath,lambda,mu,psi,r,rho,origin_true, j,min_samples, max_samples)
    result<-system(cmd, timeout=1.0)
  }
  true_rates$seed[j] =j
  true_rates$t_or[j] = origin_true
  true_rates$R0[j] = R0
  true_rates$delta[j] = delta
  true_rates$s[j] = s
  true_rates$lambda[j] = lambda
  true_rates$mu[j] = mu
  true_rates$psi[j] = psi
  true_rates$r[j] = r
  true_rates$rho[j] = rho
  true_rates$Ne[j] = Ne
  true_rates$Ne_orig[j] = Ne_orig
  true_rates$tau[j] = tau
  true_rates$botDuration[j] = botDuration
  true_rates$pairwiseCoalProbBottleneck[j] = probCoalBottleneck
  
  tree_f <- list.files(path=paste0(sim_dir), pattern="trees_and_pars", full.names = TRUE)
  tp <- readLines(tree_f)
  tree_full <- tp[2]
  tree_sampled <- tp[4]
  line_first_trait <- which(tp=="full tree traits")+1
  line_last_trait <- which(tp=="sampled tree traits")-1
  offset <- tp[which(tp=="offset")+1]
  taxon <- ""
  taxonInf <- ""
  samples <- NULL
  gene_taxon<- ""
  date<-""
  first <- T
  
  for (k in line_first_trait:line_last_trait){
    tax <- strsplit(tp[k], '=')[[1]][1]
    taxon <- paste(taxon, paste0('<taxon spec="Taxon" id="',tax,'"/>'), collapse = "\n")
    taxonInf <- paste(taxonInf, paste0('<taxon idref="',tax,'"/>'), collapse = "\n")
    date <- paste(date, tp[k], sep="")
    if (first){
      samples <- paste(samples, paste0(tax,'=',n_gene_samples), sep = "")
    } else {
      samples <- paste(samples, paste0(tax,'=',n_gene_samples), sep = ",")
    }
    gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',tax,'" spec="TaxonSet">'), collapse = "\n")
    for (kk in 1:n_gene_samples){
      gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',tax,kk,'" spec="Taxon"/>'), collapse = "\n")
    }
    gene_taxon<- paste(gene_taxon, '</taxon>', collapse = "\n")
    first<-F
  }
  
  taxonSa <- ""
  taxonSA_idref<-""
  patientTaxonSets<-""
  saRestriction_start <- which(tp=="same patient samples")+1
  saRestriction_end <- which(tp=="origin")-1
  if (saRestriction_end-saRestriction_start>=0){
    for(l in saRestriction_start:saRestriction_end){
      lineSplit<-strsplit(tp[l], ":")[[1]]
      taxonSa <- paste(taxonSa, paste0('<taxonSet spec="TaxonSet" id="taxon_', lineSplit[1], '_">'), collapse ="\n")
      patientTaxonSets<- paste(patientTaxonSets, paste0('<patientTaxonSets idref="taxon_',  lineSplit[1], '_"/>'), collapse ="\n")
      taxonSA_idref <- paste(taxonSA_idref, paste0('<taxonsets idref="taxon_',  lineSplit[1], '_"/>'), collapse="\n")
      # taxonSa <- paste(taxonSa, paste0('<taxon idref="', lineSplit[1], '_"/>'), collapse ="\n")
      for (n in 1:as.numeric(lineSplit[2])){
        taxonSa <- paste(taxonSa, paste0('<taxon idref="', lineSplit[1], "_", n, '_"/>'), collapse ="\n")
      }
      taxonSa <- paste(taxonSa, '</taxonSet>', collapse ="\n")
    }
  }
  
  psi_fullSampling<- f_psi(1, delta, r)
  #########################################################
  ################# Make simulation XML ###################
  #########################################################
  sim <- readLines( paste0(wd,"../templates/complete_gene_simulation_template.xml"))
  sim  <- gsub(pattern = "<insertTaxon/>",
               replace = taxon, x = sim)
  sim  <- gsub(pattern = "<insertTransmissionTree/>",
               replace = paste('<transmissionTreeInput spec="TreeParser"
                                newick="',tree_full,
                               '" adjustTipHeights="false" IsLabelledNewick="true"/>', ''), x = sim)
  sim  <- gsub(pattern = "<insertSampleCounts/>",
               replace = samples, x = sim)
  
  # TnT parameters
  sim  <- gsub(pattern = "<insertPopulationSizes/>",
               replace = paste('<populationSizes spec="RealParameter" value="',
                               Ne,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertPopSizeAboveOrigin/>",
               replace = paste('<popSizeAboveOrigin spec="RealParameter" value="',
                               Ne_orig,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertBottleneckStrength/>",
               replace = paste('<bottleneckStrength spec="RealParameter" value="',
                               tau*Ne,'"/>', sep=''), x = sim)
  
  # fBD parameters
  sim  <- gsub(pattern = "<insertBirthRate/>",
               replace = paste('<birthRate spec="RealParameter" value="',
                               lambda,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertDeathRate/>",
               replace = paste('<deathRate spec="RealParameter" value="',
                               mu,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertSamplingRate/>",
               replace = paste('<samplingRate spec="RealParameter" value="',
                               psi_fullSampling,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertSamplingExtantRate/>",
               replace = paste('<samplingExtantRate spec="RealParameter" value="',
                               1,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertOrigin/>",
               replace = paste('<origin spec="RealParameter" value="',
                               origin_true,'"/>', sep=''), x = sim)
  
  writeLines(sim, con=paste0(sim_dir, '/completeGeneTreeSim.xml'))
  
  cmd <- paste0('java -jar ', TnTPath, ' -seed ', j, ' -overwrite ', sim_dir, "/completeGeneTreeSim.xml")
  # cmd <- paste0('java -jar /Users/jugne/Documents/Source/tnt_21092021.jar -seed ', seeds[j], ' -overwrite ', sim_dir, "/completeGeneTreeSim.xml")
  setwd(sim_dir)
  result<-system(cmd)
  
  n_samples<- as.numeric(tp[which(tp=="total sample count")+1])
  n_hosts <- as.numeric(tp[which(tp=="number of hosts")+1])
  line_first_trait <- which(tp=="sampled tree traits")+1
  line_last_trait <- which(tp=="same patient samples")-1
  sampled_tax_list<- character(n_samples*n_gene_samples)
  taxon <- ""
  taxonInf <- ""
  samples <- NULL
  gene_taxon<- ""
  date<-""
  first <- T
  i=1
  for (k in line_first_trait:line_last_trait){
    tax <- strsplit(tp[k], '=')[[1]][1]
    taxon <- paste(taxon, paste0('<taxon spec="Taxon" id="',tax,'"/>'), collapse = "\n")
    taxonInf <- paste(taxonInf, paste0('<taxon idref="',tax,'"/>'), collapse = "\n")
    date <- paste(date, tp[k], sep="")
    if (first){
      samples <- paste(samples, paste0(tax,'=',n_gene_samples), sep = "")
    } else {
      samples <- paste(samples, paste0(tax,'=',n_gene_samples), sep = ",")
    }
    gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',tax,'" spec="TaxonSet">'), collapse = "\n")
    for (kk in 1:n_gene_samples){
      gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',tax,kk,'" spec="Taxon"/>'), collapse = "\n")
      sampled_tax_list[(i-1)*n_gene_samples+kk]<-paste0(tax, kk)
    }
    gene_taxon<- paste(gene_taxon, '</taxon>', collapse = "\n")
    first<-F
    i=i+1
  }
  
  

  full_nexus <- list.files(path=paste0(sim_dir), pattern="completeGeneTreeSim.alignment.nexus", full.names = TRUE)
  tp <- readLines(full_nexus)
  tp[4] <- paste0("\tdimensions ntax=",n_samples*n_gene_samples,";")
  tp[5] <- paste0("\ttaxlabels ",paste(sampled_tax_list, collapse= " "), ";")
  remove_lines<-(n_hosts-n_samples)*n_gene_samples
  i = 1
  for (line in 12:length(tp)){
    tmp_substr<- strsplit(tp[line], " ")[[1]][1] 
     if(substr(tmp_substr,3,nchar(tmp_substr)) %notin% sampled_tax_list){
       remove_lines[i]<- line
       i=i+1
     }
  }
  tp<- tp[-remove_lines]
  if (str_sub(tp[length(tp)], -1, -1)!=";"){
    tp[length(tp)]<- paste0(tp[length(tp)], ";")
  }
  writeLines(tp, con=paste0(sim_dir, '/sampledGeneTreeSim.alignment.nexus'))
  
  
  
  
  #########################################################
  ################# Make inference XML ####################
  #########################################################
  
  
  # MSC params
  Ne_init <- rlnorm(1,mu_ne, sigma_ne)
  probCoalBottleneck_init <- runif(1, 0.0, 1.0)
  tau_init <- -log(1-probCoalBottleneck_init)
  botDuration_init <- tau_init*Ne_init
  
  # SABD params
  R0_init <- runif(1, l_R0, h_R0)
  delta_init <- runif(1, l_delta, h_delta)
  lambda_init <- f_lambda(R0_init, delta_init)
  r_init <- runif(1, l_r, h_r)
  mu_init <- f_mu(delta_init, s, r_init)
  # lambda_init <- runif(1, l_lambda, h_lambda)
  # mu_init <- runif(1, l_mu, h_mu)
  # rho_init <- runif(1, l_rho, h_rho)
  rho_init <- 0
  
  
  # 
  # taxon <- ""
  # samples <- NULL
  # first <- T
  # date<-""
  # gene_taxon<-""
  # taxonInf <- ""
  
  a <- read.nexus.data(paste0(sim_dir,"/sampledGeneTreeSim.alignment.nexus"))
  print(as.numeric(unique(sapply(strsplit(names(a), "_"), "[", 1 ))))
  
  # for (l in as.numeric(unique(sapply(strsplit(names(a), "_"), "[", 1 )))){
  #   taxon <- paste(taxon, paste0('<taxon spec="Taxon" id="',l,'_"/>'), collapse = "\n")
  #   taxonInf <- paste(taxonInf, paste0('<taxon idref="',l,'_"/>'), collapse = "\n")
  #   if (first){
  #     samples <- paste(samples, paste0(l,'_=',n_gene_samples), sep = "")
  #     date <- paste(date, paste0(l,'_=',nodeheight(tree, node=1)-nodeheight(tree, node=which(tree$tip.label==l))) , sep="")
  #     first <- F
  #   } else {
  #     samples <- paste(samples, paste0(l,'_=',n_gene_samples), sep = ",")
  #     date <- paste(date, paste0(l,'_=',nodeheight(tree, node=1)-nodeheight(tree, node=which(tree$tip.label==l))) , sep=",")
  #   }
  #   gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',l,'_" spec="TaxonSet">'), collapse = "\n")
  #   for (kk in 1:n_gene_samples){
  #     gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',l,"_",kk,'" spec="Taxon"/>'), collapse = "\n")
  #   }
  #   gene_taxon<- paste(gene_taxon, '</taxon>', collapse = "\n")
  #   
  # }
  
  if(saRestriction_end-saRestriction_start<0){
  inf <- readLines(paste0(wd,"../templates/inference_template.xml"))
  inf  <- gsub(pattern = "<trait>",
               replace = paste0("<trait id='tipDates' spec='beast.evolution.tree.TraitSet' traitname='date-backward' units='year' value='",date,"'>"), x = inf)
  inf  <- gsub(pattern = "<insertTaxonSuperSet/>",
               replace = gene_taxon, x = inf)
  inf  <- gsub(pattern = "<insertTaxon/>",
               replace = taxonInf, x = inf)
  
  inf  <- gsub(pattern = "<insertOrigin/>",
               replace = paste0('<origin id="origin" spec="RealParameter" value="',n_samples,'." lower="0."/>'), x = inf)
  inf  <- gsub(pattern = "<insertR0/>",
               replace = paste('<skylineValues id="R0" spec="RealParameter" value="',
                               R0_init,'"/>', sep=''), x = inf)
  inf  <- gsub(pattern = "<insertBecominUninfectiousRate/>",
               replace = paste('<skylineValues id="becominUninfectiousRate" spec="RealParameter" value="',
                               delta_init,'"/>', sep=''), x = inf)
  inf  <- gsub(pattern = "<insertsamplingProportion/>",
               replace = paste('<skylineValues id="samplingProportion" spec="RealParameter" value="',
                               s,'"/>', sep=''), x = inf)
  inf  <- gsub(pattern = "<insertRemovalProb/>",
               replace = paste('<skylineValues id="removalProb" spec="RealParameter" value="',
                               r_init,'"/>', sep=''), x = inf)
  inf  <- gsub(pattern = "<insertRhoValues/>",
               replace = paste('<values id="rho" spec="RealParameter" value="',
                               rho_init,'"/>', sep=''), x = inf)
  
  inf  <- gsub(pattern = "<insertFinalSampleOffset/>",
               replace = paste('<finalSampleOffset id="finalSampleOffset" spec="RealParameter" value="',
                               offset,'"/>', sep=''), x = inf)
  
  inf  <- gsub(pattern = "<insertSampleCounts/>",
               replace = samples, x = inf)
  inf  <- gsub(pattern = "<insertBirthRate/>",
               replace = paste('<birthRate spec="RealParameter" value="',
                               lambda_init,'"/>', sep=''), x = inf)
  inf  <- gsub(pattern = "<insertDeathRate/>",
               replace = paste('<deathRate spec="RealParameter" value="',
                               mu_init,'"/>', sep=''), x = inf)
  inf  <- gsub(pattern = "<insertSamplingRate/>",
               replace = paste('<samplingRate spec="RealParameter" value="',
                               psi,'"/>', sep=''), x = inf)
  inf<- gsub(pattern = "<insertBottleneckStrength/>",
             replace = paste('<bottleneckStrength spec="RealParameter" value="',
                             botDuration_init ,'"/>', sep=''), x = inf)
  # TnT parameters
  inf  <- gsub(pattern = "<insertPopulationSizes/>",
               replace = paste('<parameter id="popSize" name="stateNode" value="',
                               Ne_init,'"/>', sep=''), x = inf)
  inf  <- gsub(pattern = "<insertPairwiseCoalProbBottleneck/>",
               replace = paste('<parameter id="pairwiseCoalProbBottleneck" name="stateNode" value="',
                               probCoalBottleneck_init,'"/>', sep=''), x = inf)
  
  # SABD Priors
  inf  <- gsub(pattern = "<insertR0Prior/>",
               replace = paste0('<distr spec="beast.math.distributions.Uniform"
                          lower="',l_R0,'" upper="',h_R0,'" offset="0."/>'), x = inf)
  inf  <- gsub(pattern = "<insertBecominUninfectiousPrior/>",
               replace = paste0('<distr spec="beast.math.distributions.Uniform"
                          lower="',l_delta,'" upper="',h_delta,'" offset="0."/>'), x = inf)
  inf  <- gsub(pattern = "<insertRemovalPrior/>",
               replace = paste0('<distr spec="beast.math.distributions.Uniform"
                          lower="',l_r,'" upper="',h_r,'" offset="0."/>'), x = inf)
  inf  <- gsub(pattern = "<insertOriginPrior/>",
               replace = paste0('<distr spec="beast.math.distributions.Uniform"
                          lower="0." upper="1000." offset="0."/>'), x = inf)
  inf  <- gsub(pattern = "<insertRhoPrior/>",
               replace = paste0('<distr spec="beast.math.distributions.Uniform"
                          lower="',l_rho,'" upper="',h_rho,'" offset="0."/>'), x = inf)
  
  # TnT Priors
  inf  <- gsub(pattern = "<insertBottleneckPrior/>",
               replace = '<distr spec="beast.math.distributions.Uniform"
                          lower="0.0" upper="1.0" offset="0."/>', x = inf)
  inf  <- gsub(pattern = "<insertPopSizesPrior/>",
               replace = paste0('<distr spec="beast.math.distributions.LogNormalDistributionModel" M="',mean_ne,'" S="',sigma_ne,'" meanInRealSpace="true"/>'), x = inf)
  
  
  writeLines(inf, con=paste0(inf_dir, '/inference.xml'))
  }
  #########################################################
  ############ Make restricted inference XML ##############
  #########################################################
  
  if(saRestriction_end-saRestriction_start>=0){
    inf <- readLines(paste0(wd,"../templates/inference_template_saRestricted.xml"))
    inf  <- gsub(pattern = "<trait>",
                 replace = paste0("<trait id='tipDates' spec='beast.evolution.tree.TraitSet' traitname='date-backward' units='year' value='",date,"'>"), x = inf)
    inf  <- gsub(pattern = "<insertTaxonSuperSet/>",
                 replace = gene_taxon, x = inf)
    inf  <- gsub(pattern = "<insertTaxon/>",
                 replace = taxonInf, x = inf)
  
    inf  <- gsub(pattern = "<insertSARestriction/>",
                 replace = taxonSa, x = inf)
  
    inf  <- gsub(pattern = "<insertTaxonRestrictionsInit/>",
                 replace = patientTaxonSets, x = inf)
  
    inf  <- gsub(pattern = "<insertSamplingConstraintDistribution/>",
                 replace = taxonSA_idref, x = inf)
  
    inf  <- gsub(pattern = "<insertOrigin/>",
                 replace = paste0('<origin id="origin" spec="RealParameter" value="',n_samples,'." lower="0."/>'), x = inf)
    inf  <- gsub(pattern = "<insertR0/>",
                 replace = paste('<skylineValues id="R0" spec="RealParameter" value="',
                                 R0_init,'"/>', sep=''), x = inf)
    inf  <- gsub(pattern = "<insertBecominUninfectiousRate/>",
                 replace = paste('<skylineValues id="becominUninfectiousRate" spec="RealParameter" value="',
                                 delta_init,'"/>', sep=''), x = inf)
    inf  <- gsub(pattern = "<insertsamplingProportion/>",
                 replace = paste('<skylineValues id="samplingProportion" spec="RealParameter" value="',
                                 s,'"/>', sep=''), x = inf)
    inf  <- gsub(pattern = "<insertremovalProb/>",
                 replace = paste('<skylineValues id="removalProb" spec="RealParameter" value="',
                                 r_init,'"/>', sep=''), x = inf)
    inf  <- gsub(pattern = "<insertRhoValues/>",
                 replace = paste('<values id="rho" spec="RealParameter" value="',
                                 0,'"/>', sep=''), x = inf)
  
    inf  <- gsub(pattern = "<insertFinalSampleOffset/>",
                 replace = paste('<finalSampleOffset id="finalSampleOffset" spec="RealParameter" value="',
                                 offset,'"/>', sep=''), x = inf)
  
    inf  <- gsub(pattern = "<insertSampleCounts/>",
                 replace = samples, x = inf)
    inf  <- gsub(pattern = "<insertBirthRate/>",
                 replace = paste('<birthRate spec="RealParameter" value="',
                                 lambda_init,'"/>', sep=''), x = inf)
    inf  <- gsub(pattern = "<insertDeathRate/>",
                 replace = paste('<deathRate spec="RealParameter" value="',
                                 mu_init,'"/>', sep=''), x = inf)
    inf  <- gsub(pattern = "<insertSamplingRate/>",
                 replace = paste('<samplingRate spec="RealParameter" value="',
                                 psi,'"/>', sep=''), x = inf)
    inf<- gsub(pattern = "<insertBottleneckStrength/>",
               replace = paste('<bottleneckStrength spec="RealParameter" value="',
                               botDuration_init,'"/>', sep=''), x = inf)
    # TnT parameters
    inf  <- gsub(pattern = "<insertPopulationSizes/>",
                 replace = paste('<parameter id="popSize" name="stateNode" value="',
                                 Ne_init,'"/>', sep=''), x = inf)
    inf  <- gsub(pattern = "<insertPairwiseCoalProbBottleneck/>",
                 replace = paste('<parameter id="pairwiseCoalProbBottleneck" name="stateNode" value="',
                                 probCoalBottleneck_init,'"/>', sep=''), x = inf)
  
    # SABD Priors
    inf  <- gsub(pattern = "<insertR0Prior/>",
                 replace = paste0('<distr spec="beast.math.distributions.Uniform"
                            lower="',l_R0,'" upper="',h_R0,'" offset="0."/>'), x = inf)
    inf  <- gsub(pattern = "<insertBecominUninfectiousPrior/>",
                 replace = paste0('<distr spec="beast.math.distributions.Uniform"
                            lower="',l_delta,'" upper="',h_delta,'" offset="0."/>'), x = inf)
    inf  <- gsub(pattern = "<insertRemovalPrior/>",
                 replace = paste0('<distr spec="beast.math.distributions.Uniform"
                            lower="',l_r,'" upper="',h_r,'" offset="0."/>'), x = inf)
    inf  <- gsub(pattern = "<insertOriginPrior/>",
                 replace = paste0('<distr spec="beast.math.distributions.Uniform"
                            lower="0." upper="1000." offset="0."/>'), x = inf)
    inf  <- gsub(pattern = "<insertRhoPrior/>",
                 replace = paste0('<distr spec="beast.math.distributions.Uniform"
                            lower="',l_rho,'" upper="',h_rho,'" offset="0."/>'), x = inf)
  
    # TnT Priors
    inf  <- gsub(pattern = "<insertBottleneckPrior/>",
                 replace = '<distr spec="beast.math.distributions.Uniform"
                            lower="0.0" upper="1.0" offset="0."/>', x = inf)
    inf  <- gsub(pattern = "<insertPopSizesPrior/>",
                 replace = paste0('<distr spec="beast.math.distributions.LogNormalDistributionModel" M="',mean_ne,'" S="',sigma_ne,'" meanInRealSpace="true"/>'), x = inf)
  
  
    writeLines(inf, con=paste0(inf_dir, '/inference.xml'))
  }
  # cmd <- paste0('java -jar /Users/jugne/Documents/Source/tnt_07072021.jar -seed ', seeds[j], ' -overwrite ', inf_dir, "/inference.xml")
  # setwd(inf_dir)
  # result<-system(cmd)
  
  
}

write.csv(true_rates, paste0(wd,"true_rates.csv"))
