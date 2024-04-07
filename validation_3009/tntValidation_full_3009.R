# requires bdmm-prime package to be installed on beast: https://github.com/tgvaughan/BDMM-Prime
library(beastio)
library(ape)
library(phytools)
library(ggplot2)

f_lambda <- function(R0, delta){
  return(R0*delta)
}

f_mu <- function(delta, s, r){
  return(delta*(1-s)/(1-(1-r)*s))
}

f_psi <- function(s, delta, r){
  return(s*delta/(1-(1-r)*s))
}

wd<- "~/Documents/Source/TnT-material/validation_3009/"
setwd(wd)

# Source here https://github.com/jugne/TnT
tnt_path <- "../TnT.jar"
dirrectAncestorAnnotator <- "../DirectAncestorAnnotator.jar"


n_runs = 200
n_gene_samples = 5
seed = 64418
set.seed(seed)
seeds <- sample(0:2147483647, n_runs)

# define params of the lognormal distribution of the Ne
mean_ne <- 1;
sigma_ne <- 0.5;
mu_ne <- log(mean_ne) - sigma_ne^2/2

origin_true = 3.0;
min_samples=5;
max_samples=30;


true_rates<-data.frame(seed=numeric(n_runs), t_or=numeric(n_runs), R0=numeric(n_runs), delta=numeric(n_runs), s=numeric(n_runs), r=numeric(n_runs), rho=numeric(n_runs), Ne=numeric(n_runs),
                       Ne_orig=numeric(n_runs), tau=numeric(n_runs), botDuration=numeric(n_runs), pairwiseCoalProbBottleneck=numeric(n_runs))

for (j in 1:n_runs){
  run_dir <- paste0(wd, "full/run_",j)
  dir.create(run_dir)
  sim_dir <- paste0(run_dir,"/sim")
  dir.create(sim_dir)
  inf_dir <- paste0(run_dir,"/inf")
  dir.create(inf_dir)
  
  
  while(T){
  
  # MSC params
  Ne <- rlnorm(1,mu_ne, sigma_ne)
  Ne_orig <- rlnorm(1,mu_ne, sigma_ne)
  probCoalBottleneck <- runif(1, 0.0, 1.0)
  tau <- -log(1-probCoalBottleneck)
  botDuration <- tau*Ne
  
  
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
  rho <- runif(1, l_rho, h_rho)
  
  print(paste(lambda, mu, psi, r, rho, origin_true, seeds[j]))
  
  #########################################################
  ########### Make transmission simulation XML ############
  #########################################################

  true_rates$seed[j] =seeds[j]
  true_rates$t_or[j] = origin_true
  true_rates$R0[j] = R0
  true_rates$delta[j] = delta
  true_rates$s[j] = s
  true_rates$r[j] = r
  true_rates$rho[j] = rho
  true_rates$Ne[j] = Ne
  true_rates$Ne_orig[j] = Ne_orig 
  true_rates$tau[j] = tau
  true_rates$botDuration[j] = botDuration
  true_rates$pairwiseCoalProbBottleneck[j] = probCoalBottleneck
  
  
  tr_sim <- readLines( paste0(wd,"/templates/trTreeTemplate.xml"))
  tr_sim  <- gsub(pattern = "insertMinSamples",
                  replace = paste('minSamples="',
                                  min_samples,'"',sep=''), x = tr_sim)
  tr_sim  <- gsub(pattern = "insertMaxSamples>",
                  replace = paste('maxSamples="',
                                  max_samples,'">', sep=''), x = tr_sim)
  tr_sim  <- gsub(pattern = "<insertOrigin/>",
               replace = paste('<origin id="origin" spec="RealParameter" value="',
                               origin_true,'"/>', sep=''), x = tr_sim)
  tr_sim  <- gsub(pattern = "<insertR0/>",
                  replace = paste('<skylineValues spec="RealParameter" value="',
                                  R0,'"/>', sep=''), x = tr_sim)
  tr_sim  <- gsub(pattern = "<insertBecominUninfectiousRate/>",
                  replace = paste('<skylineValues spec="RealParameter" value="',
                                  delta,'"/>', sep=''), x = tr_sim)
  tr_sim  <- gsub(pattern = "<insertsamplingProportion/>",
                  replace = paste('<skylineValues spec="RealParameter" value="',
                                  s,'"/>', sep=''), x = tr_sim)
  tr_sim  <- gsub(pattern = "<insertRemovalProb/>",
                  replace = paste('<skylineValues spec="RealParameter" value="',
                                  r,'"/>', sep=''), x = tr_sim)
  tr_sim  <- gsub(pattern = "<insertRhoValues/>",
                  replace = paste('<values spec="RealParameter" value="',
                                  rho,'"/>', sep=''), x = tr_sim)
  
  writeLines(tr_sim, con=paste0(sim_dir, '/trTreeSim.xml'))
  
  cmd <- paste0('java -jar ',tnt_path,' -seed ', seeds[j], ' -overwrite ', sim_dir, "/trTreeSim.xml")
  setwd(sim_dir)
  result<-system(cmd, timeout=1)
  if (result!=124){
    break
  }
  }
  
  tree_f <- list.files(path=paste0(run_dir,"/sim"), pattern="trTreeSim.trees", full.names = TRUE)
  tp <- readLines(tree_f)
  newick <- strsplit(tp[length(tp)-1], '= ')[[1]][2]
  tree<- read.tree(text=newick)
  
  matches = regmatches(newick,gregexpr('[0-999]:', newick))
  for(i in 1:length(matches)) {
    
    for(found in matches[[i]]){
      
      newick = sub(pattern = found, 
                   replacement = paste(strsplit(found, ":")[[1]][1],  "_:", sep=""),
                   x=newick)
      
    }
  }
  
  taxon <- ""
  samples <- NULL
  first <- T
  date<-""
  gene_taxon<-""
  taxonInf <- ""
  
  node_order<-data.frame(idx=numeric(length(tree$tip.label)),tipLabel=numeric(length(tree$tip.label)), nodeHeight=numeric(length(tree$tip.label)))
  for(l in 1:length(tree$tip.label)){
    node_order$idx[l]<-l
    node_order$tipLabel[l]<-as.numeric(tree$tip.label[l])
    node_order$nodeHeight[l]<-nodeheight(tree, node=1)-nodeheight(tree, node=l)
  }
  
  node_order<- node_order[order(node_order$nodeHeight, node_order$tipLabel),]
  
  cmd <- paste0('java -jar ',directAncestorAnnotator, " ", sim_dir, "/trTreeSim.trees")
  setwd(sim_dir)
  result<-system(cmd)
  taxonSa <- ""
  taxonSA_idref<-""
  patientTaxonSets<-""
  
  saRestriction <- list.files(path=paste0(run_dir,"/sim"), pattern="directAncestors.log", full.names = TRUE)
  saRestriction <- readLines(saRestriction)
  if (length(saRestriction)>1){
  for(l in 2:length(saRestriction)){
    lineSplit<-strsplit(saRestriction[l], ":")[[1]]
    taxonSa <- paste(taxonSa, paste0('<taxonSet spec="TaxonSet" id="taxon_', as.numeric(strsplit(lineSplit[1], "_")[[1]][1])+1,'_">'), collapse ="\n")
    patientTaxonSets<- paste(patientTaxonSets, paste0('<patientTaxonSets idref="taxon_',  as.numeric(strsplit(lineSplit[1], "_")[[1]][1])+1,'_"/>'), collapse ="\n")
    taxonSA_idref <- paste(taxonSA_idref, paste0('<taxonsets idref="taxon_',  as.numeric(strsplit(lineSplit[1], "_")[[1]][1])+1,'_"/>'), collapse="\n")
    taxonSa <- paste(taxonSa, paste0('<taxon idref="', as.numeric(strsplit(lineSplit[1], "_")[[1]][1])+1, '_"/>'), collapse ="\n")
    lineSplit_ <- strsplit(lineSplit[2], ",")[[1]]
    for (n in 1:length(lineSplit_)){
      taxonSa <- paste(taxonSa, paste0('<taxon idref="', as.numeric(strsplit(lineSplit_[n], "_")[[1]][1])+1, '_"/>'), collapse ="\n")
    }
    taxonSa <- paste(taxonSa, '</taxonSet>', collapse ="\n")
  }
  }
  
  
  for (l in node_order$idx){
    taxon <- paste(taxon, paste0('<taxon spec="Taxon" id="',tree$tip.label[l],'_"/>'), collapse = "\n")
    taxonInf <- paste(taxonInf, paste0('<taxon idref="',tree$tip.label[l],'_"/>'), collapse = "\n")
    if (first){
      samples <- paste(samples, paste0(tree$tip.label[l],'_=',n_gene_samples), sep = "")
      date <- paste(date, paste0(tree$tip.label[l],'_=',nodeheight(tree, node=1)-nodeheight(tree, node=l)) , sep="")
      first <- F
    } else {
      samples <- paste(samples, paste0(tree$tip.label[l],'_=',n_gene_samples), sep = ",")
      date <- paste(date, paste0(tree$tip.label[l],'_=',nodeheight(tree, node=1)-nodeheight(tree, node=l)) , sep=",")
    }
    gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',tree$tip.label[l],'_" spec="TaxonSet">'), collapse = "\n")
    for (kk in 1:n_gene_samples){
      gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',tree$tip.label[l],"_",kk,'" spec="Taxon"/>'), collapse = "\n")
    }
    gene_taxon<- paste(gene_taxon, '</taxon>', collapse = "\n")
    
  }
  
  
  log <- list.files(path=paste0(run_dir,"/sim"), pattern="trTreeSim.log", full.names = TRUE)
  t <- read.table(log[[1]], header=TRUE, sep="\t")

  #########################################################
  ############## Make  gene simulation XML ################
  #########################################################
  sim <- readLines( paste0(wd,"/templates/geneSimTemplate.xml"))
  sim  <- gsub(pattern = "<insertTaxon/>",
               replace = taxon, x = sim)
  sim  <- gsub(pattern = "<insertTransmissionTree/>",
               replace = paste('<transmissionTreeInput spec="TreeParser"
                                newick="',newick,
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
                               botDuration,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertFinalSampleOffset/>",
               replace = paste('<finalSampleOffset spec="RealParameter" value="',
                               t$finalSampleOffset,'"/>', sep=''), x = sim)



  # fBD parameters
  sim  <- gsub(pattern = "<insertBirthRate/>",
               replace = paste('<birthRate spec="RealParameter" value="',
                               lambda,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertDeathRate/>",
               replace = paste('<deathRate spec="RealParameter" value="',
                               mu,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertSamplingRate/>",
               replace = paste('<samplingRate spec="RealParameter" value="',
                               psi,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertSamplingExtantRate/>",
               replace = paste('<samplingExtantRate spec="RealParameter" value="',
                               rho,'"/>', sep=''), x = sim)
  sim  <- gsub(pattern = "<insertOrigin/>",
               replace = paste('<origin spec="RealParameter" value="',
                               origin_true,'"/>', sep=''), x = sim)

  writeLines(sim, con=paste0(sim_dir, '/geneTreeSim.xml'))

  cmd <- paste0('java -jar ',tnt_path,' -seed ', seeds[j], ' -overwrite ', sim_dir, "/geneTreeSim.xml")
  setwd(sim_dir)
  result<-system(cmd)

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
  rho_init <- runif(1, l_rho, h_rho)
  

  
  taxon <- ""
  samples <- NULL
  first <- T
  date<-""
  gene_taxon<-""
  taxonInf <- ""
  
  a <- read.nexus.data(paste0(sim_dir,"/geneTreeSim.alignment.nexus"))
  print(as.numeric(unique(sapply(strsplit(names(a), "_"), "[", 1 ))))
  
  for (l in as.numeric(unique(sapply(strsplit(names(a), "_"), "[", 1 )))){
    taxon <- paste(taxon, paste0('<taxon spec="Taxon" id="',l,'_"/>'), collapse = "\n")
    taxonInf <- paste(taxonInf, paste0('<taxon idref="',l,'_"/>'), collapse = "\n")
    if (first){
      samples <- paste(samples, paste0(l,'_=',n_gene_samples), sep = "")
      date <- paste(date, paste0(l,'_=',nodeheight(tree, node=1)-nodeheight(tree, node=which(tree$tip.label==l))) , sep="")
      first <- F
    } else {
      samples <- paste(samples, paste0(l,'_=',n_gene_samples), sep = ",")
      date <- paste(date, paste0(l,'_=',nodeheight(tree, node=1)-nodeheight(tree, node=which(tree$tip.label==l))) , sep=",")
    }
    gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',l,'_" spec="TaxonSet">'), collapse = "\n")
    for (kk in 1:n_gene_samples){
      gene_taxon <- paste(gene_taxon,  paste0('<taxon id="',l,"_",kk,'" spec="Taxon"/>'), collapse = "\n")
    }
    gene_taxon<- paste(gene_taxon, '</taxon>', collapse = "\n")
    
  }

  inf <- readLines(paste0(wd,"/templates/inference_template.xml"))
  inf  <- gsub(pattern = "<trait>",
               replace = paste0("<trait id='tipDates' spec='beast.evolution.tree.TraitSet' traitname='date-backward' units='year' value='",date,"'>"), x = inf)
  inf  <- gsub(pattern = "<insertTaxonSuperSet/>",
               replace = gene_taxon, x = inf)
  inf  <- gsub(pattern = "<insertTaxon/>",
               replace = taxonInf, x = inf)

  inf  <- gsub(pattern = "<insertOrigin/>",
               replace = paste0('<origin id="origin" spec="RealParameter" value="',length(tree$tip.label)*2,'." lower="0."/>'), x = inf)
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
                               t$finalSampleOffset,'"/>', sep=''), x = inf)

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
  
  #########################################################
  ############ Make restricted inference XML ##############
  #########################################################
  
  inf <- readLines(paste0(wd,"/templates/inference_template_saRestricted.xml"))
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
               replace = paste0('<origin id="origin" spec="RealParameter" value="',length(tree$tip.label)*2,'." lower="0."/>'), x = inf)
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
                               t$finalSampleOffset,'"/>', sep=''), x = inf)
  
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
  
  
  writeLines(inf, con=paste0(inf_dir, '/inference_saRestricted.xml'))
}

write.csv(true_rates, paste0(wd,"full/true_rates.csv"))


nSamples=numeric(n_runs)
nSA = numeric(n_runs)
nRejected = numeric(n_runs)
for (j in 1:n_runs){
  run_dir <- paste0(wd, "full/run_",j)
  log <- list.files(path=paste0(run_dir,"/sim"), pattern="trTreeSim.log", full.names = TRUE)
  t <- read.table(log[[1]], header=TRUE, sep="\t")
  nSamples[j] = t$NonSASampleCount + t$SACount
  nSA[j] = t$SACount
  nRejected[j] = t$numRejectedTrees
}

summarySamples <- data.frame(nSamples=nSamples, nSA=nSA, nRejectedTrTreeSims=nRejected) 
write.csv(summarySamples, paste0(wd,"full/summarySamples.csv"))
