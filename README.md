# TnT-material

This repository contains the code, needed to perform analyses and make figures presented in manuscript "Integrating Transmission Dynamics and Pathogen Evolution Through a
Bayesian Approach" by U. Stolz, T.G. Vaughan and T. Stadler.

To run all the scripts you will need Java (at least 8), R, python 3, and BEAST2.6.x. You may need change working directory paths in R scripts, found at the beggining of each code file.

## Validation and comparison

### Figure 2
1. Run [validation_3009/tntValidation_full_3009.R](validation_3009/tntValidation_full_3009.R)
2. Run all analyses found in validation_3009/run_*/inf/inference.xml with TnT.jar (or installation of TnT on Beast2)
3. Make plots by running [validation_3009/validationPlots.R](validation_3009/validationPlots.R)
4. Figures at figures/qq*.pdf makes up the Figure 2 of the main text.

### Figure 3(a), 4 and Table 2
1. Run [comparisons_180822/makeTnTSims.R](comparisons_180822/makeTnTSims.R)
2. Run [comparisons_180822/makeSCOTTI.R](comparisons_180822/makeSCOTTI.R)
3. Run comparisons_180822/tnt/run_*/inference.xml with TnT.jar (or installation of TnT on Beast2)
4. Run SCOTTI/run_\*/scotti_run_\*.xml with SCOTTI (used version was 2.0.2)
5. Run [comparisons_180822/make_figures.R](comparisons_180822/make_figures.R)
6. This will make [Figure 3(a)](comparisons_180822/figures/transmissionTimes.pdf) from the main text, [Figure 4 left](comparisons_180822/figures/pr_roc/roc_direct.svg), [Figure 4 right](comparisons_180822/figures/pr_roc/pr_direct.svg) and [Table 2](/comparisons_180822/figures/pr_roc/summary.csv) from the main text.

## HIV analysis

### Figure 3(b), 5, Table 3 and Supplementary Figure 1
1. Add appropriate env and pol sequences to [env_pol_saConstrained_s1_sameNe_run0_noSeqs.xml](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/env_pol_saConstrained_s1_sameNe_run0_noSeqs.xml)
2. Make three identical coppies of this file, changing "noSeqs.xml" to "run0.xml", "run1.xml", "run2.xml".
3. Run these files with different starting seed with TnT (use either the standalone .jar file or TnT installation on Beast2.6.x)
4. Repeat steps 1-3 for [env_pol_saConstrained_s1_sameNe_noBot_run0_noSeqs.xml](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/env_pol_saConstrained_s1_sameNe_noBot_run0_noSeqs.xml)
5. Run [combine_logs1](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/combine_tnt.R) and [combine_logs2](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/combine_tnt.R)
6. Run [analyse transmission 1](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/analyse_transmission.R) and [analyse transmission 2](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/analyse_transmission.R)
7. This will make [Figure 3(b)](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/combined/infectionTimes2.pdf), [Figure 5](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/combined/inferred_transmission.svg) and [Table 3](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/combined/param_summary.csv) from the main text. from the main text, as well as [Supplementary Figure 1 (a)](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_botPrior/combined/comb_log_trace.pdf) and [Supplementary Figure 1 (b)](hiv/inf/final/pol_RT/full_sameNe_SW_Ne_aboveOrigin_strictClock_noBot/combined/comb_log_trace.pdf).  





