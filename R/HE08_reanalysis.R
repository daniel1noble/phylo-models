#Chamberlain et al. meta-analysis re-run
#model:HE08
#Honnay et al. 2008 Oikos
#meta-analysis: Genetic variability (Fst) in fragmented versus non-fragmented plant populations
#meta-regression categories: persistant (pe/1) and transient (tr/0)

#####
#setup
#####
#load packages
pacman::p_load(metafor, MCMCglmm, kableExtra, tidyverse, phytools, ape)
#load dataset
Honnay<- read.table("datasets/HonnayEtal2008Oikos_data.txt")
colnames(Honnay) <- c("spp", "es", "sv", "inc")
# Read the trees
Honnay_tree <- read.tree("datasets/HonnayEtal2008Oikos_tree.txt")
# Plot the trees to have a look
plot(Honnay_tree, type = "phylogram", main = "Honnay et al. 2010")
# We need to check trees are ultmetric
is.ultrametric(Honnay_tree) # TRUE
#covert trees to correlation matrices
# For metafor, we use a correlation matrix
A_metafor_Honnay <- ape::vcv(Honnay_tree, corr = TRUE)
# For MCMCglmm we use the inverse of the (co)variance matrix. Need to be scaled. 
A_MCMCglmm_Honnay <- MCMCglmm::inverseA(Honnay_tree, scale =TRUE)$Ainv



#####
#RE models
#####

# Metafor base RE model
summary(Honnay_RE <- rma(es ~ 1, vi = sv, data = Honnay))
#effect 0.16 (0.09, 0.22)


#MCMCglmm base RE model
prior = list(R = list(V = 0.0001, fix = 1), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 100)))
HonnayMCMC_RE <- MCMCglmm(es ~ 1, mev = Honnay$sv, random = ~spp, data = Honnay, verbose = FALSE, prior = prior)
summary(HonnayMCMC_RE)
#ES 0.16 (0.09,0.22)
#our results in metafor and MCMCglmm match Chamberlain

#meta-regression
#categories are : persistant (pe/1) and transient (tr/0)

#metafor
summary(Honnay_REmod <- rma(es ~factor(inc)-1, vi = sv, data = Honnay))

#MCMCglmm
HonnayMCMC_REmod <- MCMCglmm(es ~ factor(inc)-1, mev = Honnay$sv, random = ~spp,
                             data = Honnay, verbose = FALSE, prior = prior)
summary(HonnayMCMC_REmod)
#fac 0(tr) = 0.18 (0.10,0.27), fac 1(pe) = 0.12 (0.02,0.23) 
#our results in metafor and MCMCglmm match Chamberlain



#####
#controlling for phylo
#####

#overall model 
#metafor
summary(Honnay_phy <- rma.mv(es ~ 1, V = sv, random = list(~1|spp), R = list(spp = A_metafor_Honnay),
                             data = Honnay))
#MCMCglmm
prior = list(R = list(V = 0.0001, fix = 1), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 100)))
HonnayMCMC_phy <- MCMCglmm(es ~ 1, mev = Honnay$sv, random = ~spp, data = Honnay, verbose = FALSE, 
                           ginverse = list(spp = A_MCMCglmm_Honnay), prior = prior)
summary(HonnayMCMC_phy)

#again, we're still getting comparable results between metafor, MCMC and Chamberlain
#Chamberlain did not find a difference when controlling for phylo and neither did we


#meta-regression
#categories are : persistant (pe/1) and transient (tr/0)
#metafor
summary(Honnay_phymod <- rma.mv(es ~factor(inc)-1, V = sv,random = list(~1|spp), R = list(spp = A_metafor_Honnay)
                            , data = Honnay))
#MCMCglmm
HonnayMCMC_phymod <- MCMCglmm(es ~ factor(inc)-1, mev = Honnay$sv, random = ~spp, data = Honnay, verbose = FALSE, 
                              ginverse = list(spp = A_MCMCglmm_Honnay), prior = prior)
summary(HonnayMCMC_phymod)
# fac 0 (tr): 0.18 (0.10,0.27), fac1 (pe): 0.12 (0.02,0.23)

#metafor and MCMC models roughly match
#But, this is where we differ from Chamberlain
#Chamberlain found that the ES of fac0(tr) switched direction when phylo controlled for
#however, we find controlling for phylo does not change the effect size here

#####
#model comparison
#####

#overall
AIC(Honnay_phy,Honnay_RE)
#same AIC

#metaregression
AIC(Honnay_phymod,Honnay_REmod)
#same AIC


#conclusion
#we find controlling for phylogeny has no influence on ES for any group, and does not affect model fit either
#Chamberlain found a significantly different effect in the tr group when phylo controlled for


