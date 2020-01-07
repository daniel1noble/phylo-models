#Chamberlain et al. meta-analysis re-run
#model:PeE10
#meta-analysis: Strength of relationship (slope) between P:C ratio in an organism and its resource
#Persson et al. 2010 Oikos
#calculated overall ES and ES for a moderator (tax group) with 4 categories
#categories are : Zooplankton (zo/1), Aquatic inverts (ai/2),Terrestrial inverts (ti/3), Bacteria (ba/4)

#####
#setup
#####
#load packages
pacman::p_load(metafor, MCMCglmm, kableExtra, tidyverse, phytools, ape)
#load dataset
Persson<- read.table("datasets/PerssonEtal2010Oikos_data.txt")
colnames(Persson) <- c("spp", "es", "sv", "inc")
table(Persson$es)
table(Persson$sv)
#note that 19/28 of the obs have an ES of 0
#one outlier: Mexithauma_quadripaludium that has an ES about an order of mag. higher than any other
#looked back in the lit and seems the outlier ES is probably accurate
#unusual circumstance in a snail that has too much P in diet 
#another note: Chamberlain graph for RE and FE models are the same for this one

# Read the tree
Persson_tree <- read.tree("datasets/PerssonETal2010Oikos_tree.txt")
# Plot the tree to have a look
plot(Persson_tree, type = "phylogram", main = "Persson et al. 2010")
# We need to check trees are ultmetric
is.ultrametric(Persson_tree) # TRUE

#covert trees to correlation matrices
# For metafor, we use a correlation matrix
A_metafor_Persson <- ape::vcv(Persson_tree, corr = TRUE)
# For MCMCglmm we use the inverse of the (co)variance matrix. Need to be scaled. 
A_MCMCglmm_Persson <- MCMCglmm::inverseA(Persson_tree, scale =TRUE)$Ainv


#####
# RE models
#####

#overall effect

#metafor
Persson_RE <- rma(es ~ 1, vi = sv, data = Persson)
#model did not converge

#because metafor did not converge, we will try same model with outlier removed
PerssonSub <- subset(Persson, es!=1.26)
PerssonSub_RE <- rma(es ~ 1, vi = sv, data = PerssonSub)
summary(PerssonSub_RE)
#model runs but tau is estimated as 0
#estimate = 0.059 (-0.09,0.21)
#we get an es almost an order of magnitude lower than Chamberlain
#possibly because we excluded the outlier


#meta-regression

#metafor
#original dataset
Persson_REmod <- rma(es , vi = sv,mods= ~ factor(inc)-1, data = Persson)
#model did not converge

#outlier removed
PerssonSub$inc<-as.factor(PerssonSub$inc)
PerssonSub_REmod <- rma(es , vi = sv,mods= ~ factor(inc)-1, data = PerssonSub)
summary(PerssonSub_REmod)
#cannot estimate heterogeneity (tau=0)
#NS results for each category
#this is in contrast to Chamberlain who found a highly significant result for category 2
#the large es for ai (2) in Chamberlain's seems to be driven by the estimate we removed due to non-convergance 


#####
#controlling for phylo
#####

#overall effect

#metafor
Persson_phy <- rma.mv(es ~ 1, V = sv, random = list(~1|spp), R = list(spp = A_metafor_Persson), data = Persson)
summary(Persson_phy)
#effect of 0.13 (-0.6, 0.9)
#higher than  our non-phylo metafor estimate because we include the outlier here
#our estimate is much smaller than Chamberlain (~0.7), and has much larger CI...making it NS

#removing outlier (for comparison with non-phylo model)
Persson_phySub <- rma.mv(es ~ 1, V = sv, random = list(~1|spp), R = list(spp = A_metafor_Persson), data = PerssonSub)
summary(Persson_phySub)
#estimate = 0.06 (-0.09,0.21)
#we get a similar result whether include phylo or not 
#either way we are getting an estimate about 7 times lower and with a much higher SE then they are getting
#our phylo analysis and non-phylo have non-overlapping confidence intervals

#meta-regression
#categories are : Zooplankton (zo/1), Aquatic inverts (ai/2),Terrestrial inverts (ti/3), Bacteria (ba/4)

#metafor - entire dataset
summary(Persson_phymod <- rma.mv(es , V = sv,random = list(~1|spp), R = list(spp = A_metafor_Persson),
                                  mods= ~ factor(inc)-1, data = Persson))
#the main difference from non-phylo is factor 2 (ai) where we see large increase in ES (although effect still NS)
#this is driven by the inclusion of outlier driving factor 2 estimate here (whereas removed from other meta-a)
#our estimate for factor 2 (0.67, NS) is much lower than theirs (1.4) and has larger CI

#metafor - removing outlier
summary(PerssonSub_phymod <- rma.mv(es , V = sv,random = list(~1|spp), R = list(spp = A_metafor_Persson),
                                 mods= ~ factor(inc)-1, data = PerssonSub))


#####################################
#MCMCglmm
########################################

#####
#RE model
#####

#full dataset
prior = list(R = list(V = 0.0001, fix = 1), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 100)))
PerssonMCMC_RE <- MCMCglmm(es ~ 1, mev = Persson$sv, random = ~spp, data = Persson, verbose = FALSE, prior = prior)
summary(PerssonMCMC_RE)
#ES 0.16

#outlier removed (to compare to metafor results)
PerssonMCMC_REsub <- MCMCglmm(es ~ 1, mev = PerssonSub$sv, random = ~spp, data = PerssonSub, verbose = FALSE, prior = prior)
summary(PerssonMCMC_REsub)
#ES is 0.059 (-0.1, 0.2) - identical to metafor results


#meta-regression
#categories are : Zooplankton (zo/1), Aquatic inverts (ai/2),Terrestrial inverts (ti/3), Bacteria (ba/4)
Persson$inc<-as.factor(Persson$inc)
PerssonMCMC_REmod <- MCMCglmm(es ~ factor(inc)-1, mev = Persson$sv, random = ~spp,
                              data = Persson, verbose = FALSE, prior = prior)
summary(PerssonMCMC_REmod)
#here, unlike metafor but like Chamberlain we find a sig effect for group 2 (0.57, CI: 0.13, 1.04)
#in metafor we excluded outlier though, so not comparable
#magnitude of effect is about half what Chamberlain found with much bigger CI still

#meta-reg with outlier removed (to compare to metafor results)
PerssonMCMC_REmodSub <- MCMCglmm(es ~ factor(inc)-1, mev = PerssonSub$sv, random = ~spp,
                              data = PerssonSub, verbose = FALSE, prior = prior)
summary(PerssonMCMC_REmodSub)
#ES in group two changes from 0.6 to 0.0005 with the removal of that single point
#roughly matches metafor result
#single point is driving significance of the group

#####
#phylo model
#####

#full dataset
prior = list(R = list(V = 0.0001, fix = 1), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 100)))
PerssonMCMC_phy <- MCMCglmm(es ~ 1, mev = Persson$sv, random = ~spp, data = Persson, verbose = FALSE, 
                            ginverse = list(spp = A_MCMCglmm_Persson), prior = prior)
summary(PerssonMCMC_phy)
#ES is 0.13 (-0.8, 1.0), same as metafor with slightly higher CI
#larger estimate for phylo than RE (matches Chamberlain)
#effect still 7 times smaller than Chamberlain and NS


#meta-regression
PerssonMCMC_phymod <- MCMCglmm(es ~ factor(inc)-1, mev = Persson$sv, random = ~spp, data = Persson, verbose = FALSE, 
                            ginverse = list(spp = A_MCMCglmm_Persson), prior = prior)
summary(PerssonMCMC_phymod)
#slightly larger ES for group 2 (0.64) when phylo included but becomes NS due to large CI


#maybe the y-axis in Chamberlain fig4(u) is just mislabeled and so all their estimates appear larger than they are?
#the CI we estimate are much larger than those estimated in Chamberlain
#the large CI we estimate results in no sig difference between RE vs. phylo in our analysis

#####
#model comparison
#####

#overall
AIC(PerssonSub_phy,PerssonSub_RE)
#same AIC

#metaregression
AIC(PerssonSub_phymod,PerssonSub_REmod)
#same AIC










