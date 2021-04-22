# R version 3.6.1  

library(tidyverse)
library(dplyr)
library(furrr)
library(Metrics)
library(lme4)
library(ggeffects)
library(ARTool)
library(emmeans)
library(ggpubr)
library(boot)

source("helperFuncs.R")
source("plotFuncs.R")

plan(multisession)

# load data from simulations  
# iteration number; trial number; well (1 correct, 0 not found, -2 wrong); ACh; DA; wmax 
sim.dat.raw <- readRDS("../../data/simulations/simdat.rds")

fit.phase <- NULL # to fit to both phases. 0 to fit to initial learning only, 1 for reversal learning only 
tr.day <- 10 # number of trials per day 
switch.tr <- 80 # reward was shifted to new location after 80 trials 

# simulated data 
sim.dat <- formatSimDat(sim.dat.raw, switch.tr, tr.day) # group data according to each parameter combi, 561 in total 
trial.success.percent <- getTrialSuccess(sim.dat) # % success per trial 
day.success.percent <- getModDaySuccess(sim.dat) # % success per day 

# load experimental data 
datalist <- readRDS("../../data/experiments/processed_trial_data.rds")

exp.dat <- datalist$combined.dat # outcome of each trial for each mouse  
ids <- datalist$ids # subject ids and groups  
m.success.perc <- datalist$mouse.success.p # % of successful trials per day for each mouse 
n.obs <- nrow(exp.dat %>% unnest())

# catergorize mice according to how quickly an 80% success rate was attained 
m.hit.80 <- datalist$hit.80.ids %>% 
	spread(phase, day.hit) %>% 
	mutate(slow.learn = ifelse(`0` >= 6, 1, 0), slow.rev = ifelse(`1` >= 8, 0.5, 0)) %>% 
	mutate(cat = slow.learn + slow.rev)

m.hit.80$cat <- recode(m.hit.80$cat, `0` = "fast", `1.5` = "slow", `1` = "slow.learn", `0.5` = "slow.rev")

# Figure 3c
plotLearnCats(m.hit.80)
# ==============================================================
# prep model data for fitting 

agent.trial.success <- getAgentTrialPerf(sim.dat) # 

agent.ids <- getAgentIds(sim.dat) 

agent.day.success <- getAgentDaySuccess(agent.trial.success)

## ====== Fig 2C: reducing ACh leads to poorer reversal performance 

group.eg.dat <- agent.day.success %>% 
	left_join(., agent.ids, by="ID") %>% 
	mutate(ACh.DA = round(ACh/DA, 2)) %>% 
	filter(DA == 0.00115, ACh.DA %in% c(0.3, 0.16)) %>% # ratio of ACh:DA 
	mutate(condition = ifelse(ACh.DA == 0.16, "ON", "OFF")) %>% 
	unnest()

plotGpEstDayPerf(group.eg.dat, c("OFF", "ON"))

group.dat <- datalist$group.success.p %>%
	mutate(n.day = ifelse(phase == 0, 8, 12)) %>%  
	group_by(condition) %>% 
	nest() %>% 
	rename(sub = condition)

# ======================== fitting to individual mice  
# there are 561 * 100 agents (param combi x 100 iter)
# sub mses returns a df of the rmse for each agent
# sub.params: for each iteration, pick the agent with the best rmse

# fit the model to the full task - learning and reversal stages 
sub.mses <- getAgentMse(agent.day.success, m.success.perc, agent.ids, fit.phase)
sub.params <- pickTopAgent(sub.mses, sub, iter)

# fit the model to reversal learning stage only - fit.phase = 1
sub.mses.rev <- getAgentMse(agent.day.success, m.success.perc, agent.ids, 1)
sub.params.rev <- pickTopAgent(sub.mses.rev, sub, iter)

# fit the model to initial learning stage only - fit.phase = 0
sub.mses.learn <- getAgentMse(agent.day.success, m.success.perc, agent.ids, 0)
sub.params.learn <- pickTopAgent(sub.mses.learn, sub, iter)
 
# average 100 agents for each subject
sub.avg.agent <- left_join(sub.params, agent.day.success, by="ID") %>%
	left_join(., ids, by = "sub") %>% 
	unnest() %>%   
 	group_by(sub, phase, day, condition) %>% 
 	rename(agent.pred = predicted) %>% 
	summarize(predicted = mean(agent.pred), sem=sd(agent.pred)/sqrt(n()))

# Can the simulations reproduce behavioural outcomes measured in the experiment?
# Supplementary Figure 8a
plotGpEstDayPerf(sub.avg.agent, c("GFP", "OFF", "ON"))

est.hit80 <- getEstHit80(sub.avg.agent)

# Supplementary Figure 8b
plotHit80Sig(est.hit80)

# comparing number of days taken agents to attain an 80% success rate
m1 <- art(hit.80 ~ condition*phase, est.hit80)
anova(m1)
contrast(emmeans(artlm(m1, "condition:phase"), ~ condition:phase), method="pairwise", interaction=TRUE)

# Analysis of Variance of Aligned Rank Transformed Data

# Table Type: Anova Table (Type III tests) 
# Model: No Repeated Measures (lm)
# Response: art(hit.80)

#                   Df Df.res F value     Pr(>F)    
# 1 condition        2     84  1.6676   0.194887    
# 2 phase            1     84 36.7059 3.7433e-08 ***
# 3 condition:phase  2     84  3.1384   0.048475   *
# ---
# Signif. codes:   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# > contrast(emmeans(artlm(h, "condition:phase"), ~ condition:phase), method
# ="pairwise", interaction=TRUE)
#  condition_pairwise phase_pairwise estimate   SE df t.ratio p.value
#  GFP - OFF          0 - 1              4.38 15.9 84 0.276   0.7833 
#  GFP - ON           0 - 1             30.17 15.2 84 1.983   0.0507 
#  OFF - ON           0 - 1             25.79 12.2 84 2.122   0.0368

# applying logistic regression from experimental data analysis to simulated data 
# group data according to iteration number 
agent.trial.success.params <- left_join(sub.params, agent.trial.success %>% unnest(), by=c("ID", "iter", "ACh", "DA")) %>% 
	left_join(., ids, by="sub") %>% 
	mutate(condition = fct_relevel(condition, c("ON", "OFF", "GFP"))) %>% 
  	group_by(iter) %>% 
  	nest()

# fit the fixed-effects glm 
ps <- agent.trial.success.params %>% 
	mutate(p = future_map(data, fitGlm))

p.vals <- as_tibble(do.call(rbind, ps$p)) %>% 
	  mutate(off.sig = conditionOFF < 0.05, off.reversal.sig = `conditionOFF:phase1` < 0.05, 
    GFP.sig = conditionGFP < 0.05, GFP.reversal.sig = `conditionGFP:phase1` < 0.05, 
    off.effect = ifelse(off.sig | off.reversal.sig , 1, 0), 
    GFP.effect = ifelse(GFP.sig | GFP.reversal.sig, 1, 0))

sum(p.vals$off.reversal.sig) # 81/100 iterations produce a significant effect in the reversal stage between light-on and light-off
sum(p.vals$GFP.reversal.sig) # 71/100

# final parameter estimates for each mouse by averaging over iterations 
est.sub.params <- averageOverIters(sub.params, ids)
plotEstParams(est.sub.params, condition, ACh)
kruskal.test(ACh ~ condition, data=est.sub.params)
#p = 0.45 # no significant difference in estimated ACh between control and light-on groups

# calculate confidence intervals of parameter estimates 
sub.param.nest <- sub.params %>% 
	group_by(sub) %>% 
	nest()

est.params.ci <- sub.param.nest %>% 
	mutate(param.ci = map(data, bo))

cis <- do.call("rbind", est.params.ci$param.ci)
colnames(cis) <- c("ach.low", "ach.hi", "da.low", "da.hi", "achda.low", "achda.hi")

sub.cis <- bind_cols(est.sub.params, as.tibble(cis))

# individual parameter estimates if model was fit to reversal stage only 
est.params.rev <- averageOverIters(sub.params.rev, ids)
plotEstParams(est.params.rev, condition, ACh)

# individual parameter estimates if model was fit to learning stage only 
est.params.learn <- averageOverIters(sub.params.learn, ids)
plotEstParams(est.params.learn, condition, DA)

# merged dataframe of parameters fit separately to learning and reversal stages 
est.params.phase <- left_join(est.params.learn, est.params.rev, by=c("sub", "condition")) %>% 
	select(sub, ACh.x, DA.x, ACh.DA.x, condition, ACh.y, DA.y, ACh.DA.y) %>% 
	rename(ACh0 = ACh.x, DA0 = DA.x, ACh.DA0 = ACh.DA.x, ACh1 = ACh.y, DA1 = DA.y, ACh.DA1 = ACh.DA.y)

# plot change in parameters across stages, faceted by learner category 
est.params.phase.cat <- left_join(est.params.phase, m.hit.80 %>% select(sub, condition, cat, slow.learn, slow.rev))
plotParamsByCat(est.params.phase.cat, ACh.DA0, ACh.DA1, cat)

dat1 <- est.params.rev %>%
	select(-err) %>%  
	mutate(phase = 1, DA = DA + 0.0024)

dat2 <- est.params.learn %>%
	select(-err) %>%  
	mutate(phase = 0)

dat3 <- rbind(dat1, dat2) %>% 
	left_join(., m.hit.80) %>% 
	mutate(phase =factor(phase))

# Supplementary Figure 11
plotParamsByPhase(dat3)

m.dat <- m.success.perc %>% 
	unnest() %>% 
	select(-sem) %>% 
	left_join(., sub.avg.agent, by=c("phase", "day", "sub")) %>% 
	select(-over.80, -n.day)

m.dat.con <- m.dat %>%
	mutate(day = ifelse(phase == 0, day, day+8)) %>% 
	left_join(., dat3, by=c("sub", "condition", "phase")) %>%
	filter(condition == "ON") %>% # change condtion "OFF", "GFP", "ON" 
	mutate(`1` = str_pad(as.character(`1`), 2, "left", 0)) %>% 
	unite("orderSub", c(`1`, sub))

# for supplementary figures 5, 6, 7  
main.plot <- plotMainSubFits(m.dat.con)
insets <- plotInsetSubFits(m.dat.con)

pl2 <- main.plot + insets

# Figure 3d 
select.subs <- c("C1", "C7", "53BR", "C6", "53LR","57BL")	

m.dat.select <- m.dat %>% 
	left_join(., sub.cis, by=c("sub", "condition")) %>%
	select(-err) %>% 
	filter(sub %in% select.subs) %>% 
	mutate(sub = factor(sub, levels=c("57BL", "C1", "C7", "C6", "53BR", "53LR")))

main.plot2 <- plotMainSubFits2(m.dat.select)
insets2 <- plotInsetSubFits2(m.dat.select)

pl3 <- main.plot2 + insets2

# learning rates predicted by the model, for the whole parameter space   
sim.learning.rate <- sim.dat %>% 
  mutate(coefs = future_map(data, fitLearningRateToSim)) %>% 
  unnest(coefs) %>% 
  mutate(ACh.DA = round(ACh/DA, 2))

# Supplementary Figure 3a
# number of days taken to reach 80% as predicted by the model, for the whole parameter space 
sims.hit.80 <- trial.success.percent %>% 
	mutate(h80 = future_map(data, find80InMod)) %>% 
	unnest(h80) %>% 
	mutate(ACh.DA = round(ACh/DA, 2),hit.80=na_if(hit.80, -1))	

#sub.params2 <- sub.params %>% left_join(., ids, by="sub") %>% mutate(ACh.DA = round(ACh/DA, 2))
#plotSubOnLearningRate(sim.learning.rate, sub.params2)
#plotSubOnLearningRate(sim.learning.rate, est.sub.params)

p1 <- plotSubOnHit80(sims.hit.80, 0, sub.cis, 50)
p2 <- plotSubOnHit80(sims.hit.80, 1, sub.cis, 70)
g <- arrangeGrob(p1, p2, ncol=2) 

# =================== To test how likely it is to detect a difference in ACh values
# =================== for our current sample sizes 

# Supplmentary Figure 10a showing simulated draws from a constrained parameter space  
p1b <- plotModHit80(sims.hit.80, 0, 50) +
	geom_rect(aes(xmin=0.00055, xmax=0.00295, ymin=0.1, ymax=0.47), color="black", fill=NA) + 
	geom_rect(aes(xmin=0.00065, xmax=0.00305, ymin=-0.01, ymax=0.38), color="black", fill=NA, linetype="longdash")

p2b <- plotModHit80(sims.hit.80, 1, 70) + 
	geom_rect(aes(xmin=0.00055, xmax=0.00295, ymin=0.1, ymax=0.47), fill=NA, color="black") + 
	geom_rect(aes(xmin=0.00065, xmax=0.00305, ymin=-0.01, ymax=0.38), fill=NA, color="black", linetype="longdash")

g2 <- arrangeGrob(p1b, p2b, ncol=2) 

# simulated samples 
sim.dat1 <- agent.day.success %>% 
	left_join(., agent.ids, by=c("ID")) %>% 
	mutate(ACh.DA = round(ACh/DA, 2)) %>% 
	ungroup() %>%  
	filter(ACh.DA < 0.47)

sim.dat.control<- sim.dat1 %>% 
	filter(ACh.DA > 0.1) 

sim.dat.suppressed <- sim.dat1 %>% 
	filter(ACh.DA < 0.38)

test.param.ps <- 0
test.hit80.ps <- 0

for (i in 1:1000){
	sim.off <- sim.dat.control %>% 
		sample_n(., 16) %>% # we draw 16 parameter sets for light-off
		mutate(con = factor("off"))

	sim.gfp <- sim.dat.control %>% 
		sample_n(., 8) %>% # 8 for gfp 
		mutate(con = factor("gfp"))

	sim.on <- sim.dat.suppressed %>% 
		sample_n(., 21) %>% # 21 for light-on 
		mutate(con = factor("on"))
	
	sim.dat.frame <- rbind(sim.off, sim.gfp, sim.on)
	# any difference in ACh across the different conditions in this sample? 
	test.param.ps[i] <- kruskal.test(ACh ~ con, sim.dat.frame)$p.value 

}

sum(test.param.ps < 0.05) # 566 of 1000

# an example of 10 such simulations is shown in Supplementary Figure 10b  
sim.params <- NULL

for (i in 1:10){
	sim.off <- sim.dat.control %>% 
		sample_n(., 16) %>% 
		mutate(group = "off")

	sim.gfp <- sim.dat.control %>% 
		sample_n(., 8) %>% 
		mutate(group = "gfp")

	sim.on <- sim.dat.suppressed %>% 
		sample_n(., 21) %>% 
		mutate(group = "on")
	
	sim.dat.frame <- rbind(sim.off, sim.gfp, sim.on) %>% 
		mutate(iter = paste("sample_", i))
	sim.params <- rbind(sim.params, sim.dat.frame) %>% 
		mutate(group = factor(group, levels = c("gfp", "off", "on")))
}

my_comparisons <- list(c("gfp", "off"), c("gfp", "on"), c("on", "off"))
p5 <- ggboxplot(sim.params, x="group", y="ACh", facet.by="iter", add="dotplot") + 
	stat_compare_means(label = "p.signif", comparisons = my_comparisons, size=3)

ggpar(p5, ylim=c(0, 0.002))

#================================= Parameter identifiability  
flatten.sim.dat <- agent.day.success %>% 
	  unnest() 

id.list <- unique(flatten.sim.dat$ID)
set.seed(189)
sample.ids <- sample(id.list, 200)

sample.sim.frame <-	flatten.sim.dat %>% 
	filter(ID %in% sample.ids) %>%
	rename(sub = ID, success.p = predicted) %>% 
	mutate(n.day = ifelse(phase == 0, 8, 12))  %>% 
	group_by(sub) %>% 
	nest()

sample.mses <- getAgentMse(agent.day.success, sample.sim.frame, agent.ids, fit.phase)
recovered.params <- pickTopAgent(sample.mses, ID, iter) %>% 
	rename(DA.rec = DA, ACh.rec = ACh) %>% 
	filter(rmse != 0)

ground.truth <- agent.ids %>% 
	filter(ID %in% sample.ids) %>% 
	select(-iter)

est.recovered.params <- recovered.params %>% 
	group_by(sub) %>% 
	summarize(ACh.r = mean(ACh.rec), DA.r=mean(DA.rec)) %>% 
	mutate(AChDA.r = ACh.r/DA.r, sub = as.integer(sub))

recovered.from.sim <- left_join(ground.truth, est.recovered.params, by =c("ID" = "sub")) %>% 
	mutate(AChDA = ACh/DA) %>% 
	gather(var, estVal, -ID) %>% 
	separate(var, c("param", "recovered")) %>% 
	mutate(recovered = replace_na(recovered, "o")) %>% 
	spread(recovered, estVal) %>% 
	rename("recovered" = "r", "groundTruth" = "o")

# Supplementary Figure 9
plotRecoveredParams(recovered.from.sim)

# factors influencing estimated ACh  
est.params.cat <- left_join(est.sub.params, m.hit.80 %>% 
			    select(sub, condition, cat, slow.learn, slow.rev)) %>% 
				mutate(condition = fct_relevel(condition, c("ON", "OFF", "GFP")))

m2 <- (lm(ACh ~ slow.rev*DA, est.params.cat))
anova(m2)


