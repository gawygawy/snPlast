# R version 3.6.1 

rm(list = ls())

library(tidyverse)
library(lme4)
library(lmerTest)
library(furrr)
library(ggeffects)
library(ARTool)
library(emmeans)

source("plotFuncs.R")

# load experimental data 
datalist <- readRDS("../../data/experiments/processed_trial_data.rds")

#Figure 1b
plotGpPerf(datalist$group.success.p)

#Figure 1c
days.to.80 <- datalist$hit.80 %>% 
  mutate(phase = factor(phase), condition = factor(condition))

plotHit80Sig(days.to.80)

# comparing number of days to reach 80% success rate, during initial learning ...
learn.dat <- days.to.80 %>% filter(phase == "Learning") %>% droplevels()

learn.m <- art(hit.80 ~ condition, data = learn.dat )
learn.aov <- anova(learn.m) # n.s.

# and during reversal learning 
reversal.dat <- days.to.80 %>% filter(phase == "Reversal") %>% droplevels()

rev.mod <- art(hit.80 ~ condition, data = reversal.dat )
rev.aov <- anova(rev.mod)
rev.lm <- artlm(rev.mod, "condition")
rev.emm <- emmeans(rev.lm, ~ condition)
contrast(rev.emm , method="pairwise", adjust="none")
# contrast  estimate   SE df t.ratio p.value
# GFP - OFF    0.531 5.21 42  0.102  0.9193 
# GFP - ON   -10.717 5.00 42 -2.143  0.0380 
# OFF - ON   -11.249 3.99 42 -2.816  0.0074 

eff_size(rev.emm, sigma = sigma(rev.lm), edf = df.residual(rev.lm))
# contrast  effect.size    SE df lower.CL upper.CL
# GFP - OFF      0.0441 0.433 42    -0.83   0.9180
# GFP - ON      -0.8903 0.427 42    -1.75  -0.0292
# OFF - ON      -0.9344 0.347 42    -1.63  -0.2338

#Figure 1c
plotGpPerf(datalist$group.success.p)

trial.dat <- datalist$combined.dat %>% 
  filter(!(sub %in% c("GFP", "OFF", "ON"))) %>% 
  unnest()

trial.dat <- trial.dat %>% 
  mutate(condition = fct_relevel(condition, c("ON", "OFF", "GFP")))

# Building the fixed effects model 
fit1 <- glm(rew.found ~ condition, data=trial.dat, family=binomial(link="logit"))
fit2 <- glm(rew.found ~ condition + phase, data=trial.dat, family=binomial(link="logit"))
fit3 <- glm(rew.found ~ condition + phase + trial, data=trial.dat, family=binomial(link="logit"))
fit4 <- glm(rew.found ~ condition + phase + trial +
  phase:trial, data=trial.dat, family=binomial(link="logit"))
fit4_2 <- glm(rew.found ~ condition + phase + trial +
  condition:phase, data=trial.dat, family=binomial(link="logit"))
fit5 <- glm(rew.found ~ condition + phase + trial +
  condition:phase + phase:trial, data=trial.dat, family=binomial(link="logit"))

summary(fit5) # final model - equation 1 in the manuscript 

# Building the mixed effects model. include mouse-specific terms 
fit1b <- glmer(rew.found ~ condition +
    (0 + trial|sub) + (1|sub), data=trial.dat, family=binomial(link="logit"))

fit2b <- glmer(rew.found ~ condition + phase +
    (0 + trial|sub) + (1|sub), data=trial.dat, family=binomial(link="logit"))

fit3b <- glmer(rew.found ~ condition + phase + trial +
    (0 + trial|sub) + (1|sub), data=trial.dat, family=binomial(link="logit"))

fit4b <- glmer(rew.found ~ condition + phase + trial +
  condition:phase + 
    (0 + trial|sub) + (1|sub), data=trial.dat, family=binomial(link="logit"))

fit5b <- glmer(rew.found ~ condition + phase + trial +
  phase:trial + condition:phase +
    (0 + trial|sub) + (1|sub), data=trial.dat, family=binomial(link="logit"))

anova(fit1b, fit2b, fit3b, fit4b, fit5b) # model comparison 

summary(fit5b) # final model - equation 2 in the manuscript 

anova(fit5, fit5b) # random effects model is much better fit to experimental data 

# fixed effects 
fix.df <- ggpredict(
  fit5b, 
  terms=c("trial[all]", "condition", "phase"), 
  type="fe"
) %>%
  dplyr::select(x, group, predicted, facet) %>%
  filter(!((x > 9) & (facet==0))) %>%
  rename(trial=x, condition=group, phase=facet) %>%
  mutate_at(c("phase", "condition"), as.factor) 

# random effects 
re.df <- ggpredict(
  fit5b, 
  terms=c("trial[all]", "sub", "phase"), 
  type="re"
) %>%
  dplyr::select(x, group, predicted, facet) %>%
  filter(!((x > 9) & (facet==0))) %>%
  rename(trial=x, sub=group, phase=facet, indiv.pred=predicted) %>%
  mutate_at(c("phase", "sub"), as.factor) %>% 
  arrange(sub, phase, trial) %>% 
  left_join(., datalist$ids)

# Figure 3b 
plotLogisticRegression(re.df, fix.df)

# random intercepts and slopes estimated for each subject 
sub.re <- as.data.frame(ranef(fit5b))
sub.re <- left_join(datalist$ids, sub.re, by=c("sub"="grp"))

# Supplementary Figure 2 
plotRandomEffects(sub.re)


