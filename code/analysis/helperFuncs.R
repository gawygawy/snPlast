formatSimDat <- function(sim.dat.raw, switch.tr, tr.day){
  sim.dat <- sim.dat.raw %>%
    mutate(rew.found=ifelse((well == 1 | well == 2), 1, 0),
           ID = group_indices(sim.dat.raw, DA, ACh, wmax, iter),
           phase = ifelse(trial <= switch.tr, 0, 1)) %>%
    group_by(ID, phase) %>%
    mutate(ntrial = seq(from=1.1,length.out=length(unique(trial)), by=1/tr.day),
           day=floor(round(ntrial-0.1, 1)), ACh=as.numeric(ACh), DA=as.numeric(DA), phase=as.factor(phase)) %>%
    select(-trial) %>% 
    rename(trial=ntrial) %>% 
    mutate(trial=round(trial, 1)) %>% 
    ungroup() %>%
    group_by(wmax, DA, ACh) %>%
    nest()
}

getTrialSuccess <- function(sim.dat){
  trial.success.percent <- sim.dat %>% 
    unnest() %>% 
    group_by(DA, ACh, phase, trial) %>%
    summarize(predicted=mean(rew.found)) %>%
    mutate(predicted=recode(predicted, `0`=1e-15, `1`=1-1e-15)) %>%
    ungroup() %>%
    group_by(DA, ACh) %>% 
    nest()
}

getModDaySuccess <- function(sim.dat){
  sim.dat %>% 
    unnest() %>% 
    group_by(DA, ACh, phase, day) %>%
    summarize(predicted=mean(rew.found), sem=sqrt((predicted*(1-predicted))/n())) %>%
    mutate(predicted=recode(predicted, `0`=1e-15, `1`=1-1e-15)) %>%
    ungroup() %>%
    group_by(DA, ACh) %>% 
    nest() 
}

getAgentTrialPerf <- function(sim.dat){
agent.trial.success <- sim.dat %>% 
    unnest() %>% 
    select(-well) %>% 
    rename(rew.found.mod = rew.found) %>% 
    group_by(ID) %>% 
    nest()
}

getAgentIds <- function(sim.dat){
  agent.ids <- sim.dat %>%
    unnest() %>% 
    ungroup() %>% 
    distinct(ID, ACh, DA, iter)
  return(agent.ids)
} 

getAgentDaySuccess <- function(agent.trial.success){
  agent.day.success <- agent.trial.success %>%                                
    unnest() %>%
    group_by(ID, phase, day) %>%
    summarize(predicted=mean(rew.found.mod)) %>% 
    ungroup() %>% 
    group_by(ID) %>% 
    nest()
  return(agent.day.success)  
}

calcMse <- function(obs.success, agent.success){
  
  df_new <- left_join(agent.success, obs.success, by = c("phase", "day")) %>% drop_na() 
  err <- df_new %>% 
    mutate(sqErr = (predicted - success.p)^2, wts = 10/n.day) %>% # weight initial learning (8 days) and reversal learning (12 days) equally  
    summarize(w.rmse = sqrt(sum(sqErr*wts)/sum(wts)))

  return(err$w.rmse)
}

addMse <- function(sim.df, obs.df, phase.no=NULL){
  
  if (!is.null(phase.no)){
    sim.df <- sim.df %>% filter(phase==phase.no)
    obs.df <- obs.df %>% mutate(data=map(data, ~ .x[.x$phase==phase.no, ]))
  }

  p1 <- obs.df %>%
    mutate(mse=map(data, calcMse, sim.df))
  
  return(unlist(p1$mse))
}

getAgentMse <- function(trial.success.percent, exp.dat, agent.ids, fit.phase){
  mses <- future_map(trial.success.percent$data, addMse, exp.dat, phase.no=fit.phase)

  mouse.mse <- do.call(rbind, mses)
  colnames(mouse.mse) <- exp.dat$sub

  agent.mses<- bind_cols(agent.day.success %>% select(-data), as_tibble(mouse.mse)) %>% 
    left_join(., agent.ids, by="ID") 

  return(agent.mses)
}

pickTopAgent <- function(.agent.mses, ...){
# either sub, iter or sub, iter, DA 
  agent.params <- .agent.mses %>% 
    gather(key="sub", value="rmse", -c("ACh", "DA", "ID", "iter")) %>%
    group_by(...) %>%
    group_modify(~{
      .x %>%
        top_n(-1, rmse) %>% 
        slice_head() # to deal with ties 
      }) %>% 
    ungroup() %>% 
    mutate(ACh.DA = round(ACh/DA, 2))

    return(agent.params)
} 

pickTopAvgAgent <- function(.agent.mses, ...){
# either sub, iter or sub, iter, DA 
  agent.params <- .agent.mses %>% 
    gather(key="sub", value="rmse", -c("ACh", "DA", "ID", "iter")) %>%
    group_by(sub, DA, ACh) %>%
    summarize(avg.rmse = mean(rmse)) %>%
    group_by(...) %>%
    group_modify(~{
      .x %>%
        top_n(-1, avg.rmse) %>% 
        slice_head() # to deal with ties 
      }) %>% 
    ungroup() %>% 
    mutate(ACh.DA = round(ACh/DA, 2))

    return(agent.params)
} 

find80InMod <- function(data){
 data %>% group_by(phase) %>%
    mutate(hit.80=if_else(any(predicted >= 0.8), which.max(predicted>= 0.8), as.integer(-1))) %>% 
    distinct(hit.80)
}

getEstHit80 <- function(.sub.avg.agent){

  est.sub.hit80 <- .sub.avg.agent %>%  
    group_by(sub) %>% 
    nest() %>% 
    mutate(hit.80.init=map(data, find80InMod)) %>% 
    unnest(hit.80.init) %>% 
    select(-data) %>% 
    ungroup() %>% 
    mutate(hit.80=na_if(hit.80, -1)) %>% 
    mutate(hit.80 = ifelse(phase == 0, replace_na(hit.80, 9), replace_na(hit.80, 13))) %>% 
    left_join(., ids, by="sub")

  return(est.sub.hit80)
}
  
fitGlm <- function(est.trial.dat){
  fit6 <- glm(rew.found.mod ~ condition + phase + trial +
    phase:trial + condition:phase, data=est.trial.dat, family=binomial(link="logit"))

  p.val <- summary(fit6)[[12]][, 4]

  return(p.val)
}

meanfun <- function(data, i){
  d <- data[i, ]
  c(mean(d$ACh), mean(d$DA), mean(d$ACh.DA)) 
}

bo <- function(data){
  b <- boot(data, statistic=meanfun, R=1000)
  ach.ci <- boot.ci(b, conf=0.95, type="bca", index=1)$bca[, 4:5]
  da.ci <- boot.ci(b, conf=0.95, type="bca", index=2)$bca[, 4:5]
  achda.ci <- boot.ci(b, conf=0.95, type="bca", index=3)$bca[, 4:5]

  c(ach.ci, da.ci, achda.ci)
}

fitLearningRateToSim <- function(df){

  models <- plyr::dlply(df, "phase", function(x) glm(rew.found ~ trial, data = x, family=binomial(link="logit")))

  coefs <- plyr::ldply(models, coef)

  return(coefs)
}

averageOverIters <- function(paramSamples, ids){
  avgParams <- paramSamples %>% 
    group_by(sub) %>% 
    rename(iterACh = ACh, iterDA = DA) %>% 
    summarize(ACh = mean(iterACh), DA = mean(iterDA), err=mean(rmse)) %>% 
    mutate(ACh.DA = ACh/DA) %>% 
    left_join(., ids, by="sub")

  return(avgParams)
}
