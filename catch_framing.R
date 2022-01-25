library(tidyverse)
library(baseballr)
library(devtools)
library(retrosheet)
library(readr)
library(Lahman)
library(lme4)
library(doParallel)
library(rstanarm)
library(stringr)
library(knitr)
library(merTools)

###### Re-creating the Baseball Prospectus Method
## for catcher framing; Using STAN ##

### way of learning STAN & advance Bayesian knowledge ###

### https://www.baseballprospectus.com/news/article/38289/bayesian-bagging-generate-uncertainty-intervals-catcher-framing-story/

setwd("~/Desktop/bayesian_baseball")
df<-read.csv("mlb_2020_statcast_pitcher.csv")
player_ref<-read.csv("SFBB-Player-ID-Map.csv")

unique(df$description)

df<-df %>%
  mutate(cs = ifelse(description == "called_strike",1,0),
         pitch_diff_sc = fld_score - bat_score,
         batter = as.factor(batter),
         pitcher = as.factor(pitcher),
         fielder_2 = as.factor(fielder_2),
         umpire = as.factor(umpire))

sum(is.na(df$pitch_diff_sc))

frame_mod_mle<-lme4::glmer(factor(cs) ~ pitch_diff_sc +
                       (1|fielder_2) + ### catcher ID
                       (1|pitcher) + 
                       (1|batter) ,
                     data = df,
                     family = binomial(link = "probit")
                       )

summary(frame_mod_mle)


##### mcmc multi-level ####
frame.mod.mcmc<-stan_glmer(cs ~
                        pitch_diff_sc + 
                          (1|fielder_2) +
                          (1|pitcher) +
                          (1|batter),
                        data = df,
                        family=binomial(link="probit"),
                        chains = 4,
                        prior_intercept = normal(0,10),
                        prior = normal(0,1),
                        prior_aux = exponential(1),
                        prior_covariance = decov(1,1,1,1),
                        adapt_delta = .8,
                        iter = 100,
                        QR = FALSE)

saveRDS(frame.mod.mcmc,"catching_frame_mcmc.rds")
### warning messages ###
# 1. largest r-hat is 1.19 indicating chains have not mixed
# 2. Bulk Effective Sample Sizes (ESS) is too low, indicating posterior means and medians may be unreliable
# 3. Tail Effective Sample Size is too low, indicating posterior variances and tail quantiles may be unreliable


##### Interrogate MCMC model ####
summ.vals<-as.data.frame(summary(frame.mod.mcmc))
summ.vals$parms<-row.names(summ.vals)

catcher.vals<-dplyr::filter(
  summ.vals,stringr::str_detect(
    parms,'fielder_2'
  ))

catcher.vals$catcher<-stringr::str_sub(
  catcher.vals$parms,25,30)

player_ref$MLBID<-as.character(player_ref$MLBID)


catcher.vals<-catcher.vals %>%
  dplyr::left_join(player_ref,by=c("catcher" = "MLBID")) %>%
  dplyr::select(-c(parms)) %>%
  dplyr::select(PLAYERNAME,catcher,mean,mcse,sd,`10%`,`50%`,`90%`,n_eff,Rhat) %>%
  dplyr::arrange(desc(mean))


##### Query MCMC model ######
frame.mcmc.preds<-frame.mod.mcmc$data

SDs<-data.frame(frame.mod.mcmc$ses)
SDs$parms<-row.names(SDs)

catcher.SD<-dplyr::filter(SDs,str_detect(
  parms,'fielder_2'
))

catcher.coefs<-ranef(frame.mod.mcmc)[[3]]
catcher.coefs$catcher<-row.names(catcher.coefs)

names(catcher.coefs)[1]<-"catcher.eff"
catcher.coefs$catcher.sd<-catcher.SD[,1]

frame.mcmc.preds$fitted.values<-frame.mod.mcmc$fitted.values
frame.mcmc.preds$full.lin.preds<-frame.mod.mcmc$linear.predictors

###### Convert group-level intercepts to probability scale ####
frame.mcmc.preds<-frame.mcmc.preds %>%
  dplyr::left_join(catcher.coefs,by=c("fielder_2" = "catcher"))

frame.mcmc.preds$wo_catcher<-with(frame.mcmc.preds,
                                pnorm(full.lin.preds - catcher.eff)  
                                  )

frame.mcmc.preds$wo_catcher_plus1SD<-with(frame.mcmc.preds,
                                          pnorm(full.lin.preds - catcher.eff + catcher.sd))


##### summarize predictor values ####
catcher.mcmc.final<-frame.mcmc.preds %>%
  dplyr::group_by(fielder_2)%>%
  dplyr::summarise(
    mcmc.mean = mean(fitted.values - wo_catcher),
    mcmc.sd = mean(wo_catcher_plus1SD - wo_catcher),
    chances = n()
  ) %>%
  dplyr::left_join(player_ref,by=c("fielder_2" = "MLBID")) %>%
  dplyr::select(PLAYERNAME,chances,mcmc.mean,mcmc.sd) %>%
  arrange(desc(mcmc.mean))

catcher.mcmc.final %>%
  head(10L)


###### simulate Posterior #####
RE.sims<-merTools::REsim(frame_mod_mle,n.sims=1000,seed=1234)

catcher.sims<-dplyr::filter(RE.sims,groupFctr == "fielder_2")

catcher.sims<-catcher.sims %>%
  rename(catcher = groupID) %>%
  dplyr::select(catcher,mean,median,sd)

##### extra coefficients ###
frame.sim.preds<-frame.mcmc.preds
frame.sim.preds<-frame.sim.preds %>%
  dplyr::left_join(catcher.sims,by=c("fielder_2" = "catcher"))


frame.sim.preds$all_frame<-predict(frame_mod_mle,type="response")
frame.sim.preds$all_catcher_mean<-pnorm(predict(frame_mod_mle,
                                                type='link',
                                                re.form = ~
                                                  (1|pitcher) +
                                                  (1|batter)) +
                                          frame.sim.preds$mean)

frame.sim.preds$wo_catcher_plus1SD<-pnorm(predict(
  frame_mod_mle,
  type='link',
  re.form = ~ 
    (1|pitcher) +
    (1|batter)) +
    frame.sim.preds$catcher.sd)

frame.sim.preds$wo_catcher<-predict(frame_mod_mle,
                                    type="response",
                                    re.form = ~ 
                                      (1|pitcher) +
                                      (1|batter))

##### Compile and compare simulation results to MCMC ######
catcher.sim.final<-frame.sim.preds %>%
  dplyr::group_by(fielder_2) %>%
  dplyr::summarise(
    sim.mean = mean(all_frame - wo_catcher),
    sim.sd = mean(wo_catcher_plus1SD - wo_catcher),
    chances = n()
  ) %>%
  dplyr::left_join(player_ref,by=c("fielder_2" = "MLBID")) %>%
  dplyr::select(PLAYERNAME,chances,sim.mean,sim.sd)

catcher.sim.mcmc<-catcher.sim.final %>%
  dplyr::inner_join(catcher.mcmc.final,
                    by=c('PLAYERNAME','chances')) %>%
  dplyr::select(PLAYERNAME,chances,
                mcmc.mean,sim.mean,mcmc.sd,sim.sd) %>%
  dplyr::arrange(desc(sim.mean))

catcher.sim.mcmc


###### Bayesian Block Bootstrap ######



