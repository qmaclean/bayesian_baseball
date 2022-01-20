library(tidyverse)
library(data.table)
library(MASS)
#devtools::install_github("dgrtwo/ebbr")
library(ebbr)
library(splines)
library(purrr)
library(broom)
library(broomExtra)
library(bayesrules)
library(rstan)
library(janitor)
library(bayesplot)
library(rstanarm)
library(doParallel)
library(cowplot)
library(brms)

options(digits = 3)

# simulated data,
# generate a sequence of numbers for each combination of a and b
# to plot the probability density function.
# "\u03B1" unicode for the greek letter alpha
sim <- data.table( a = c( 81, 82, 81 + 100 ),
                   b = c( 219, 219, 219 + 200 ) )
sim <- sim[ , .( x = seq( 0, 1, by = 0.002 ) ), by = .( a, b ) ]

sim[ , `:=`( y = dbeta( x, a, b ),
             parameters = paste0( "\u03B1 = ", a, ", \u03B2 = ", b ) ) ]
sim[ , parameters := factor( parameters, levels = unique(parameters) ) ]
sim

## plot of the distribution
# plot of the distribution
PlotBeta <- function(sim)
{
  ggplot( sim, aes( x, y, color = parameters ) ) + geom_line() +
    xlim( 0, .5 ) + ylab("Density of beta") + theme_bw()
}
PlotBeta( sim = sim[ a == 81, ] )
PlotBeta( sim = sim[a %in% c(81,82)])
PlotBeta(sim = sim)


#### loading in historical data
library(Lahman)
data("Master")
data("Batting")
data("Pitching")




pitchers <- Pitching %>%
  group_by(playerID) %>%
  summarize(gamesPitched = sum(G)) %>%
  filter(gamesPitched > 3)

# Add player names
player_names <- Master %>%
  tbl_df() %>%
  dplyr::select(playerID, nameFirst, nameLast, bats) %>%
  unite(name, nameFirst, nameLast, sep = " ")

# include the "bats" (handedness) and "year" column for later
career_full <- Batting %>%
  filter(AB > 0) %>%
  anti_join(pitchers, by = "playerID") %>%
  group_by(playerID) %>%
  summarize(H = sum(H), AB = sum(AB), year = mean(yearID)) %>%
  inner_join(player_names, by = "playerID") %>%
  filter(!is.na(bats))

# We don't need all this data for every step
career <- career_full %>%
  mutate(avg = H / AB) %>%
  filter(AB > 500)

ggplot(career,aes(avg)) +
  geom_histogram(binwidth = .005)

###
m<-MASS::fitdistr(career$avg,densfun = dbeta,
                  start = list(shape1 = 70,shape2 = 200))

alpha0<-m$estimate[1]
beta0<-m$estimate[2]


ggplot(career) +
  geom_histogram(aes(x=avg,y=..density..),binwidth = .005) +
  xlab("Batting Avg") +
  stat_function(fun = function(x) dbeta(x,alpha0,beta0),
                color="red",size=1)


career<-career %>%
  mutate(est = (alpha0 + H) / (alpha0 + beta0 + AB)) %>%
  arrange(desc(est))


ggplot(career,aes(avg,est,color=AB)) +
  geom_point() +
  geom_hline(yintercept = alpha0 / (alpha0 + beta0),
             color = "red",lty = 2) +
  labs(x = "Batting Avg",y = "Empirical Bayes batting avg") +
  geom_abline(intercept = 0,slope=1,color="red") +
  scale_color_gradient(trans = "log",breaks = 10^(1:4))

plot_beta_binomial(alpha = alpha0,beta = beta0,y= .275,n = 50)

## horizontal dash = 0.259

#### Understanding Credible Intervals
binom.test(x = 1,n=3)$conf.int

## Posterior Distribution
career<-career %>%
  mutate(alpha1 = alpha0 + H,
         beta1 = beta0 + AB - H) %>%
  mutate(eb_estimate = (H + alpha0) / (AB + alpha0 + beta0))

career<-career %>%
  add_ebb_estimate(H,AB,prior_subset = AB >= 500)


padres_2007<-c('bardjo01','blumge01','branyru01',
               'gilesma01','gonzaad01','greenkh01',
               'camermi01','bradlmi01','gilesbr02',
               'sledgte01')


career %>%
  filter(playerID %in% padres_2007) %>%
  mutate(name = reorder(name,.fitted)) %>%
  ggplot(aes(.fitted,name)) +
  geom_point() +
  geom_errorbarh(aes(xmin = .low,xmax = .high)) +
  labs(x = "Estimated batting average (w/ 95% confidence interval)",
       y = "Player")
  

### Hierarchical modeling
career <- career_full %>%
  mutate(avg = H / AB) %>%
  filter(AB > 500)

eb_career_ab <- career %>%
  add_ebb_estimate(H, AB, method = "gamlss",
                   mu_predictors = ~ log10(AB))


eb_career_ab %>%
  filter(AB > 10) %>%
  rename(Raw = .raw,Shrunken = .fitted) %>%
  gather(type,estimate,Raw,Shrunken) %>%
  ggplot(aes(AB,estimate)) +
  geom_jitter(alpha=1/10) +
facet_wrap(~type) +
  scale_x_log10()




eb_career_prior<-career_full %>%
  ebb_fit_prior(H,AB,method = "gamlss",
                mu_predictors = ~ 0 + ns(year,df=5) * bats + log(AB))

fd<-crossing(H = 300,
             AB = 1000,
             year = seq(1885,2021),
             bats = c("L","R"))

augment(eb_career_prior, newdata = fd) %>%
  mutate(prior = .alpha0 / (.alpha0 + .beta0),
         prior.low = qbeta(.025, .alpha0, .beta0),
         prior.high = qbeta(.975, .alpha0, .beta0)) %>%
  ggplot(aes(year, prior, color = bats)) +
  geom_line() +
  geom_ribbon(aes(ymin = prior.low, ymax = prior.high), alpha = .1, lty = 2) +
  ylab("Prior distribution (mean + 95% quantiles)")

### Hypothesis Testing 
test_250<-career_full
test_250 <- test_250 %>%
  add_ebb_estimate(H, AB, method = "gamlss",
                   mu_predictors = ~ log10(AB)) %>%
  add_ebb_prop_test(.250,sort=TRUE)

### player-player A/B test

## Mitch Haniger
haniger<-eb_career_ab %>%
  filter(name == "Mitch Haniger")



haniger_params<-c(haniger$.alpha1,haniger$.beta1)
haniger_params

compare_haniger <- eb_career_ab %>%
  add_ebb_prop_test(haniger_params, approx = TRUE, sort = TRUE)
compare_haniger

#### Mixture Models
career_w_pitchers <- Batting %>%
  filter(AB >= 25, lgID == "NL", yearID >= 1980) %>%
  group_by(playerID) %>%
  summarize(H = sum(H), AB = sum(AB), year = mean(yearID)) %>%
  mutate(isPitcher = playerID %in% pitchers$playerID) %>%
  inner_join(player_names, by = "playerID")

ggplot(career_w_pitchers,aes(H/AB)) +
  geom_histogram()


###### simulations #########
## calculate empirical bayes prior
prior<-career_full %>%
  ebb_fit_prior(H,AB)

prior

alpha0<-tidy(prior)$alpha
beta0<-tidy(prior)$beta

qplot(rbeta(10000,alpha0,beta0)) +
  theme_bw()

ggplot(career,aes(AB)) +
  geom_histogram() +
  scale_x_log10() +
  theme_bw()

set.seed(2017)
career_sim<-career  %>%
  mutate(p = rbeta(n(),alpha0,beta0),
         H = rbinom(n(),AB,p))

career_sim

career_sim_eb<-career_sim %>%
  add_ebb_estimate(H, AB)

career_sim_eb

career_sim_gathered<-career_sim_eb %>%
  rename(Shrunken = .fitted,Raw = .raw) %>%
  gather(type,estimate,Shrunken,Raw)

career_sim_gathered %>%
  filter(AB >= 10) %>%
  ggplot(aes(p, estimate, color = AB)) +
  geom_point() +
  geom_abline(color = "red") +
  geom_smooth(method = "lm", color = "white", lty = 2, se = FALSE) +
  scale_color_continuous(trans = "log", breaks = c(10, 100, 1000, 10000)) +
  facet_wrap(~ type) +
  labs(x = "True batting average (p)",
       y = "Raw or shrunken batting average",
       title = "Empirical Bayes shrinkage reduces variance, but causes bias",
       subtitle = "Red line is x = y; dashed white line is a linear fit")

### mean-squared error and bias relative to AB
career_sim_gathered %>%
  group_by(type) %>%
  summarize(mse = mean((estimate - p)^ 2))

## mse of shrunken is lower than raw, which means the
## method has succeeded

metric_by_bin<- career_sim_gathered %>%
  group_by(type, AB = 10 ^ (round(log10(AB)))) %>%
  summarize(mse = mean((estimate - p) ^ 2))

ggplot(metric_by_bin,aes(AB,mse, color = type)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Number of at-bats (AB)",
       y = "Mean-squared-error within this bin (note log scale)",
       title = "Mean squared error is higher with raw estimate, especially for low AB") +
  theme_bw()


set.seed(2017)


#credible intervals
career_sim_eb %>%
  sample_n(20) %>%
  mutate(name = reorder(name, .fitted)) %>%
  ggplot(aes(.fitted, name)) +
  geom_point() +
  geom_point(aes(x = p), color = "red") +
  geom_errorbarh(aes(xmin = .low, xmax = .high)) +
  labs(x = "Estimated batting average (w/ 95% credible interval)",
       y = "Player",
       title = "Credible intevals for 20 randomly selected players",
       subtitle = "The true batting average of each player is shown in red")

career_sim_eb %>%
  summarize(coverage = mean(.low <= p & p <= .high))

# fit priors
sim_prior<-ebb_fit_prior(career_sim,H,AB)


## FDR control 
pt<-career_sim_eb %>%
  add_ebb_prop_test(.3,sort=TRUE)

hall_of_fame_ba <- pt %>%
  filter(.qvalue <= .1) %>%
  arrange(desc(year))
hall_of_fame_ba

# expect less than 10% of true ba avgs
mean(hall_of_fame_ba$p < .3)

### false discovery rate - fdr
### posterior error probability - PeP
pt %>%
  mutate(true_fdr = cummean(p < .3)) %>%
  ggplot(aes(.qvalue,true_fdr)) +
  geom_line() +
  geom_abline(color="orange") +
  labs(x = "q-value",
       y = "True FDR at this q-value threshold")


### beta-binomial regression 
bb_reg<-career %>%
  ebb_fit_prior(H, AB, method = "gamlss",mu_predictors = ~ log10(AB))

tidy(bb_reg)

set.seed(2017)


career_sim_ab <- augment(bb_reg, career) %>%
 # select(playerID, AB, .alpha0, .beta0) %>%
  mutate(p = rbeta(n(), .alpha0, .beta0),
         H = rbinom(n(), AB, p))

career_ab_prior <- career_sim_ab %>%
  ebb_fit_prior(H, AB, method = "gamlss", mu_predictors = ~ log10(AB))

tidy(career_ab_prior) 

### looks pretty close! ###

##### Replications ####
set.seed(2017)


### this adds 50 simulations to run our method
sim_replications <- career %>%
  crossing(replication = 1:50) %>%
  mutate(p = rbeta(n(), alpha0, beta0),
         H = rbinom(n(), AB, p))

sim_replications

sim_replication_models<-sim_replications %>%
  nest(-replication) %>%
  mutate(prior = map(data, ~ ebb_fit_prior(.,H,AB)))

sim_replication_priors <- sim_replication_models %>%
  unnest(map(prior, tidy), .drop = TRUE)

sim_replication_priors


#### MCMC work
min_ba<-min(career$avg)
max_ba<-max(career$avg)
len<-(max_ba - min_ba) * 1000






