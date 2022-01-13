library(tidyverse)
library(data.table)
library(MASS)

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

career<-Batting %>%
  filter(AB > 0) %>%
  group_by(playerID) %>%
  summarize(H = sum(H),
            AB = sum(AB)) %>%
  mutate(avg = H / AB)

career<-career %>%
  left_join(Master,by=c("playerID")) %>%
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

## horizontal dash = 0.259

#### Understanding Credible Intervals
binom.test(x = 1,n=3)$conf.int

## Posterior Distribution
career<-career %>%
  mutate(alpha1 = alpha0 + H,
         beta1 = beta0 + AB - H) %>%
  mutate(eb_estimate = (H + alpha0) / (AB + alpha0 + beta0))

padres_2007<-c('bardjo01','blumge01','branyru01',
               'gilesma01','gonzaad01','greenkh01',
               'camermi01','bradlmi01','gilesbr02',
               'sledgte01')

career_pads_07<-career %>%
  filter(playerID %in% c('bardjo01','blumge01','branyru01',
                         'gilesma01','gonzaad01','greenkh01',
                         'camermi01','bradlmi01','gilesbr02',
                         'sledgte01'))

### create beta distribution 
career_pads_07<-career_pads_07 %>%
  mutate(PEP = pbeta(.3,alpha1,beta1))



