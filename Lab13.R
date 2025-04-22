library(tidyverse)
library(ggplot2)
#Load in Data
dat.finch = read.csv("zebrafinches.csv")

#Question 1
library(moments)#used for calculating statistics

n <- length(dat.finch$further)
x_bar <- mean(dat.finch$further)
s <- sd(dat.finch$further)
t_obs <- x_bar / (s / sqrt(n)) #t-value
skew <- skewness(dat.finch$further)

#Could use this instead shows same thing
#(ttest <- t.test(x = dat.finch$further,
                # mu = 0,
                # alternative = "less")) 
# Gaussian PDF and CDF at t
fz <- dnorm(t_obs)
Fz <- pnorm(t_obs)

# Edgeworth approximation error
edgeworth_error <- (skew / sqrt(n)) * ((2 * t_obs^2 + 1) / 6) * fz


