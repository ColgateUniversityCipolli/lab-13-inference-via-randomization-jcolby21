library(tidyverse)
library(ggplot2)
#Load in Data
dat.finch = read.csv("zebrafinches.csv")

#Question 1
#Part A
library(moments)#used for calculating statistics (skewness)

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

#Part B
t_vals <- seq(-10, 10, length.out = 1000)
fz_vals <- dnorm(t_vals)
error_vals <- (skew / sqrt(n)) * ((2 * t_vals^2 + 1) / 6) * fz_vals

error_df <- data.frame(t = t_vals, error = error_vals)

ggplot(error_df, aes(x = t, y = error)) +
  geom_line(color = "blue") +
  labs(title = "Edgeworth Approximation Error across t-values",
       x = "t", y = "Error in P(T ≤ t)") +
  theme_minimal()

#Part C
alpha <- 0.05
target_error <- 0.10 * alpha  # 10% of alpha
t_alpha <- qnorm(alpha)  # for left-tailed test
fz_alpha <- dnorm(t_alpha)

# Solve for n
numerator <- skew * (2 * t_alpha^2 + 1) * fz_alpha
n_required <- (numerator / (6 * target_error))^2


