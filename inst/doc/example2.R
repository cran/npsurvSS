## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(npsurvSS)
library(tidyverse)
library(ggplot2)

## ------------------------------------------------------------------------
# Active
arm1 <- create_arm(size=240,
                   accr_time=14,
                   accr_dist="truncexp",
                   accr_param=0.1,
                   surv_scale=per2haz(9),
                   loss_scale=log(1/0.99)^(1/2)/25,
                   loss_shape=2,
                   total_time=25)

# Control
arm0 <- create_arm(size=120,
                   accr_time=14,
                   accr_dist="truncexp",
                   accr_param=0.1,
                   surv_scale=per2haz(6),
                   loss_scale=log(1/0.99)^(1/2)/25,
                   loss_shape=2,
                   total_time=25)

## ---- fig.show='hold'----------------------------------------------------
# Accrual
tibble(
  x = seq(0, 14, 0.1),
  y = paccr(x, arm0)
) %>%
  ggplot(aes(x, y)) +
  geom_line() +
  labs(x = "Time from first patient in (months)",
       y = "Accrual CDF")

# Loss to follow-up
tibble(
  x = seq(0, 25, 0.1),
  y = ploss(x, arm0)
) %>%
  ggplot(aes(x, y)) +
  geom_line() +
  labs(x = "Time from study entry (months)",
       y = "Loss to follow-up CDF")

## ---- fig.width=5--------------------------------------------------------
# Calculate survival curves
x.vec     <- seq(0, 25, 0.1)  # vector of unique x-coordinates
tau1.vec  <- seq(0, 4.5, 1.5) # vector of unique changepoints
arm1t     <- arm1             # initialize active arm
y         <- c()              # vector or all y-coordinates
for (tau1 in tau1.vec) {
  # Update scale and interval parameters for the active arm
  arm1t$surv_scale    <- c(arm0$surv_scale[1], per2haz(9-tau1, 1-0.5/psurv(tau1, arm0, lower.tail=F)))
  arm1t$surv_interval <- c(0, tau1, Inf)
  
  # Calculate y-coordinates
  y <- c(y, psurv(x.vec, arm1t, lower.tail=F))
}

# Visualize survival curves
tibble(
  tau1 = rep(tau1.vec, each=length(x.vec)), 
  x = rep(x.vec, length(tau1.vec)),
  y = y
) %>%
  ggplot(aes(x, y)) +
  geom_line(aes(color=factor(tau1), lty=factor(tau1))) +
  labs(x = "Time from study entry (months)",
       y = "Survival function",
       color = "Changepoint",
       lty = "Changepoint")

## ------------------------------------------------------------------------
tau1.vec <- seq(0, 5.5, 0.5) # vector of changepoints
table_4a <- data.frame(matrix(0, nrow=length(tau1.vec), ncol=7)) # initialize results table
for (r in 1:length(tau1.vec)) {
  
  tau1 <- tau1.vec[r]
  
  # Update scale and interval parameters for the active arm
  arm1t$surv_scale    <- c(arm0$surv_scale[1], per2haz(9-tau1, 1-0.5/psurv(tau1, arm0, lower.tail=F)))
  arm1t$surv_interval <- c(0, tau1, Inf)
  
  # Calculate power and store results
  table_4a[r,] <- c(tau1,
                    power_two_arm(arm0, arm1t,
                                  test = list(list(test="weighted logrank"),
                                              list(test="weighted logrank", weight="n"),
                                              list(test="weighted logrank", weight="FH_p1_q1"),
                                              list(test="survival difference", milestone=11),
                                              list(test="rmst difference", milestone=11),
                                              list(test="percentile difference", percentile=0.5)))$power)
  
}

# Convert table to long format and re-label tests
table_4a <- gather(table_4a, "test", "power", 2:7) %>%
  mutate(test = recode(test, X2 = "wlogrank (1)",
                       X3 = "wlogrank (GB)",
                       X4 = "wlogrank (FH)",
                       X5 = "11-month surv diff",
                       X6 = "11-month rmst diff",
                       X7 = "median diff")) %>%
  as_tibble()
names(table_4a)[1] <- "tau1"

## ------------------------------------------------------------------------
table_4a

## ---- fig.width=5--------------------------------------------------------
ggplot(table_4a, aes(x=tau1, y=power)) +
  geom_line(aes(color=test, lty=test)) +
  labs(x = "tau1",
       y = "Power",
       color = "Test",
       lty = "Test")

