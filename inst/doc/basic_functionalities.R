## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(npsurvSS)
active <- create_arm(size=120, 
                     accr_time=6, 
                     surv_scale=0.0462, 
                     loss_scale=0.00578, 
                     follow_time=12)
control <- create_arm(size=120, 
                     accr_time=6, 
                     surv_scale=0.0578, 
                     loss_scale=0.00578, 
                     follow_time=12)

## ------------------------------------------------------------------------
active <- create_arm(size=120, 
                     accr_time=6, 
                     surv_scale=per2haz(15), # corresponds to 15 month median
                     loss_scale=per2haz(120), # corresponds to 120 month median
                     follow_time=12)
control <- create_arm(size=120, 
                     accr_time=6, 
                     surv_scale=per2haz(12), # corresponds to 12 month median
                     loss_scale=per2haz(120), 
                     follow_time=12)
per2haz(15) # convert median survival to hazard rate
per2haz(0.0462) # convert hazard rate to median survival

## ------------------------------------------------------------------------
active <- create_arm_lachin(size=120, 
                            accr_time=6, 
                            surv_median=15,
                            loss_milestone=c(120, 0.5), # corresponds to 120 month median
                            follow_time=12)
control <- create_arm_lachin(size=120, 
                             accr_time=6, 
                             surv_milestone=c(12, 0.5), # corresponds to 12 month median
                             loss_median=120, 
                             follow_time=12)
class(active)

## ------------------------------------------------------------------------
x <- seq(0, 6, 0.1)
plot(x, paccr(q=x, arm=control), 
     xlab="Time from first patient in (month)", 
     ylab="Accrual CDF",
     type="l")

## ------------------------------------------------------------------------
x <- seq(0, 18, 0.1)
plot(x, psurv(q=x, arm=control, lower.tail=F), 
     xlab="Time from study entry (month)", 
     ylab="Survival function",
     type="l")

## ------------------------------------------------------------------------
power_two_arm(control, active)

## ------------------------------------------------------------------------
# unweighted log-rank
power_two_arm(control, active, test=list(test="weighted logrank"))
# Gehan-Breslow weighted log-rank
power_two_arm(control, active, test=list(test="weighted logrank", weight="n"))
# difference in 12 month survival
power_two_arm(control, active, test=list(test="survival difference", milestone=12))
# ratio of 12 month RMST
power_two_arm(control, active, test=list(test="rmst ratio", milestone=12))

## ------------------------------------------------------------------------
power_two_arm(control, active, test=list(list(test="weighted logrank"),
                                         list(test="weighted logrank", weight="n"),
                                         list(test="survival difference", milestone=12),
                                         list(test="rmst ratio", milestone=12)
))

## ------------------------------------------------------------------------
size_two_arm(control, active, 
             test=list(list(test="weighted logrank"),
                       list(test="weighted logrank", weight="n"),
                       list(test="survival difference", milestone=12),
                       list(test="rmst ratio", milestone=12)
             ))

## ------------------------------------------------------------------------
control_new <- control
active_new  <- active
control_new$size  <- 1
active_new$size   <- 1
size_two_arm(control_new, active_new)

## ------------------------------------------------------------------------
active_new$size   <- 2
size_two_arm(control_new, active_new)

## ------------------------------------------------------------------------
exp_events(control, active) # expected number of events
tau <- exp_duration(control, active, d=150) # study duration for expected number of events to equal d

tau

## ------------------------------------------------------------------------
control_new <- control
control_new$total_time  <- tau
control_new$follow_time <- tau - control_new$accr_time
active_new  <- active
active_new$total_time   <- tau
active_new$follow_time  <- tau - active_new$accr_time

exp_events(control_new, active_new) # check expected number of events

## ------------------------------------------------------------------------
trial1 <- simulate_trial(control, active, duration=18)
head(trial1, 5)
table(trial1$arm, trial1$reason)
max(trial1$time.total)

trial2 <- simulate_trial(control, active, events=150)
head(trial2, 5)
sum(trial2$censor)

## ------------------------------------------------------------------------
trial3 <- simulate_trial(control, active, duration=18, events=150)
max(trial3$time.total)
sum(trial3$censor)

## ------------------------------------------------------------------------
control_sim <- simulate_arm(control)
head(control_sim, 5)
table(control_sim$arm, control_sim$reason)

