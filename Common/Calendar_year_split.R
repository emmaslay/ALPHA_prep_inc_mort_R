
library(survival)

firstyear <- format(min(dat$tstart),"%Y")
lastyear <- format(max(dat$tstop),"%Y")

cut_times <- as.numeric(as.Date(paste0(firstyear:lastyear,"-01-01")))

dat <- survival::survSplit(Surv(as.numeric(tstart),as.numeric(tstop),failure) ~ .,
                  dat, cut = cut_times, episode = "years_one", event = "failure")

dat$years_one <- as.numeric(format(as.Date(dat$tstart,origin="1970-01-01"),"%Y"))

dat$tstart <- as.Date(dat$tstart,origin="1970-01-01")
dat$tstop <- as.Date(dat$tstop, origin="1970-01-01")

