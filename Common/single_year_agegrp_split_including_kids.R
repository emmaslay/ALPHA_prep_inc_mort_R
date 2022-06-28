
library(survival); library(tidyverse)

################ SPLIT DATA INTO SINGLE YEARS AND CREATE AGE ###################

# keep dob to remerge since disappears after survSplit call
dat_dob <- dat %>% select(idno,dob) %>% group_by(idno) %>%
    filter(row_number()==1)

dat <- survSplit(Surv(as.numeric(tstart-dob),as.numeric(tstop-dob),failure) ~ .,
                   dat, cut = seq(0,90*365.25,365.25), event = "failure")
dat$age <- floor(dat$tstart/365.25)
# assign everyone over age 90 to 90
dat$age <- ifelse(dat$age>90,90,dat$age)
# filter out anyone with age under 0
dat <- dat %>% dplyr::filter(age>=0)
# merge dob back in
dat <- merge(dat,dat_dob,by="idno",all.x=TRUE)
# convert tstart and tstop back into dates
dat <- dat %>%
    dplyr::mutate(tstart = tstart + dob,
                  tstop = tstop + dob)


