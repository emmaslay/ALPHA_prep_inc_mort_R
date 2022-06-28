
library(survival); library(tidyverse)

# Splits person time, in an already stset dataset, into age in FIVE YEAR GROUPS- creates 
# agegrp variable. Don't use this if already split on age, use create_agegrp_from_age.

################ SPLIT DATA INTO FIVE YEAR AGE GROUPS AND CREATE AGEGRP ###################

# keep dob to remerge since disappears after survSplit call
dat_dob <- dat %>% select(idno,dob) %>% group_by(idno) %>%
    filter(row_number()==1)

dat <- survSplit(Surv(as.numeric(tstart-dob),as.numeric(tstop-dob),failure) ~ .,
                 dat, cut = seq(0,90*365.25*5,365.25*5), event = "failure")
dat$age <- floor(dat$tstart/365.25)
dat <- dat %>%
    dplyr::mutate(age = floor(tstart/365.25),
                  agegrp = cut(age, c(seq(0,90,5),150), include.lowest=TRUE, right=FALSE,
                               labels = 0:18)) %>%
    dplyr::select(-age)
# merge dob back in
dat <- merge(dat,dat_dob,by="idno",all.x=TRUE)
# convert tstart and tstop back into dates
dat <- dat %>%
    dplyr::mutate(tstart = tstart + dob,
                  tstop = tstop + dob)

# Label agegrp
dat$agegrp <- as.numeric(as.character(dat$agegrp))
labelled::val_labels(dat$agegrp) <- c("0-4"=0,"5-9"=1,"10-14"=2,"15-19"=3,"20-24"=4,"25-29"=5, 
                                      "30-34"=6, "35-39"=7, "40-44"=8, "45-49"=9,  
                                      "50-54"=10, "55-59"=11, "60-64"=12, "65-69"=13,
                                      "70-74"=14, "75-79"=15, "80-84"=16, "85-89"=17,"90+"=18)

labelled::var_label(dat) <- list(agegrp = "Five year age group")

