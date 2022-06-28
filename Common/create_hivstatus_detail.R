
library(tidyverse); library(survival)

# This R script creates a variable hivstatus_detail that splits a person's person time
# into the most detailed categories.  It relies on you having already created the following
# variables in your analysis file:
#   
#   timepostneg - This is the cut off you want to you for post negative time
# timeprepos - This is the cut off you want for pre positive time
# 
# Also HIV test data variables based on hiv_tests data
# 
# first_neg_date - First negative test date
# last_neg_date - Last negative test date
# first_pos_date - First positive test date
# last_pos_date - First positive test date
# 
# 
# sero_conv_date - This is the date halfway between last negative date and first positive date
# 
# 
# The global value prefix is so that if you have more that one set of hivstatus details to create in one dataset
# you can use the prefix to differentiate them.  For example a file of children with date also on mothers and fathers
# you might want hivstatus for the child and then split the child's experience into the hivstatus_detail of the parents

# Generate a date that is X years after the last negative test
dat <- dat %>%
    dplyr::mutate(after_negXyear = !!rlang::parse_expr(paste0(prefix,"last_neg_date")) +
                      (365.25*timepostneg))

# Generate a date that is X years before the first positive test
dat <- dat %>%
    dplyr::mutate(before_posXyear = !!rlang::parse_expr(paste0(prefix,"first_pos_date")) - 
                      (365.25*timeprepos))

# keep track of failure date since tmerge doesn't account for this
dat <- dat %>% 
    dplyr::mutate(failuredate = as.Date(ifelse(failure==1,exit,NA),origin="1970-01-01"))

# Splitting up follow-up time by HIV status as appropriate

# Negative splits
# for the purpose of the tmerge function, treating each survival episode as its own
# id
dat <- dat %>% dplyr::group_by(idno) %>%
    dplyr::mutate(idno_ep = paste0(idno,".",row_number()))

dat <- survival::tmerge(dat,dat, id=idno_ep,
                         tstart = entry, tstop = exit, split_on_firstneg = tdc(first_neg_date))
dat <- survival::tmerge(dat,dat,id = idno_ep,
                         split_on_lastneg = tdc(last_neg_date))
dat <- survival::tmerge(dat,dat,id=idno_ep,
                         split_on_afternegXyear = tdc(after_negXyear))

# Positive splits
dat <- survival::tmerge(dat, dat, id = idno_ep,
                            split_on_firstpos = tdc(first_pos_date))
dat <- survival::tmerge(dat, dat, id = idno_ep,
                            split_on_beforeposXyear = tdc(before_posXyear))
dat <- survival::tmerge(dat, dat, id = idno_ep,
                            split_on_lastpos = tdc(last_pos_date))

# Sero conversion split
if(sum(!is.na(dat$sero_conv_date))>0) {
    dat <- survival::tmerge(dat, dat, id = idno_ep,
                            split_at_sc = tdc(sero_conv_date))
} else {
    dat$split_at_sc <- 0
}



# Negative time
# Negative
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_firstneg==1 & split_on_lastneg==0, 1, NA))
# before first negative test
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_firstneg==0 & !is.na(first_neg_date), 3, hivstatus_detail))
# post negative within cutoff
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_lastneg==1 & split_on_afternegXyear==0, 4, hivstatus_detail))
# post negative beyond cutoff
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_lastneg==1 & split_on_afternegXyear==1, 5, hivstatus_detail))

# Positive time
# Positive
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_firstpos==1 & split_on_lastpos==0, 2, hivstatus_detail))
# after last positive test
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_lastpos==1, 6, hivstatus_detail))
# pre positive within cutoff
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_firstpos==0 & split_on_beforeposXyear==1 & !is.na(before_posXyear) 
                                                & is.na(first_neg_date), 7, hivstatus_detail))
# pre positive beyond cutoff
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_firstpos==0 & split_on_beforeposXyear==0 & 
                                                            !is.na(first_pos_date) & is.na(first_neg_date), 8, hivstatus_detail))

# Sero conversion interval
# SC interval post neg within cutoff
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_lastneg == 1 & split_at_sc==0 & split_on_afternegXyear==0 & !is.na(sero_conv_date), 9, hivstatus_detail))
# SC interval post neg beyond cutoff
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_lastneg==1 & split_at_sc==0 & split_on_afternegXyear==1 & !is.na(sero_conv_date), 10, hivstatus_detail))
# SC interval pre pos within cutoff
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_firstpos==0 & split_at_sc==1 & split_on_beforeposXyear==1, 11, hivstatus_detail))
# SC interval pre pos beyond cutoff
dat <- dat %>% mutate(hivstatus_detail = ifelse(split_on_firstpos==0 & split_at_sc==1 & split_on_beforeposXyear==0, 12, hivstatus_detail))

# Unknown
dat <- dat %>% mutate(hivstatus_detail = ifelse(is.na(first_neg_date) & is.na(first_pos_date), 13, hivstatus_detail))

# correct failure date
dat <- dat %>%
    dplyr::mutate(failure = ifelse(tstop==failuredate & !is.na(failuredate),1,0))

# Drop all variables that were created to make hivstatus_detail
dat <- dat %>% 
    dplyr::select(-c(after_negXyear, before_posXyear, split_on_firstneg, split_on_lastneg,
                     split_on_afternegXyear, split_on_firstpos, split_on_beforeposXyear,
                     split_on_lastpos, split_at_sc,failuredate))


