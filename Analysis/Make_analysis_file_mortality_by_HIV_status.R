library(tidyverse); library(labelled); library(haven)

# MAKE ALPHA MORTALITY BY HIV STATUS DATASET
# IN THE SPIRIT OF THE GATES MORTALITY ANALYSIS
# DON'T BRING IN CARE CONT VARS HERE BECAUSE NEED TO BRING IN MILLY AND CLARA'S NEW WAY OF DOING THOSE
# Starts from the residency and HIV tests specs, makes an analysis file for each site and then pools them.



################################################################################
#############  3. Merges, retain relevant records and prep for survival analysis
################################################################################


################################################################################
##             3.1 MERGE SPECS TOGETHER AND TEMPORARY SURVIVAL PREP TO MANIPULATE RECORDS
################################################################################
# merge all specs to residency, keeping only those who also have residency information 
# (appearing on 6.1)

dat <- haven::read_dta(paste0(alphapath,
                                      "/ALPHA/clean_data/",sitename,
                                      "/residency_",sitename,".dta"))

    ## For Manicaland, drop out the communities that weren't included in R6
    if(tolower(sitename)=="manicaland"){
        dat$community <- as.integer(dat$hhold_id/10000)
        dat <- dat %>% dplyr::filter(community!=1 | community!=6 | community!=11 |
                              community!=12) %>%
            select(-community)
    }
    ## For uMkhanyakude drop TasP people
    if(tolower(sitename)=="umkhanyakude") {
        dat <- dat %>% dplyr::filter(entry_type!=1 & format(entry_date,"%Y")!="2017")
    }
    ## For Kisumu, keep only Gem
    if(tolower(sitename)=="kisumu") {
        dat <- dat %>% dplyr::filter(residence==2)
    }

# HIV testing
wide_tests <- readRDS(paste0(alphapath,"/ALPHA/prepared_data/",sitename,
                             "/hiv_tests_wide_",sitename,".rds"))
wide_tests <- wide_tests %>%
    dplyr::select(study_name, idno, first_neg_date, last_neg_date, first_pos_date, 
                  last_pos_date, lasttest_date, firsttest_date)
dat <- merge(dat,wide_tests,by=c("study_name","idno"),all.x=TRUE)

# metadata
metadat <- haven::read_dta(paste0(alphapath,
                                          "/ALPHA/clean_data/alpha_metadata.dta"))
metadat <- metadat %>% dplyr::select(study_name,earlyexit_max,timepostneg)
dat <- merge(dat,metadat,
             by="study_name", all.x=TRUE)

# dates in stata have origin 1960-01-01
if(class(dat$entry_date)=="numeric") {
    dat$entry_date <- as.Date(dat$entry_date, origin="1960-01-01")
}
if(class(dat$exit_date)=="numeric") {
    dat$exit_date <- as.Date(dat$exit_date, origin="1960-01-01")
}
if(class(dat$dob)=="numeric") {
    dat$dob <- as.Date(dat$dob,origin="1960-01-01")
}

dat <- dat %>%
    dplyr::mutate(failure = ifelse(exit_type==2,1,0),
                  exit = exit_date,
                  entry = entry_date) 

# Deal w/data errors for survival data
# Take out records who have missing entry/exit times 
# missing_entry_exit <- dat %>% dplyr::filter(is.na(exit) | is.na(entry))
dat <- dat %>% dplyr::filter(!is.na(exit) & !is.na(entry))

# Take out those who have same day exit and entry
# same_day_enter_exit <- dat %>% dplyr::filter(exit==entry)
dat <- dat %>% dplyr::filter(exit>entry)

# Take out overlapping records (exit[n-1]>entry)
# overlap <- dat %>% dplyr::group_by(idno) %>%
#     arrange(entry) %>%
#     dplyr::mutate(overlap = ifelse(lag(exit)>entry,1,0)) %>%
#     filter(overlap==1)
dat <- dat %>% dplyr::group_by(idno) %>%
    arrange(entry) %>%
    filter(lag(exit)<=entry | is.na(lag(exit)))

# Take out those who are missing dob
# nodob <- dat %>% filter(is.na(dob))
dat <- dat %>% dplyr::filter(!is.na(dob))

# Take out records whose DOB is on or later than exit date
# dob_after_exit <- dat %>% filter(dob>=exit)
dat <- dat %>% dplyr::filter(dob<exit)

# Take out records that begin on or after death
dat <- dat %>% group_by(idno) %>%
    mutate(death_date = ifelse(failure==1,exit,NA),
           first_death_date = as.Date(ifelse(is.finite(min(death_date,na.rm=TRUE)),
                                             min(death_date, na.rm=TRUE),NA),origin="1970-01-01"),
           on_after_death = ifelse(entry>=first_death_date,1,0)) %>%
    filter(is.na(on_after_death) | on_after_death!=1) %>%
    select(-c(death_date,first_death_date,on_after_death))

################################################################################
##      3.2. DROP RECORDS BEFORE AGE OF 15                                    ##
################################################################################
# split at age 15 and drop all episodes before this age
#create a variable for date turned 15
dat <- dat %>%
    dplyr::mutate(d15 = dob + 365.25*15)

# for the purposes of tmerge treat each idno_ep as a unique id
dat <- dat %>% dplyr::group_by(idno) %>%
    dplyr::mutate(idno_ep = paste0(idno,".",row_number()))

dat <- survival::tmerge(dat,dat, id=idno_ep,
                        tstart = entry, tstop = exit, pre15 = tdc(d15))

dat <- dat %>% dplyr::filter(pre15==1) %>%
    dplyr::select(-pre15,-d15)


# update entry & exit since tmerge 
dat <- dat %>% 
    dplyr::mutate(entry = tstart, 
                  exit = tstop)

################################################################################
##      3.3. CHECK CONSISTENCY OF START/END DATES ACROSS SPECS AND CORRECT    ##
################################################################################

## EARLY EXIT ISSUES
dat <- dat %>% dplyr::group_by(study_name, idno) %>%
    dplyr::arrange(entry) %>%
    dplyr::mutate(episode_sequence = row_number(),
                  episode_total = n(),
                  last_episode = ifelse(episode_sequence==episode_total,1,0),
                  last_exit_original = exit[last_episode==1])

# if exit is on the same date as the last test, move the exit date to one day after
dat <- dat %>% 
    dplyr::mutate(early_exit_fixed_t = ifelse(as.character(exit)==as.character(lasttest_date) & 
                                                  !is.na(lasttest_date) &
                                                  last_episode==1, 1, NA),
                  exit = as.Date(ifelse(as.character(exit)==as.character(lasttest_date) & 
                                            !is.na(lasttest_date) & last_episode==1,
                                        exit+1, exit),origin="1970-01-01"))

# identify people whose latest 6.1 exit is before latest 6.2, 9.1 or 9.2 report & calculate difference
# exit before last test
dat <- dat %>%
    dplyr::mutate(exitgap_test = ifelse(exit<lasttest_date & !is.na(lasttest_date) & last_episode==1
                                        & !is.na(exit),
                                        lasttest_date - exit,NA),
                  early_exit_problem = ifelse(!is.na(exitgap_test),1,NA))

# change exit to one day after last test, SR or clinic report [all exit types for 
# first 2, only if not dead or out-migrated for clinic data as they could move out
# but still go to same clinic]
# May 2015, stop moving exit to after last clinic date as biased towards positive now we have more data
dat <- dat %>%
    dplyr::mutate(early_exit_fixed_t = case_when(!is.na(early_exit_fixed_t) ~ as.integer(early_exit_fixed_t),
                                                 exitgap_test <= earlyexit_max ~ 1L,
                                                 !is.na(early_exit_problem) & is.na(early_exit_fixed_t) ~ 0L,
                                                 TRUE ~ NA_integer_),
                  exit_new = ifelse(early_exit_fixed_t==1,lasttest_date + 1, NA))

dat <- dat %>% dplyr::group_by(study_name,idno) %>%
    dplyr::mutate(early_exit_fixed = max(early_exit_fixed_t,na.rm=TRUE))
dat$early_exit_fixed[!is.finite(dat$early_exit_fixed)] <- NA

dat <- dat %>%
    dplyr::mutate(exit = as.Date(ifelse(!is.na(exit_new), exit_new, exit),
                                 origin="1970-01-01"))

dat <- dat %>%
    dplyr::mutate(tstart = as.Date(as.character(entry)),
                  tstop = as.Date(as.character(exit)))


################################################################################
##      4. SPLIT at Calendar years                                            ##
################################################################################

source(paste0(alphapath,"/common/Calendar_year_split.R"))
source(paste0(alphapath,"/common/create_fouryear.R"))

dat <- dat %>% dplyr::filter(years_one>=1989)

# convert tstart and tstop back to dates (origin in R is 1970-01-01)
dat <- dat %>%
    dplyr::mutate(tstart = as.Date(tstart,origin="1970-01-01"),
                  tstop = as.Date(tstop,origin="1970-01-01"))

################################################################################
##      5. SPLIT AT SINGLE YEARS AND MAKE 5 YEAR AGEGRP                       ##
################################################################################

source(paste0(alphapath,"/common/single_year_agegrp_split_including_kids.R"))
dat <- dat %>% dplyr::filter(age>=15)
source(paste0(alphapath,"/common/Create_agegrp_from_age.R"))

# update exit & entry after splits
dat <- dat %>%
    dplyr::mutate(entry = tstart,
                  exit = tstop)


################################################################################
##      6. ASSIGN HIV STATUS                                                  ##
################################################################################

dat <- dat %>%
    dplyr::mutate(first_neg_date = as.Date(first_neg_date,origin="1970-01-01"),
                  last_neg_date = as.Date(last_neg_date,origin="1970-01-01"),
                  first_pos_date = as.Date(first_pos_date,origin="1970-01-01"),
                  last_pos_date = as.Date(last_pos_date,origin="1970-01-01"))

# # post negative times from cubic spline
# postneg<- haven::read_dta(paste0(alphapath,
#                                          "/ALPHA/Estimates_Incidence/Post_negative_times/post_neg_ages_",sitename,
#                                          ".dta"))
# dat <- merge(select(dat,-timepostneg),select(postneg,-study_name),by=c("sex","age","fouryear","agegrp"),all.x=TRUE)
# 
# # some aren't linked because that fouryear/age combination doesn't exist
# # this is usually because there is no HIV data, to solve this carry forwards estimates from earlier period
# # for sites which have stopped testing, copy forwards the old estimates
# postneg <- postneg %>%
#     dplyr::mutate(fouryear_forwards1 = fouryear + 1,
#                   fouryear_forwards2 = fouryear + 2)
# 
# dat <- merge(dat,select(postneg,sex,age,fouryear_forwards1,timepostneg),
#               by.x=c("sex","age","fouryear"),by.y=c("sex","age","fouryear_forwards1"),
#               all.x=TRUE,suffixes=c("","_offset"))
# 
# dat <- dat %>%
#     dplyr::mutate(timepostneg = ifelse(is.na(timepostneg),timepostneg_offset,timepostneg)) %>%
#     dplyr::select(-timepostneg_offset)
# 
# # do this for offset2 as well?
# # dat <- merge(dat,select(postneg,sex,age,fouryear_forwards2,timepostneg),
# #              by.x=c("sex","age","fouryear"),by.y=c("sex","age","fouryear_forwards2"),
# #              all.x=TRUE,suffixes=c("","_offset"))
# #
# # dat <- dat %>%
# #     dplyr::mutate(timepostneg = ifelse(is.na(timepostneg),timepostneg_offset,timepostneg)) %>%
# #     dplyr::select(-timepostneg_offset)

# For now using timepostneg from metadata






## HIV STATUS based on 6.2b data
# need a value here for pre-positive time so the do file can run, but will later discard
# all the pre-positive time and make it unknown
dat <- dat %>%
    dplyr::mutate(sero_conv_date = as.Date(ifelse(is.finite(last_neg_date) & is.finite(first_pos_date),
                                                  as.numeric(last_neg_date) + 
                                                      ((as.numeric(first_pos_date) - as.numeric(last_neg_date))/2), NA), origin="1970-01-01"),
                  timeprepos = 1)
# add 6 months to the time post-negative test to allow enough time for the next interview
dat$timepostneg <- dat$timepostneg + 0.5

source(paste0(alphapath,"/Common/create_hivstatus_detail.R"))

# time before a first test and after the cutoff after a negative test is unknown,
# in the seroconversion interval the years are allocated to negative up to the
# timepostneg cutoff then to unknown until the first positive test

# 1 "Negative" => 1
# 2 "Positive"  => 2
#
# 3 "Before first negative test" => 3
# 4 "Post Negative within cutoff" =>1
# 5 "Post negative beyond cutoff" => 3
#
# 6 "After last positive test" =>2
# 7 "Pre positive within cutoff" =>3
# 8 "pre positive beyond cutoff" =>3
#
# 9 "SC interval post neg within cutoff" =>1
# 10 "SC interval post neg beyond cutoff" =>3
# 11 "SC interval pre positive within cutoff" =>3 
# 12 "SC interval pre positive beyond cutoff" =>3
#
# 13 "Unknown, never tested" => 3

dat <- dat %>%
    dplyr::mutate(hivstatus_broad = case_when(hivstatus_detail %in% c(1,4,9) ~ 1,
                                              hivstatus_detail %in% c(2,6) ~ 2,
                                              hivstatus_detail %in% c(3,5,7,8,10,11,12,13) ~ 3,
                                              TRUE ~ NA_real_))

################################################################################
##      7. SAVE SITE SPECIFIC RESULTS READY FILES                             ##
################################################################################

## clean up variables for outdat
dat <- dat %>%
    dplyr::mutate(entry = tstart,
                  exit = tstop,
                  t0 = as.numeric((tstart-dob)/365.25),
                  t = as.numeric((tstop - dob)/365.25),
                  d = failure,
                  st = 1,
                  origin = dob) %>%
    select(-any_of(c("tstart","tstop","above_max_age","above95","age_last_test","failure",
           "id", "idno_ep", "no_haz_estimate","postnegage", "str_study")))

# Label variables
labelled::var_label(dat) <- list(d = "1 if failure; 0 if censored",
                                  origin = "Evaluated value of origin()",
                                  st = "1 if record is to be used; 0 otherwise",
                                  t = "Analysis time when record ends",
                                  t0 = "Analysis time when record begins",
                                  age = "Age in single years",
                                  agegrp = "Five year age group",
                                  earlyexit_max = "Acceptable interval between end of residency episode and death/HIV test., days",
                                  entry_date = "1 event_date",
                                  entry_type = "1 event_type",
                                  exit_date = "2 event_date",
                                  exit_type = "2 event_type",
                                  fouryear = "Calendar year, grouped in 4 years post 2005",
                                  hivstatus_detail = "HIV status detail",
                                  years_one = "Calendar year")

dat$agegrp <- as.numeric(as.character(dat$agegrp))
labelled::val_labels(dat$agegrp) <- c("0-4"=0,"5-9"=1,"10-14"=2,"15-19"=3,"20-24"=4,"25-29"=5, 
                                     "30-34"=6, "35-39"=7, "40-44"=8, "45-49"=9,  
                                     "50-54"=10, "55-59"=11, "60-64"=12, "65-69"=13,
                                     "70-74"=14, "75-79"=15, "80-84"=16, "85-89"=17,"90+"=18)
labelled::val_labels(dat$early_exit_fixed) <- c("No change to exit" = 0, "Exit changed to last 6.2 plus 1 day"=1,
                                               "Exit changed to last 9.1 plus 1 day" = 2,
                                               "Exit changed to last 9.2 plus 1 day" = 3)
labelled::val_labels(dat$early_exit_problem) <- c("Exit<last6.2"=1,"Exit<last9.1"=2,
                                                  "Exit<last9.2"=3,"Exit<last6.2&9.1"=4,
                                                  "Exit<last6.2&9.2"=5,"Exit<last 9.1&9.2"=6,
                                                  "Exit<last6.2&9.1&9.2"=7)
dat$fouryear <- as.numeric(as.character(dat$fouryear))
labelled::val_labels(dat$fouryear) <- c("earliest-1999"=0,"2000-04"=1,"2004-08"=2,"2009-12"=3,
                                        "2013-16"=4,"2017-20"=5,"2021-24"=6)
labelled::val_labels(dat$hivstatus_broad) <- c("Negative" = 1, "Positive" = 2, "Unknown" = 3)
labelled::val_labels(dat$hivstatus_detail) <- c("Negative" = 1,"Positive"  = 2,
                                                "Before first negative test" = 3,
                                                "Post Negative within cutoff" = 4,
                                                "Post negative beyond cutoff" = 5,
                                                "After last positive test" = 6,
                                                "Pre positive within cutoff" = 7,
                                                "pre positive beyond cutoff" = 8,
                                                "SC interval post neg within cutoff" = 9,
                                                "SC interval post neg beyond cutoff" = 10,
                                                "SC interval pre positive within cutoff" = 11, 
                                                "SC interval pre positive beyond cutoff" = 12,
                                                "Unknown, never tested" = 13)
labelled::val_labels(dat$sex) <- c("Men" = 1, "Women" =2)
labelled::var_label(dat) <- list(agegrp = "Five year age group",
                                 fouryear = "Calendar year, grouped in 4 years post 2005")




if(!dir.exists(paste0(alphapath,"/ALPHA/Ready_data_mortality/",
                      sitename))) {
    dir.create(paste0(alphapath,"/ALPHA/Ready_data_mortality/",
                      sitename),recursive=TRUE)
}
saveRDS(dat,paste0(alphapath,"/ALPHA/Ready_data_mortality/",
                   sitename,"/mortality_by_status_ready_",sitename,".rds"))





