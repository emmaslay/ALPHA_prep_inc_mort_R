##########################################################################
### ALPHA BASIC INCIDENCE ANALYSIS USING MIDPOINT AS SEROCONVERSION DATE
##########################################################################

# this is the simplest incidence R script and is the starting point for incidence 
# analyses, most of which will additionally use the multiply imputed seroconversion dates
# and the behavioural and socio-demographic data 

library(tidyverse); library(survival); library(haven)


# need an estimate for number of years after 15th birthday we make the adjustment for prevalent positives
# this could vary between sites but in the workshop everyone used 3
maxinter <- 2


    ## START WITH RESIDENCY EPISODES
    dat <- haven::read_dta(paste0(alphapath,
                                          "/ALPHA/clean_data/",sitename,
                                          "/residency_",sitename,".dta"))

    # For Manicaland, drop out the communities that weren't included in R6
    if(tolower(sitename)=="manicaland") {
        dat$community <- floor(dat$hhold_id/10000)
        dat <- dat %>% 
            dplyr::filter(!(community==1 | community==6 | community==11 | community==12)) %>%
            dplyr::select(-community)
    }
    ## For uMkhanyakude drop TasP people
    if(tolower(sitename)=="umkhanyakude") {
        dat <- dat %>% dplyr::filter(!(entry_type==1 & format(entry_date,"%Y")=="2017"))
    }
    
    # MERGE IN TEST DATES IN WIDE FORM, ONLY USING TESTS DONE WHEN PERSON WAS RESIDENT
    wide_tests <- readRDS(paste0(alphapath,"/ALPHA/prepared_data/",sitename,
                                 "/hiv_tests_wide_",sitename,".rds"))
    wide_tests <- wide_tests %>%
        dplyr::select(study_name, idno, first_neg_date, last_neg_date, first_pos_date, 
                      last_pos_date, lasttest_date, firsttest_date)
    dat <- merge(dat,wide_tests,by=c("study_name","idno"),all.x=TRUE)
    
    # merging in metadata
    dat <- merge(dat,haven::read_dta(paste0(alphapath,
                                                    "/ALPHA/clean_data/alpha_metadata.dta")),
                 by="study_name", all.x=TRUE)
    
    ## For Kisumu, keep only Gem
    if(tolower(sitename)=="kisumu") {
        dat <- dat %>% dplyr::filter(residence==2)
    }
    
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
    
    # Survival data - death (exit type 2) is failure
    # setting on mortality temporarily to organise the data - want to be able to use survSplit
    dat <- dat %>%
        dplyr::mutate(exit = exit_date,
                      entry = entry_date,
                      failure = ifelse(exit_type==2,1,0))
    
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
    
    
    
    ## EARLY EXIT ISSUES- these arise in sites where the DSS interview and the HIV test don't always take place on the same day.  It is necessary to move the exit date to after
	# the test dates, but there is an upper limit (defined in the metadata) beyond which the exit date shouldn't be moved.
    # The upper limit depends on what is known about fieldwork and how long the lag between DSS and HIV test is likely to have been.
    # If the tests are beyond the upper limit they are discarded if this is the last episode
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
    
    # identify people whose latest residency exit is before latest HIV test report & calculate difference
    # exit before last test
    dat <- dat %>%
        dplyr::mutate(exitgap_test = ifelse(exit<lasttest_date & !is.na(lasttest_date) & last_episode==1
                                            & !is.na(exit),
                                            lasttest_date - exit,NA),
                      early_exit_problem = ifelse(!is.na(exitgap_test),1,NA))
    
    # change exit to one day after last test, SR or clinic report [all exit types for 
    # first 2, only if not dead or out-migrated for clinic data as they could move out
    # but still go to same clinic]
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
    
    
    ##########################################################################
    ##  3. SPLIT USING 6.2b DATA (HIV TESTS) & CREATE HIV STATUS VARIABLE   ##
    ##########################################################################
    
    dat <- dat %>%
        dplyr::mutate(sero_conv_date = as.Date(ifelse(is.finite(last_neg_date) & is.finite(first_pos_date),
                                              as.numeric(last_neg_date) + 
                                                  ((as.numeric(first_pos_date) - as.numeric(last_neg_date))/2), NA), origin="1970-01-01"),
                      timeprepos = 0)
    
    ####### RIGHT NOW TIMEPOSTNEG IS BROUGHT IN FROM THE METADATA - IF CHANGE TO
    ####### USING THE VERSION PREPPED BY EMMA, BRING IN HERE
    
    source(paste0(alphapath,"/Common/create_hivstatus_detail.R"))
    
    # fix entry & exit date since tmerge 
    dat <- dat %>%
        dplyr::mutate(entry = tstart,
                      exit = tstop)
    
    ############### REORGANISE FOR INCIDENCE ANALYSIS ###############

    ## define date of 15th birthday
    dat$fifteen <- dat$dob + (15*365.25)
    
    # OPTION TO INCLUDE YOUNGEST RESIDENTS COMING TO FIRST SERO WHO WERE
    # RESIDENT AT EARLIER SERO BUT NOT AGE ELIGIBLE ASSUMING THEY WERE
    # NEGATIVE AT AGE 15
    # COMMENT OUT THIS BIT IF THIS OPTION NOT NEEDED
    
    maxrounds <- mean(dat$nrounds_sero)
    
    dat <- dat %>% 
        dplyr::mutate(tested_first_opportunity = ifelse(firsttest_date < (fifteen + maxinter*365.25) & 
                                                    firsttest_date > fifteen & !is.na(firsttest_date), 1, 0))
    
    # check if resident at 15
    dat <- dat %>% dplyr::group_by(idno) %>%
        dplyr::mutate(resat15 = any(fifteen>tstart & fifteen<tstop))
    
    # RETAIN RECORDS RELEVANT FOR INCIDENCE ANALYSIS
    # BETWEEN TWO NEGATIVE TESTS, AFTER A NEGATIVE TEST BUT WITHIN CUTOFF, IN THE SEROCONVERSION INTERVAL
    dat <- dat %>% 
        dplyr::mutate(tokeep = ifelse(hivstatus_detail==1 | hivstatus_detail==4 |
                                          hivstatus_detail==9 | hivstatus_detail==10 |
                                          hivstatus_detail==11 | hivstatus_detail==12 |
                                          # also keep the time of young people prior to their first test 
                                          # if they 1) were too young to have been tested in the previous round and 
                                          # 2) had their first test at the first opportunity
                                          (tested_first_opportunity==1 & resat15 & tstop<=firsttest_date), 1, NA))
    
    
    dat <- dat %>% dplyr::filter(tokeep==1) %>%
        dplyr::select(-tokeep)
    
    # for young prevalent positives, impute a seroconversion date midway between their
    # first test date and fifteenth birthday
    dat <- dat %>%
        dplyr::mutate(sero_conv_date = as.Date(ifelse(tested_first_opportunity==1 & resat15 &
                                                  tstop<=firsttest_date & firsttest_date==first_pos_date,
                                              ((first_pos_date - fifteen)*0.5) + fifteen, sero_conv_date),origin="1970-01-01"))
    
    
    dat <- dat %>% dplyr::select(idno, sex, dob, residence, entry, entry_type, exit, exit_type, 
                                 study_name, first_neg_date, last_neg_date, first_pos_date,
                                 last_pos_date, firsttest_date, lasttest_date, failure, 
                                 tstop, tstart, hivstatus_detail, nrounds_sero, fifteen, tested_first_opportunity, 
                                 sero_conv_date)
    
    #################################################################################
    #  SPLIT THE DATA BY AGE AND CALENDAR YEAR FOR AGE AND YEAR SPECIFIC ESTIMATES  #
    #################################################################################
    
    # split the data by age group and calendar year
    source(paste0(alphapath,"/common/single_year_agegrp_split_including_kids.R"))
    
    source(paste0(alphapath,"/common/Calendar_year_split.R"))
    
    # create new categorical variables
    dat <- dat %>% dplyr::filter(age>=15)
    source(paste0(alphapath,"/common/Create_agegrp_from_age.R"))
    source(paste0(alphapath,"/common/Create_birth_cohort_from_dob.R"))
    
    source(paste0(alphapath,"/common/create_fiveyear.R"))
    source(paste0(alphapath,"/common/create_fouryear.R"))
    
    # update entry and exit after tmerge calls
    dat <- dat %>%
        dplyr::mutate(entry = tstart,
                      exit = tstop)
    
    
    
    #################################################################################
    #   SET UP FOR ANALYSIS WITH MIDPOINT AS SEROCONVERSION                         #
    #################################################################################
    
    # censor everyone at the most recent test date
    dat <- dat %>% dplyr::group_by(idno) %>%
        dplyr::mutate(idno_ep = paste0(idno,".",row_number()))
    dat <- survival::tmerge(dat,dat, id=idno_ep,
                             tstart = entry, tstop = exit,
                             afterlasttest = tdc(lasttest_date))
    
    dat <- dat %>% dplyr::filter(afterlasttest==0) %>% 
        dplyr::select(-afterlasttest)
    
    ## recreate entry and exit variables - old dates no longer valid after splits 
    dat$start_ep_date <- dat$tstart
    dat$end_ep_date <- dat$tstop
    dat$entry <- dat$tstart
    dat$exit <- dat$tstop
    dat$died <- dat$failure
    dat$t0 <- as.numeric((dat$start_ep_date - dat$dob)/365.25)
    dat$t <- as.numeric((dat$end_ep_date - dat$dob)/365.25)
    
    # CAN'T USE A FILE SURVIVAL SET ON INCIDENCE AS THE BASIS FOR MI ESTIMATES OF INCIDENCE RISK FACTORS. THIS IS BECAUSE 
	# THE DATASET REQUIRES ADDITIONAL SPLITTING TO PREPARE THE DATA FOR RISK FACTOR ANALYSIS. WE NEED ALL RECORDS TO BE SPLIT AS APPROPRIATE BUT
	# IF THE DATASET IS SET ON SERCONVERSION, WITH FAILURE AT THE MIDPOINT, THE STSET WILL EXCLUDE ANY SEROCONVERTORS WHO ARE NOT RESIDENT AT THE MIDPOINT
	# THIS MEANS THEY WOULD NOT BE PROPERLY INCLUDED IN THE RISK FACTOR ANALYSIS.  TO GET AROUND THIS PROBLEM, WE CREATE THE RISK FACTOR FILE FROM ONE
	# THAT IS SET ON MORTALITY SO THAT EVERYONE IS INCLUDED.
    if(!dir.exists(paste0(alphapath,"/ALPHA/incidence_ready_data/",
                          sitename))) {
        dir.create(paste0(alphapath,"/ALPHA/incidence_ready_data/",
                          sitename),recursive=TRUE)
    }
    saveRDS(dat, paste0(alphapath,"/ALPHA/incidence_ready_data/",sitename,
                        "/incidence_temp_for_risk_factors_",sitename,".rds"))
    
    # NOW SET THE DATA FOR INCIDENCE ANALYSIS USING THE MIDPOINT AS THE SEROCONVERSION DATE
    dat$sero_conv_date <- as.Date(dat$sero_conv_date,origin="1970-01-01")
    # count the number of people observed to seroconvert - not necessarily resident on this date
    n_distinct(dat$idno[!is.na(dat$first_pos_date) & !is.na(dat$last_neg_date) & 
                            dat$first_pos_date>dat$last_neg_date])
    
    # make a variable to indicate records that contain a seroconversion
    dat <- dat %>%
        dplyr::mutate(serocon_fail = ifelse(sero_conv_date>start_ep_date & 
                                                sero_conv_date<=end_ep_date & !is.na(sero_conv_date),1,0))
    
    # Move the end date of episodes that contain a seroconversion - the episode will now end at the time of seroconversion
    dat$end_ep_date <- as.Date(ifelse(dat$serocon_fail==1, dat$sero_conv_date, dat$end_ep_date),origin="1970-01-01")
    dat$t <- as.numeric((dat$end_ep_date - dat$dob)/365.25)
    
    # Clean up the data for survival analysis
    # Remove dates that end on or before enter date
    dat <- dat %>% 
        dplyr::filter(start_ep_date<end_ep_date)
    # If using the midpoint file for analyses, need to remove observations that 
    # start on or after seroconversion date
    # dat <- dat %>% 
    #     dplyr::filter(start_ep_date<sero_conv_date | is.na(sero_conv_date))

    ## clean up variables for outdat
    dat <- dat %>%
        dplyr::mutate(entry = start_ep_date,
                      exit = end_ep_date,
                      t0 = as.numeric((start_ep_date-dob)/365.25),
                      t = as.numeric((end_ep_date - dob)/365.25),
                      d = serocon_fail,
                      st = 1,
                      origin = dob) %>%
        select(-any_of(c("tstart","tstop","above_max_age","above95","age_last_test",
                      "id", "idno_ep", "no_haz_estimate","postnegage", "str_study")))
    
    # Label variables
    labelled::var_label(dat) <- list(age = "Age in single years",
                                     agegrp = "Five year age group",
                                     birth_cohort = "Five year birth_cohort",
                                     end_ep_date = "Date this record ends on, for everyone included in midpoint file",
                                     entry_type = "1 event_type",
                                     exit_type = "2 event_type",
                                     fifteen = "Date of 15th birthday",
                                     fiveyear = "Calendar year in 5-year groups",
                                     fouryear = "Calendar year, grouped in 4 years post 2005",
                                     hivstatus_detail = "HIV status detail",
                                     sero_conv_date = "Seroconversion date: midpoint",
                                     start_ep_date = "Date this record starts on, for everyone included in midpoint file",
                                     tested_first_opportunity = "Tested for the first time at the first survey after 15th birthday",
                                     years_one = "Calendar year")
    
    dat$agegrp <- as.numeric(as.character(dat$agegrp))
    labelled::val_labels(dat$agegrp) <- c("0-4"=0,"5-9"=1,"10-14"=2,"15-19"=3,"20-24"=4,"25-29"=5, 
                                          "30-34"=6, "35-39"=7, "40-44"=8, "45-49"=9,  
                                          "50-54"=10, "55-59"=11, "60-64"=12, "65-69"=13,
                                          "70-74"=14, "75-79"=15, "80-84"=16, "85-89"=17,"90+"=18)
    dat$birth_cohort <- as.numeric(as.character(dat$birth_cohort))
    labelled::val_labels(dat$birth_cohort) <- c("1900-1904" = 1, "1905-1909" = 2, "1910-1914" =3,
                                               "1915-1919" = 4, "1920-1924"= 5, "1925-1929" = 6,
                                               "1930-1934" = 7, "1935-1939" = 8, "1940-1944" = 9,
                                               "1945-1949" = 10, "1950-1954" = 11, "1955-1959" = 12, 
                                               "1960-1964" = 13, "1965-1969" = 14, "1970-1974" = 15, 
                                               "1975-1979" = 16, "1980-1984" = 17, "1985-1989" = 18,
                                               "1990-1994" = 19, "1995-1999" = 20, "2000-2004" = 21, 
                                               "2005-2009" = 22, "2010-2014" = 23, "2015-2019" = 24,
                                               "2020-2024" = 25, "2025-2029" = 26)
    dat$fiveyear <- as.numeric(as.character(dat$fiveyear))
    labelled::val_labels(dat$fiveyear) <- c("1989-1994" = 1, "1995-1999" = 2, "2000-2004" = 3,
                                          "2005-2009" = 4, "2010-2014" = 5, "2015-2019" = 6,
                                          "2020-2024"= 7)
    dat$fouryear <- as.numeric(as.character(dat$fouryear))
    labelled::val_labels(dat$fouryear) <- c("earliest-1999"=0,"2000-04"=1,"2004-08"=2,"2009-12"=3,
                                            "2013-16"=4,"2017-20"=5,"2021-24"=6)
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
                                     birth_cohort = "Five year birth_cohort",
                                     fiveyear = "Calendar year in 5-year groups",
                                     fouryear = "Calendar year, grouped in 4 years post 2005")
    
    ############################################################################
    #     SAVE THE DATA
    ############################################################################
    saveRDS(dat, paste0(alphapath,"/ALPHA/incidence_ready_data/",sitename,
                        "/incidence_ready_midpoint_",sitename,".rds"))
    
    