
##################################################
### MAKE HIV test data WIDE FOR MERGING
##################################################

# USES SUMMARY INFORMATION ON RESIDENCY

# changed Nov 2019 to exclude clinic data from AHRI

library(haven);library(tidyverse); library(rlang)

# OPEN DATASET
dat <- haven::read_dta(paste0(alphapath,
                                      "/ALPHA/clean_data/",sitename,
                                      "/hiv_tests_",sitename,".dta"))

# Check that dates are correctly encoded:
if(class(dat$hiv_test_date)!="Date") {
  dat$hiv_test_date <- as.Date(dat$hiv_test_date,origin="1960-01-01")
}
if(class(dat$test_report_date)!="Date") {
  dat$test_report_date <- as.Date(dat$test_report_date,origin="1960-01-01")
}

# merge in the metadata to make info vars and expiry dates, and the 6.1 residency 
# summary for checking test dates
metadat <- haven::read_dta(paste0(alphapath,
                                          "/ALPHA/clean_data/alpha_metadata.dta"))
metadat <- metadat %>% dplyr::select(study_name, study_dx_after,
                                     earlyexit_max, analysis_end_date)
dat <- merge(dat,metadat,by="study_name")
dat <- merge(dat,readRDS(paste0(alphapath,"/ALPHA/prepared_data/",
                                                sitename,"/residency_summary_dates_",
                                                sitename,".rds")),
             by=c("study_name","idno"), all.x=TRUE)

dat <- dat %>% dplyr::select_if(!names(.) %in% "date")
## keep population based results only for this analysis- after discussion with
## AHRI have taken out their clinic results
# (keep self report for Kisumu [study_name 8] only for positives)
dat <- dat %>% dplyr::filter(source_of_test_information==0 | source_of_test_information==1 |
                               (source_of_test_information==4 & study_name==8 & hiv_test_result==1))

dat <- dat %>% dplyr::filter(if_any(matches("sample_type"),
                             ~ !(.==1 & study_name==8)))
dat <- dat %>% dplyr::select_if(!names(.) %in% "sample_type")

# Karonga only - drop tests from before 2000 as problems with these
dat <- dat %>% dplyr::filter(!(study_name==1 & hiv_test_date<as.Date("2000-01-01")))
# Kisumu only - drop self-reported results that are from more than 1 year before the survey
dat <- dat %>% dplyr::filter(!(study_name==8 & source_of_test_information==4 &
                                 (test_report_date - hiv_test_date)>365))
# create a sequence variable and pick up maximum number of tests to use later on in
# looping through test results
dat <- dat %>%   dplyr::group_by(study_name,idno) %>%
  dplyr::arrange(hiv_test_date) %>%
  dplyr::mutate(test_sequence=row_number())
maxtests <- max(dat$test_sequence)

# FIRST & LAST TESTS

## find first and last test dates (all tests)
# any test
dat <- dat %>% dplyr::group_by(study_name,idno) %>%
  dplyr::mutate(firsttest=min(hiv_test_date),
                firsttest_date_all = as.Date(ifelse(is.finite(firsttest), min(firsttest), NA),
                                             origin="1970-01-01"),
                lasttest = max(hiv_test_date),
                lasttest_date_all = as.Date(ifelse(is.finite(lasttest), max(lasttest), NA),
                                            origin="1970-01-01"))
# negative
dat <- dat %>% dplyr::group_by(study_name,idno) %>%
  dplyr::mutate(first_negative = min(hiv_test_date[hiv_test_result==0]),
                first_neg_date_all = as.Date(ifelse(is.finite(first_negative), min(first_negative), NA),
                                             origin="1970-01-01"),
                last_negative = max(hiv_test_date[hiv_test_result==0 & !is.na(hiv_test_date)]),
                last_neg_date_all = as.Date(ifelse(is.finite(last_negative), max(last_negative), NA),
                                            origin="1970-01-01"))
# positive
dat <- dat %>% dplyr::group_by(study_name,idno) %>%
  dplyr::mutate(first_positive = min(hiv_test_date[hiv_test_result==1]),
                first_pos_date_all = as.Date(ifelse(is.finite(first_positive), min(first_positive), NA),
                                             origin="1970-01-01"),
                last_positive = max(hiv_test_date[hiv_test_result==1 & !is.na(hiv_test_date)]),
                last_pos_date_all = as.Date(ifelse(is.finite(last_positive), max(last_positive), NA),
                                            origin="1970-01-01"))
dat <- dat %>% dplyr::select(-c(firsttest,lasttest,first_negative,last_negative,
                                first_positive,last_positive))

# identify tests that are within residency dates

# JUNE 2018 - HAVE REMOVED ALL THIS AS DON'T WANT TO CUT OFF AT analysis_end_date
# (exit_temp has to be created which is the analysis end date as everything after this 
# is dropped later)
# get number of residency episodes
nep <- substr(grep("entry_date",names(dat),value=TRUE),11,13)

dat$test_tag <- NA_integer_
for(i in nep) {
  dat <- dat %>% 
    dplyr::mutate(exit_temp = !!rlang::parse_expr(paste0("exit_date",i)) + earlyexit_max,
                  test_tag = case_when(!is.na(test_tag) ~ test_tag,
                                       is.na(test_tag) & 
                                         hiv_test_date>=!!rlang::parse_expr(paste0("entry_date",i)) &
                                         hiv_test_date<=exit_temp & !is.na(hiv_test_date) & 
                                         !is.na(!!rlang::parse_expr(paste0("entry_date",i))) & 
                                         !is.na(exit_temp) ~ 1L,
                                       TRUE ~ NA_integer_))
}

## find first and last test dates (within the residency dates only)
# any test
dat <- dat %>% dplyr::group_by(study_name,idno) %>%
  dplyr::mutate(testind = ifelse(test_tag==1,hiv_test_date,NA),
                firsttest_date = as.Date(ifelse(is.finite(min(testind, na.rm=TRUE)), 
                                                min(testind, na.rm=TRUE), NA), origin="1970-01-01"),
                lasttest_date = as.Date(ifelse(is.finite(max(testind, na.rm=TRUE)),
                                               max(testind,na.rm=TRUE), NA), origin="1970-01-01"))

# negative
dat <- dat %>% dplyr::group_by(study_name,idno) %>%
  dplyr::mutate(negind = ifelse(test_tag==1 & hiv_test_result==0, hiv_test_date, NA),
                first_neg_date = as.Date(ifelse(is.finite(min(negind, na.rm=TRUE)), 
                                                min(negind,na.rm=TRUE), NA), origin="1970-01-01"),
                last_neg_date = as.Date(ifelse(is.finite(max(negind,na.rm=TRUE)), 
                                               max(negind,na.rm=TRUE), NA), origin="1970-01-01"))

# positive
dat <- dat %>% dplyr::group_by(study_name,idno) %>%
  dplyr::mutate(posind = ifelse(test_tag==1 & hiv_test_result==1, hiv_test_date, NA),
                first_pos_date = as.Date(ifelse(is.finite(min(posind, na.rm=TRUE)),
                                                min(posind, na.rm=TRUE), NA), origin="1970-01-01"),
                last_pos_date = as.Date(ifelse(is.finite(max(posind, na.rm=TRUE)),
                                               max(posind, na.rm=TRUE), NA), origin="1970-01-01"))

dat <- dat %>% dplyr::select(-c(testind, negind,posind,test_tag))


# reshape to wide format, retain all test dates and details
dat <- dat %>%
  tidyr::pivot_wider(names_from=test_sequence,names_sep="",
                     values_from=c(hiv_test_date, hiv_test_result, test_assumption,
                                   test_report_date, source_of_test_information, 
                                   survey_round_name, informed_of_result)) 

# save the dataset
if(!dir.exists(paste0(alphapath,"/ALPHA/prepared_data/",sitename))) {
  dir.create(paste0(alphapath,"/ALPHA/prepared_data/",sitename),recursive=TRUE)
} 
saveRDS(dat, paste0(alphapath,"/ALPHA/prepared_data/",sitename,
                    "/hiv_tests_wide_",sitename,".rds"))

